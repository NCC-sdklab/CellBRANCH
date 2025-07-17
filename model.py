import torch
import torch.nn as nn
import torch.nn.functional as F
import math
from torch.utils.data import Dataset
from typing import Dict

class MultiLevelDataset(Dataset):
    """Dataset for multi-level cell type annotation"""
    def __init__(self, X_df, y_arr, ann_cols, score_cols):
        self.cat = torch.tensor(X_df[[f"{c}_idx" for c in ann_cols]].values, dtype=torch.long)
        self.scores = torch.tensor(X_df[score_cols].values, dtype=torch.float32)
        self.y = torch.tensor(y_arr, dtype=torch.long)
    
    def __len__(self):
        return len(self.y)
    
    def __getitem__(self, idx):
        return self.cat[idx], self.scores[idx], self.y[idx]

class PositionalEncoding(nn.Module):
    """Positional encoding for transformer"""
    def __init__(self, d_model, dropout=0.1, max_len=10):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)
        
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)
        
    def forward(self, x):
        seq_len = x.size(1)
        x = x + self.pe[:seq_len, :].unsqueeze(0)
        return self.dropout(x)

class SimpleBidirectionalTransformer(nn.Module):
    """Bidirectional Hierarchical Transformer for cell type annotation"""
    
    def __init__(self, enc_in, enc_out, ont_df, score_dim, terminal_nodes_idx, 
                 d_model=128, nhead=4, num_layers=2, dropout=0.1, device='cpu'):
        super().__init__()
        self.d_model = d_model
        self.terminal_nodes_idx = terminal_nodes_idx
        self.device = device
        
        # Input feature projection
        total_emb_dim = sum(len(le.classes_) for le in enc_in.values())
        self.input_projection = nn.Sequential(
            nn.Linear(total_emb_dim + score_dim, d_model),
            nn.LayerNorm(d_model),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # Phase 1: Reverse direction embeddings (5→4→3→2→1)
        self.reverse_level_embeddings = nn.ModuleList([
            nn.Embedding(len(enc_out[f'ann_level_{i}'].classes_), d_model) for i in range(5, 0, -1)
        ])
        
        # Phase 2: Forward direction embeddings (1→2→3→4→5)
        self.forward_level_embeddings = nn.ModuleList([
            nn.Embedding(len(enc_out[f'ann_level_{i}'].classes_), d_model) for i in range(1, 6)
        ])
        
        # Transformer layers
        self.reverse_transformer = nn.TransformerDecoder(
            nn.TransformerDecoderLayer(d_model, nhead, d_model*2, dropout=dropout, batch_first=True),
            num_layers=num_layers
        )
        
        self.forward_transformer = nn.TransformerDecoder(
            nn.TransformerDecoderLayer(d_model, nhead, d_model*2, dropout=dropout, batch_first=True),
            num_layers=num_layers
        )
        
        # Cross-attention
        self.cross_attention = nn.MultiheadAttention(d_model, nhead//2, dropout=dropout, batch_first=True)
        self.cross_norm = nn.LayerNorm(d_model)
        
        # Output heads
        self.reverse_output_heads = nn.ModuleList([
            nn.Sequential(
                nn.Linear(d_model, len(enc_out[f'ann_level_{i}'].classes_)),
                nn.Dropout(dropout)
            ) for i in range(5, 0, -1)
        ])
        
        self.forward_output_heads = nn.ModuleList([
            nn.Sequential(
                nn.Linear(d_model, len(enc_out[f'ann_level_{i}'].classes_)),
                nn.Dropout(dropout)
            ) for i in range(1, 6)
        ])
        
        # Register unknown indices
        for i in range(1, 6):
            unk_idx = enc_out[f'ann_level_{i}'].transform(["Unknown"])[0]
            self.register_buffer(f"unknown_{i}", torch.tensor(unk_idx, device=device))
        
        # Create hierarchical masks
        self._create_hierarchical_masks(ont_df, enc_out)
    
    def _create_hierarchical_masks(self, ont_df, enc_out):
        """Create hierarchical constraint masks"""
        for i in range(1, 5):
            p, c = f"ann_level_{i}", f"ann_level_{i+1}"
            rel = ont_df[[p, c]].fillna("Unknown").groupby(p)[c].unique().to_dict()
            p_size, c_size = len(enc_out[p].classes_), len(enc_out[c].classes_)
            mask = torch.zeros(p_size, c_size, dtype=torch.bool, device=self.device)
            
            for p_cls, children in rel.items():
                try:
                    p_idx = enc_out[p].transform([p_cls])[0]
                    for ch in children:
                        try:
                            c_idx = enc_out[c].transform([ch])[0]
                            mask[p_idx, c_idx] = True
                        except:
                            pass
                except:
                    pass
            
            self.register_buffer(f"mask_{i}_{i+1}", mask)
    
    def forward(self, cat_idxs, scores, parent_idxs=None, enc_in=None):
        batch_size = cat_idxs.size(0)
        
        # Input encoding
        cat_features = []
        for i, col in enumerate(enc_in.keys()):
            one_hot = F.one_hot(cat_idxs[:, i], num_classes=len(enc_in[col].classes_)).float()
            cat_features.append(one_hot)
        
        input_features = torch.cat(cat_features + [scores], dim=-1)
        memory = self.input_projection(input_features).unsqueeze(1)
        
        # Phase 1: Forward prediction (1→2→3→4→5)
        forward_outputs = {}
        forward_features = []
        
        for level_idx, level in enumerate([1, 2, 3, 4, 5]):
            if level_idx == 0:
                tgt = torch.zeros(batch_size, 1, self.d_model, device=cat_idxs.device)
            else:
                prev_level = level - 1
                if parent_idxs is not None and self.training:
                    prev_pred = parent_idxs[:, prev_level-1]
                else:
                    prev_pred = forward_outputs[f'ann_level_{prev_level}'].argmax(-1)
                
                prev_emb = self.forward_level_embeddings[level_idx-1](prev_pred)
                tgt = prev_emb.unsqueeze(1)
            
            # Forward transformer
            output = self.forward_transformer(tgt, memory)
            forward_features.append(output)
            
            # Output projection
            logits = self.forward_output_heads[level_idx](output.squeeze(1))
            
            # Apply hierarchical constraints
            if level > 1:
                prev_level = level - 1
                if parent_idxs is not None and self.training:
                    prev_pred = parent_idxs[:, prev_level-1]
                else:
                    prev_pred = forward_outputs[f'ann_level_{prev_level}'].argmax(-1)
                
                mask = getattr(self, f"mask_{prev_level}_{level}")[prev_pred]
                neg_inf = torch.full_like(logits, -1e4)
                logits = neg_inf.masked_scatter(mask, logits.masked_select(mask))
            
            forward_outputs[f'ann_level_{level}'] = logits
        
        # Phase 2: Reverse prediction (5→4→3→2→1)
        reverse_outputs = {}
        forward_context = torch.cat(forward_features, dim=1)
        
        for level_idx, level in enumerate([5, 4, 3, 2, 1]):
            if level_idx == 0:
                tgt = torch.zeros(batch_size, 1, self.d_model, device=cat_idxs.device)
            else:
                prev_level = [5, 4, 3, 2, 1][level_idx-1]
                if parent_idxs is not None and self.training:
                    prev_pred = parent_idxs[:, prev_level-1]
                else:
                    prev_pred = reverse_outputs[f'reverse_ann_level_{prev_level}'].argmax(-1)
                
                prev_emb = self.reverse_level_embeddings[level_idx-1](prev_pred)
                tgt = prev_emb.unsqueeze(1)
            
            # Cross-attention with Forward context
            try:
                attended_tgt, _ = self.cross_attention(tgt, forward_context, forward_context)
                tgt = self.cross_norm(tgt + attended_tgt * 0.2)
            except:
                pass
            
            # Reverse transformer
            output = self.reverse_transformer(tgt, memory)
            
            # Output projection
            logits = self.reverse_output_heads[level_idx](output.squeeze(1))
            reverse_outputs[f'reverse_ann_level_{level}'] = logits
        
        return forward_outputs, reverse_outputs

def stable_bidirectional_loss(forward_outputs, reverse_outputs, targets, terminal_nodes_idx, model,
                              target_cols, forward_weight=3.0, reverse_weight=0.3, 
                              consistency_weight=0.1, device='cpu'):
    """Bidirectional loss function"""
    
    total_loss = 0.0
    level_weights = [1.0, 1.2, 1.5, 2.0, 2.5]
    
    # 1. Forward direction loss (MAIN LOSS)
    forward_loss = 0.0
    for idx, lvl in enumerate(target_cols):
        unk = getattr(model, f"unknown_{idx+1}")
        valid = targets[:, idx] != unk
        
        if valid.any():
            valid_targets = targets[valid, idx]
            valid_outputs = forward_outputs[lvl][valid]
            
            # Class weights for imbalanced classes
            unique_classes, class_counts = torch.unique(valid_targets, return_counts=True)
            total_samples = len(valid_targets)
            
            class_weights = torch.ones(forward_outputs[lvl].size(1), device=device)
            for class_idx, count in zip(unique_classes, class_counts):
                freq = count.float() / total_samples
                weight = torch.clamp(1.0 / torch.sqrt(freq + 1e-4), min=1.0, max=8.0)
                class_weights[class_idx] = weight
            
            # Weighted cross entropy
            ce_loss = F.cross_entropy(valid_outputs, valid_targets, 
                                    weight=class_weights, reduction='mean')
            
            # NaN check
            if torch.isnan(ce_loss) or torch.isinf(ce_loss):
                ce_loss = F.cross_entropy(valid_outputs, valid_targets, reduction='mean')
            
            forward_loss += level_weights[idx] * ce_loss
    
    total_loss += forward_weight * forward_loss
    
    # 2. Reverse direction loss (AUXILIARY LOSS)
    reverse_loss = 0.0
    for idx, level in enumerate([5, 4, 3, 2, 1]):
        reverse_key = f'reverse_ann_level_{level}'
        if reverse_key in reverse_outputs:
            unk = getattr(model, f"unknown_{level}")
            valid = targets[:, level-1] != unk
            
            if valid.any():
                valid_targets = targets[valid, level-1]
                valid_outputs = reverse_outputs[reverse_key][valid]
                
                ce_loss = F.cross_entropy(valid_outputs, valid_targets, reduction='mean')
                
                if not (torch.isnan(ce_loss) or torch.isinf(ce_loss)):
                    reverse_loss += ce_loss
    
    total_loss += reverse_weight * reverse_loss
    
    # 3. Consistency loss
    consistency_loss = 0.0
    for level in range(1, 6):
        forward_key = f'ann_level_{level}'
        reverse_key = f'reverse_ann_level_{level}'
        
        if forward_key in forward_outputs and reverse_key in reverse_outputs:
            forward_probs = F.softmax(forward_outputs[forward_key], dim=-1)
            
            kl_loss = F.kl_div(
                F.log_softmax(reverse_outputs[reverse_key], dim=-1),
                forward_probs.detach(),
                reduction='batchmean'
            )
            
            if not (torch.isnan(kl_loss) or torch.isinf(kl_loss)):
                consistency_loss += kl_loss
    
    total_loss += consistency_weight * consistency_loss
    
    # Final NaN/inf check
    if torch.isnan(total_loss) or torch.isinf(total_loss):
        print("Warning: Loss is NaN/inf, using simple forward loss")
        simple_loss = 0.0
        for idx, lvl in enumerate(target_cols):
            unk = getattr(model, f"unknown_{idx+1}")
            valid = targets[:, idx] != unk
            if valid.any():
                simple_loss += F.cross_entropy(forward_outputs[lvl][valid], targets[valid, idx])
        return simple_loss / len(target_cols)
    
    return total_loss