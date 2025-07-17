#!/usr/bin/env python3
"""
Bidirectional Hierarchical Transformer Training Script
"""

import os
import argparse
import numpy as np
import pandas as pd
import torch
import torch.optim as optim
from torch.utils.data import DataLoader
from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import contextlib
from matplotlib.ticker import MaxNLocator

# Import custom modules
from utils import (
    set_seed, get_all_possible_values, identify_terminal_nodes,
    create_label_encoders, save_encoders, get_terminal_nodes_idx,
    preprocess_data, verify_hierarchy_consistency
)
from model import (
    MultiLevelDataset, SimpleBidirectionalTransformer, 
    stable_bidirectional_loss
)

class TrainingHistory:
    """Class to track training history"""
    def __init__(self, target_cols):
        self.train_loss = []
        self.val_loss = []
        self.consistency_rate = []
        self.zero_precision_cells = []
        self.stagnant_epochs = []
        self.level_f1_scores = {lvl: [] for lvl in target_cols}
        self.level_accuracy = {lvl: [] for lvl in target_cols}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Train Bidirectional Transformer')
    parser.add_argument('--data-dir', type=str, default='.', help='Data directory')
    parser.add_argument('--output-dir', type=str, default='./results', help='Output directory')
    parser.add_argument('--seed', type=int, default=9401, help='Random seed')
    parser.add_argument('--epochs', type=int, default=80, help='Number of epochs')
    parser.add_argument('--batch-size', type=int, default=256, help='Batch size')
    parser.add_argument('--lr', type=float, default=1e-5, help='Learning rate')
    parser.add_argument('--d-model', type=int, default=256, help='Model dimension')
    parser.add_argument('--nhead', type=int, default=8, help='Number of attention heads')
    parser.add_argument('--num-layers', type=int, default=2, help='Number of transformer layers')
    parser.add_argument('--dropout', type=float, default=0.5, help='Dropout rate')
    parser.add_argument('--weight-decay', type=float, default=5e-3, help='Weight decay')
    parser.add_argument('--patience', type=int, default=15, help='Early stopping patience')
    parser.add_argument('--ablation', type=str, nargs='*', default=None, 
                       help='Tools to ablate (e.g., sctype scATOMIC)')
    
    return parser.parse_args()

def load_data(data_dir):
    """Load training and test data"""
    df_train = pd.read_csv(os.path.join(data_dir, "merged_annlevel_with_new_fin.tsv"), sep='\t')
    df_test = pd.read_csv(os.path.join(data_dir, "merged_annlevel_with_new_test_fin.tsv"), sep='\t')
    ont = pd.read_csv(os.path.join(data_dir, "hlca_celltype_tree/unique_annotations_ontology_JHL_with_SW.csv"))
    
    return df_train, df_test, ont

def run_epoch(model, loader, optimizer=None, is_train=False, device='cpu', 
              enc_in=None, terminal_nodes_idx=None, target_cols=None):
    """Run one epoch of training or validation"""
    model.train() if is_train else model.eval()
    total_loss = 0.0
    all_preds = {lvl: [] for lvl in target_cols}
    all_trues = {lvl: [] for lvl in target_cols}
    
    with torch.no_grad() if not is_train else contextlib.nullcontext():
        for cat, scores, targets in tqdm(loader, desc=("Train" if is_train else "Val")):
            cat, scores, targets = cat.to(device), scores.to(device), targets.to(device)
            
            with torch.amp.autocast('cuda'):
                forward_outputs, reverse_outputs = model(cat, scores, 
                                                        parent_idxs=targets if is_train else None,
                                                        enc_in=enc_in)
                
                loss = stable_bidirectional_loss(forward_outputs, reverse_outputs, targets, 
                                                terminal_nodes_idx, model, target_cols,
                                                forward_weight=1.5, reverse_weight=0.8, 
                                                consistency_weight=0.1, device=device)
                
                # Collect predictions
                for idx, lvl in enumerate(target_cols):
                    preds = forward_outputs[lvl].argmax(1).cpu().tolist()
                    all_preds[lvl].extend(preds)
                    all_trues[lvl].extend(targets[:, idx].cpu().tolist())
            
            if is_train:
                optimizer.zero_grad(set_to_none=True)
                if device.type == 'cuda':
                    scaler = torch.amp.GradScaler('cuda')
                    scaler.scale(loss).backward()
                    scaler.unscale_(optimizer)
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=0.05)
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=0.05)
                    optimizer.step()
            
            total_loss += float(loss)
    
    return total_loss / len(loader), all_preds, all_trues

def calculate_level_metrics(val_preds, val_trues, enc_out, target_cols):
    """Calculate metrics for each level"""
    level_f1 = {}
    level_acc = {}
    zero_precision_count = 0

    for lvl in target_cols:
        unk = enc_out[lvl].transform(["Unknown"])[0]
        y_true = [t for t in val_trues[lvl] if t != unk]
        y_pred = [p for t, p in zip(val_trues[lvl], val_preds[lvl]) if t != unk]
        
        if not y_true:
            level_f1[lvl] = 0.0
            level_acc[lvl] = 0.0
            continue
            
        try:
            rpt = classification_report(
                y_true, y_pred,
                labels=[i for i, c in enumerate(enc_out[lvl].classes_) if c != "Unknown"],
                target_names=[c for c in enc_out[lvl].classes_ if c != "Unknown"],
                zero_division=0,
                output_dict=True
            )
            
            macro_f1 = rpt.get('macro avg', {}).get('f1-score', 0.0)
            level_f1[lvl] = macro_f1
            
            if 'accuracy' in rpt:
                acc = rpt['accuracy']
            else:
                acc = sum(1 for t, p in zip(y_true, y_pred) if t == p) / len(y_true)
            level_acc[lvl] = acc
            
            # Count zero precision classes
            for class_name, metrics in rpt.items():
                if isinstance(metrics, dict) and 'precision' in metrics:
                    if metrics['precision'] == 0.0:
                        zero_precision_count += 1
                        
        except Exception as e:
            level_f1[lvl] = 0.0
            level_acc[lvl] = 0.0
    
    return level_f1, level_acc, zero_precision_count

def save_training_plots(training_history, res_dir, best_epoch=None):
    """Save training visualization plots"""
    if best_epoch is None:
        best_epoch = len(training_history.train_loss)
    
    plt.figure(figsize=(20, 15))
    
    # Training/Validation Loss
    plt.subplot(3, 3, 1)
    epochs = range(1, len(training_history.train_loss) + 1)
    plt.plot(epochs, training_history.train_loss, 'b-', label='Train Loss')
    plt.plot(epochs, training_history.val_loss, 'r-', label='Val Loss')
    if best_epoch:
        plt.axvline(x=best_epoch, color='g', linestyle='--', label=f'Best Model (Epoch {best_epoch})')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Training & Validation Loss (Bidirectional Transformer)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    
    # Hierarchy Consistency
    plt.subplot(3, 3, 2)
    plt.plot(epochs, training_history.consistency_rate, 'g-', label='Hierarchy Consistency')
    if best_epoch:
        plt.axvline(x=best_epoch, color='b', linestyle='--', label=f'Best Model (Epoch {best_epoch})')
    plt.xlabel('Epoch')
    plt.ylabel('Consistency Rate')
    plt.title('Hierarchy Consistency Rate')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(0.5, 1.0)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    
    # Additional plots can be added here...
    
    plt.tight_layout()
    plt.savefig(os.path.join(res_dir, 'bidirectional_transformer_curves.png'), dpi=200)
    plt.close()

def main():
    args = parse_args()
    
    # Setup
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    SEED = set_seed(args.seed)
    
    # Constants
    TOOLS = ["singleR", "sctype", "scATOMIC", "scimilarity", "celltypist"]
    SCORE_COLS = [f"{t}_score" for t in TOOLS]
    ANN_COLS = [f"{t}_ann{i}" for t in TOOLS for i in range(1, 6)]
    ANN_LEVEL = [f"ann_level_{i}" for i in range(1, 6)]
    TARGET_COLS = [f"ann_level_{i}" for i in range(1, 6)]
    
    # Create results directory
    ablation_str = f"_ablation_{'_'.join(args.ablation)}" if args.ablation else ""
    res_dir = os.path.join(args.output_dir, f"bidirectional_transformer{ablation_str}_seed_{SEED}")
    os.makedirs(res_dir, exist_ok=True)
    
    print(f"Results will be saved to: {res_dir}")
    
    # Load data
    print("Loading data...")
    df_train, df_test, ont = load_data(args.data_dir)
    print(f"Training data: {df_train.shape}")
    print(f"Test data: {df_test.shape}")
    
    # Preprocess data
    print("Preprocessing data...")
    df_train = preprocess_data(df_train, SCORE_COLS, ANN_COLS, ANN_LEVEL)
    df_test = preprocess_data(df_test, SCORE_COLS, ANN_COLS, ANN_LEVEL)
    
    # Apply ablation if specified
    if args.ablation:
        print(f"Applying ablation to: {args.ablation}")
        from utils import ablation_data
        df_train = ablation_data(df_train, args.ablation)
        df_test = ablation_data(df_test, args.ablation)
    
    # Get all possible values and create encoders
    print("Creating label encoders...")
    all_cell_types = get_all_possible_values(ont, df_train, df_test)
    enc_in, enc_out = create_label_encoders(df_train, df_test, all_cell_types)
    
    # Save encoders
    save_encoders(enc_in, enc_out, res_dir)
    
    # Get terminal nodes
    terminal_nodes = identify_terminal_nodes(ont)
    terminal_nodes_idx = get_terminal_nodes_idx(terminal_nodes, enc_out)
    
    # Prepare data
    feature_cols = [f"{c}_idx" for c in ANN_COLS] + SCORE_COLS
    target_cols = [f"{lvl}_idx" for lvl in TARGET_COLS]
    
    X_train_full = df_train[feature_cols]
    y_train_full = df_train[target_cols].values
    train_X, val_X, train_y, val_y = train_test_split(
        X_train_full, y_train_full, test_size=0.4, 
        stratify=y_train_full[:, 0], random_state=SEED
    )
    
    test_X = df_test[feature_cols]
    test_y = df_test[target_cols].values
    
    # Create datasets
    train_ds = MultiLevelDataset(train_X, train_y, ANN_COLS, SCORE_COLS)
    val_ds = MultiLevelDataset(val_X, val_y, ANN_COLS, SCORE_COLS)
    test_ds = MultiLevelDataset(test_X, test_y, ANN_COLS, SCORE_COLS)
    
    # Create data loaders
    dt_args = dict(num_workers=4, pin_memory=True if device.type == 'cuda' else False, 
                   persistent_workers=True, drop_last=True)
    data_loaders = {
        'train': DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, **dt_args),
        'val': DataLoader(val_ds, batch_size=args.batch_size * 2, shuffle=False, **dt_args),
        'test': DataLoader(test_ds, batch_size=args.batch_size * 2, shuffle=False, **dt_args)
    }
    
    # Initialize model
    print("Initializing model...")
    model = SimpleBidirectionalTransformer(
        enc_in, enc_out, ont, len(SCORE_COLS), terminal_nodes_idx,
        d_model=args.d_model, nhead=args.nhead, num_layers=args.num_layers,
        dropout=args.dropout, device=device
    ).to(device)
    
    if hasattr(torch, 'compile'):
        model = torch.compile(model, dynamic=True, mode='reduce-overhead')
    
    optimizer = optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    
    # Training
    print("Starting training...")
    training_history = TrainingHistory(TARGET_COLS)
    best_val, patience = float('inf'), 0
    best_epoch_num = 0
    stagnant_count = 0
    
    for ep in range(1, args.epochs + 1):
        tr_loss, _, _ = run_epoch(model, data_loaders['train'], optimizer, 
                                is_train=True, device=device, enc_in=enc_in,
                                terminal_nodes_idx=terminal_nodes_idx, target_cols=TARGET_COLS)
        
        vl_loss, val_preds, val_trues = run_epoch(model, data_loaders['val'], 
                                                is_train=False, device=device, enc_in=enc_in,
                                                terminal_nodes_idx=terminal_nodes_idx, 
                                                target_cols=TARGET_COLS)
        
        training_history.train_loss.append(tr_loss)
        training_history.val_loss.append(vl_loss)
        
        print(f"Epoch {ep}: train={tr_loss:.4f}, val={vl_loss:.4f}")
        
        # Calculate level metrics
        level_f1, level_acc, zero_precision_count = calculate_level_metrics(
            val_preds, val_trues, enc_out, TARGET_COLS
        )
        
        for lvl in TARGET_COLS:
            training_history.level_f1_scores[lvl].append(level_f1[lvl])
            training_history.level_accuracy[lvl].append(level_acc[lvl])
        
        training_history.zero_precision_cells.append(zero_precision_count)
        
        # Validation metrics per level
        print("Validation Metrics per Level:")
        for lvl in TARGET_COLS:
            print(f"  {lvl} - acc: {level_acc[lvl]:.2f}, macro F1: {level_f1[lvl]:.2f}")
        
        # Early stopping
        if vl_loss < best_val:
            best_val, patience = vl_loss, 0
            best_epoch_num = ep
            stagnant_count = 0
            torch.save(model.state_dict(), os.path.join(res_dir, f"best_multilevel_seed{SEED}.pth"))
            print(f"New best model saved with val loss: {vl_loss:.4f}")
        else:
            patience += 1
            stagnant_count += 1
            if patience >= args.patience:
                print("Early stopping.")
                break
        
        training_history.stagnant_epochs.append(stagnant_count)
    
    # Save training plots
    save_training_plots(training_history, res_dir, best_epoch_num)
    
    # Test evaluation
    print("Evaluating on test set...")
    model.load_state_dict(torch.load(os.path.join(res_dir, f"best_multilevel_seed{SEED}.pth")))
    test_loss, test_preds, test_trues = run_epoch(model, data_loaders['test'], 
                                                is_train=False, device=device, enc_in=enc_in,
                                                terminal_nodes_idx=terminal_nodes_idx, 
                                                target_cols=TARGET_COLS)
    print(f"Test loss: {test_loss:.4f}")
    
    # Save test results
    test_f1, test_acc, _ = calculate_level_metrics(test_preds, test_trues, enc_out, TARGET_COLS)
    
    print("\n=== Final Results ===")
    for lvl in TARGET_COLS:
        print(f"{lvl} - Test Acc: {test_acc[lvl]:.3f}, Test F1: {test_f1[lvl]:.3f}")
    
    print(f"\nTraining completed. Results saved to {res_dir}")

if __name__ == "__main__":
    main()