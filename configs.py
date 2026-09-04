# configs.py
"""
Configuration file for Bidirectional Hierarchical Transformer
"""

import os
import torch

# Data configuration
TOOLS = ["singleR", "sctype", "scATOMIC", "scimilarity", "celltypist"]
SCORE_COLS = [f"{t}_score" for t in TOOLS]
ANN_COLS = [f"{t}_ann{i}" for t in TOOLS for i in range(1, 6)]
ANN_LEVEL = [f"ann_level_{i}" for i in range(1, 6)]
TARGET_COLS = [f"ann_level_{i}" for i in range(1, 6)]

# Model hyperparameters
class ModelConfig:
    def __init__(self):
        # Model architecture
        self.d_model = 256
        self.nhead = 8
        self.num_layers = 2
        self.dropout = 0.4
        
        # Training parameters
        self.lr = 1e-6
        self.weight_decay = 1e-3
        self.batch_size = 1024
        self.epochs = 80
        self.patience = 15
        
        # Loss function weights
        self.forward_weight = 2.0
        self.reverse_weight = 0.2
        self.consistency_weight = 0.15
        
        # Data split
        self.val_size = 0.4
        
        # Device settings
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.num_workers = 4

# File paths
class PathConfig:
    def __init__(self, data_dir="."):
        self.data_dir = data_dir
        self.train_file = os.path.join(data_dir, "merged_annlevel_with_new_fin.tsv")
        self.test_file = os.path.join(data_dir, "merged_annlevel_with_new_test_fin.tsv")
        self.ont_file = os.path.join(data_dir, "hlca_celltype_tree/unique_annotations_ontology_JHL_with_SW.csv")

# Default configuration
def get_default_config():
    return {
        'model': ModelConfig(),
        'paths': PathConfig(),
        'tools': TOOLS,
        'score_cols': SCORE_COLS,
        'ann_cols': ANN_COLS,
        'target_cols': TARGET_COLS
    }
