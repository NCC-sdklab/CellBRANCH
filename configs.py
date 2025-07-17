# config.py
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
        self.forward_weight = 1.5
        self.reverse_weight = 0.8
        self.consistency_weight = 0.1
        
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

---

# setup.py
"""
Setup script for Bidirectional Hierarchical Transformer package
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="bidirectional-hierarchical-transformer",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Bidirectional Hierarchical Transformer for Cell Type Annotation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/bidirectional-hierarchical-transformer",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "bht-train=train:main",
            "bht-predict=predict:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.md", "*.txt"],
    },
)

---

# requirements.txt
torch>=2.0.0
numpy>=1.21.0
pandas>=1.3.0
scikit-learn>=1.0.0
matplotlib>=3.5.0
tqdm>=4.62.0