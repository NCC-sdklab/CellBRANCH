# MVCell
Automatic cell type annotation for single cell RNA data

# Bidirectional Hierarchical Transformer for Cell Type Annotation

A PyTorch implementation of a bidirectional hierarchical transformer model for automated cell type annotation using multiple annotation tools.

## Overview

This model addresses the challenge of cell type annotation by:
- **Bidirectional Processing**: Forward (1→2→3→4→5) and reverse (5→4→3→2→1) prediction paths
- **Hierarchical Constraints**: Enforcing ontological relationships between annotation levels
- **Multi-tool Integration**: Combining predictions from multiple annotation tools (singleR, sctype, scATOMIC, scimilarity, celltypist)
- **Consistency Enforcement**: Ensuring predictions maintain biological hierarchy

## Architecture

```
Input: [Tool Annotations + Confidence Scores]
    ↓
Feature Projection
    ↓
┌─────────────────┬─────────────────┐
│  Forward Path   │  Reverse Path   │
│   (1→2→3→4→5)   │   (5→4→3→2→1)   │
│  + Hierarchy    │  + Cross-Attn   │
│  Constraints    │  with Forward   │
└─────────────────┴─────────────────┘
    ↓
Combined Loss: Forward(3.0x) + Reverse(0.3x) + Consistency(0.1x)
```

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/bidirectional-hierarchical-transformer.git
cd bidirectional-hierarchical-transformer

# Install dependencies
pip install -r requirements.txt

# Install package
pip install -e .
```

## Quick Start

### 1. Data Preparation

Ensure your data files are organized as follows:
```
data/
├── merged_annlevel_with_new_fin.tsv          # Training data
├── merged_annlevel_with_new_test_fin.tsv     # Test data
└── hlca_celltype_tree/
    └── unique_annotations_ontology_JHL_with_SW.csv  # Ontology file
```

**Expected data format:**
- Training/Test files should contain:
  - Tool annotation columns: `{tool}_ann{level}` (e.g., `singleR_ann1`, `sctype_ann2`)
  - Confidence score columns: `{tool}_score` (e.g., `singleR_score`, `sctype_score`)
  - Ground truth columns: `ann_level_{level}` (e.g., `ann_level_1`, `ann_level_2`)

### 2. Training

```bash
# Basic training
python train.py --data-dir ./data --output-dir ./results

# Training with custom parameters
python train.py \
    --data-dir ./data \
    --output-dir ./results \
    --epochs 100 \
    --batch-size 512 \
    --lr 1e-6 \
    --d-model 512 \
    --nhead 16 \
    --num-layers 4

# Training with ablation (removing specific tools)
python train.py \
    --data-dir ./data \
    --output-dir ./results \
    --ablation sctype scATOMIC
```

### 3. Prediction and Evaluation

```bash
# Test set evaluation
python predict.py \
    --model-dir ./results/bidirectional_transformer_seed_9401 \
    --data-dir ./data \
    --mode test

# Validation set evaluation
python predict.py \
    --model-dir ./results/bidirectional_transformer_seed_9401 \
    --data-dir ./data \
    --mode val
```

## Command Line Arguments

### Training Script (`train.py`)

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--data-dir` | str | `.` | Data directory path |
| `--output-dir` | str | `./results` | Output directory for models and results |
| `--seed` | int | `9401` | Random seed for reproducibility |
| `--epochs` | int | `80` | Number of training epochs |
| `--batch-size` | int | `256` | Training batch size |
| `--lr` | float | `1e-5` | Learning rate |
| `--d-model` | int | `256` | Model dimension |
| `--nhead` | int | `8` | Number of attention heads |
| `--num-layers` | int | `2` | Number of transformer layers |
| `--dropout` | float | `0.5` | Dropout rate |
| `--weight-decay` | float | `5e-3` | Weight decay for optimizer |
| `--patience` | int | `15` | Early stopping patience |
| `--ablation` | str* | None | Tools to ablate (space-separated) |

### Prediction Script (`predict.py`)

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--model-dir` | str | Required | Directory containing trained model |
| `--data-dir` | str | `.` | Data directory path |
| `--mode` | str | `test` | Evaluation mode: 'test' or 'val' |
| `--batch-size` | int | `1024` | Prediction batch size |
| `--output-dir` | str | None | Output directory (defaults to model-dir) |

## Programmatic Usage

```python
from utils import set_seed, load_encoders
from model import SimpleBidirectionalTransformer
import torch

# Set random seed
set_seed(9401)

# Load trained model and encoders
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
enc_in, enc_out = load_encoders("./results/model_directory")

# Initialize model (parameters should match training)
model = SimpleBidirectionalTransformer(
    enc_in=enc_in,
    enc_out=enc_out,
    ont_df=ont_df,  # Ontology dataframe
    score_dim=5,    # Number of tools
    terminal_nodes_idx=terminal_nodes_idx,
    d_model=256,
    nhead=8,
    num_layers=2,
    device=device
).to(device)

# Load trained weights
model.load_state_dict(torch.load("./results/best_model.pth"))
model.eval()

# Make predictions
with torch.no_grad():
    forward_outputs, reverse_outputs = model(cat_features, score_features, enc_in=enc_in)
    predictions = {lvl: outputs.argmax(-1) for lvl, outputs in forward_outputs.items()}
```

## Output Files

### Training Output
- `best_multilevel_seed{SEED}.pth`: Best model weights
- `encoders_in.pkl`: Input label encoders
- `encoders_out.pkl`: Output label encoders
- `bidirectional_transformer_curves.png`: Training curves visualization
- `loss_curve.png`: Loss curve plot

### Prediction Output
- `comparison_metrics_for_{mode}.csv`: Detailed metrics comparison
- `comparison_ann_level_{i}_for_{mode}.png`: Performance comparison plots

## Model Architecture Details

### Bidirectional Design
1. **Forward Path (1→2→3→4→5)**: 
   - Primary prediction path with hierarchical constraints
   - Ensures biological validity through ontology-based masking
   - Higher weight in loss function (3.0x)

2. **Reverse Path (5→4→3→2→1)**:
   - Auxiliary prediction path for consistency
   - Uses cross-attention with forward features
   - Lower weight in loss function (0.3x)

### Loss Function Components
- **Forward Loss**: Weighted cross-entropy with class balancing
- **Reverse Loss**: Standard cross-entropy for auxiliary predictions
- **Consistency Loss**: KL divergence between forward and reverse predictions

### Hierarchical Constraints
- Enforced only in forward path during training and inference
- Based on parent-child relationships in cell type ontology
- Prevents biologically impossible predictions

## Performance Metrics

The model is evaluated using:
- **Accuracy**: Overall prediction accuracy
- **Precision/Recall/F1**: Macro-averaged across all classes
- **Matthews Correlation Coefficient (MCC)**: Balanced metric for imbalanced classes
- **Hierarchy Consistency**: Percentage of predictions following ontological rules

## Troubleshooting

### Common Issues

1. **CUDA Out of Memory**:
   - Reduce `--batch-size`
   - Use gradient accumulation
   - Enable mixed precision training

2. **Unknown Categories During Prediction**:
   - Ensure test data categories are present in training data
   - Check encoder compatibility

3. **Poor Convergence**:
   - Adjust learning rate (`--lr`)
   - Increase model capacity (`--d-model`, `--num-layers`)
   - Check data quality and preprocessing

### Memory Optimization

For large datasets:
```python
# Enable gradient checkpointing
model.gradient_checkpointing_enable()

# Use mixed precision
from torch.cuda.amp import autocast, GradScaler
scaler = GradScaler()
```

## Citation

If you use this code in your research, please cite:

```bibtex
@article{bidirectional_transformer_2024,
  title={Bidirectional Hierarchical Transformer for Cell Type Annotation},
  author={Your Name},
  journal={Journal Name},
  year={2024}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Acknowledgments

- Thanks to the developers of the cell type annotation tools: singleR, sctype, scATOMIC, scimilarity, and celltypist
- Human Lung Cell Atlas (HLCA) for providing the ontological framework
