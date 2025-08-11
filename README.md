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

### 1. Clone the repository
```bash
git clone https://github.com/NCC-sdklab/MVCell.git
cd MVCell
```

---

### 2. Create the required conda environments

MVCell uses **three** conda environments:

- **`MVCell_R`** → R environment to run R-based celltype annotations (SingleR, sctype, scATOMIC) 
- **`r-reticulate-mvcell`** → Python environment for R’s `reticulate` interface for scATOMIC
- **`MVCell_Python`** → Main environment for MVCell + celltypist + scimilarity


Create them from the provided YAML files:

```bash
# Create MVCell_base
cd r-isntall
mamba env create -f r-install/MVCell_R_env.yaml   # or: conda env create -f ...

# Create r-reticulate
mamba env create -f r-reticulate-mvcell_env.yaml
```

> **Why two environments?**  
> - Keeping `r-reticulate-mvcell` separate avoids Python version/package conflicts with the main MVCell environment.  
> - This separation also keeps `r-reticulate` alone for other R projects.

---

### 3. Install remaining R packages

Activate the **`MVCell_R`** environment and run the R install script:

```bash
#cd r-install
conda activate MVCell_R
Rscript MVCell_R_install.R
```

### 4. Install MVCell python (To be completed)

```
### Install dependencies
```
pip install -r requirements.txt
```

### Install package
```
pip install -e .
```

# Quick Start

## Usage Instructions for Cell Type Annotation in MVCell

Although these scripts are not strictly part of MVCell's core code, MVCell requires users to run
**ScType**, **SingleR**, and **scATOMIC** in R to prepare annotation outputs.
To simplify this process, two R scripts and dedicated conda environments are provided.

---

### **Requirements**

- A `.h5ad` file containing the **count matrix**.
- The provided conda environments:
  - `MVCell_base`
  - `r-reticulate`
  - Any additional environments as specified in the installation instructions.

---

### **Step 1 — Convert `.h5ad` to `SingleCellExperiment` (SCE)**

Run the `MVCell_V1_R_init.R` script to:
- Load the `.h5ad` file
- Convert it to an `SingleCellExperiment` (SCE) object
- Generate log-normalized counts (`logcounts`)

```bash
conda activate MVCell_R
Rscript MVCell_V1_R_init.R path/to/input_file.h5ad 8
```

> The output `.rds` file will contain the processed SCE object for downstream annotation.

---

### **Step 2 — Run Cell Type Annotation**

Run the `MVCell_V1_Run_celltype.R` script to execute all three annotation methods:
- **ScType**
- **SingleR**
- **scATOMIC**

This script will:
- Read the processed `.rds` file from Step 1
- Run each annotation algorithm
- Save a combined annotation report compatible with MVCell

```bash
conda activate MVCell_base
Rscript MVCell_V1_Run_celltype.R path/to/input_file.h5ad 8
```

---

### **Output**

The final annotation report (`.csv`) can be fed directly into MVCell for further processing and analysis.

---

### **Summary Workflow**

1. **Prepare environments** — Follow MVCell installation instructions.
2. **Convert `.h5ad` → `.rds` (SCE)** — `MVCell_V1_R_init.R`
3. **Run cell type annotation** — `MVCell_V1_Run_celltype.R`
4. **Feed output to MVCell**

```text
.h5ad (counts) → [init script] → SCE (.rds) → [annotation script] → annotations.csv → MVCell
```


### 1. Data Preparation

Ensure your data files are organized as follows:
```
data/
├── train_set.tsv    # Training data
├── test_set.tsv     # Test data
└── HLCA_tree.csv    # Ontology file
```

**Expected data format:**
- Training/prediction files should contain:
  - Tool annotation columns: `{tool}_ann{level}` (e.g., `singleR_ann1`, `sctype_ann2`)
  - Confidence score columns: `{tool}_score` (e.g., `singleR_score`, `sctype_score`)
  - Ground truth columns: `ann_level_{level}` (e.g., `ann_level_1`, `ann_level_2`) - just for training /  Do not require when prediction

### 2. Prediction 

```bash

python predict.py \
    --model-dir ./pretrained_model \
    --data-dir ./data \
    --test-file ./data/test_set.tsv \ # you can replace your own data
    --output-dir ./results


## Command Line Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--model-dir` | str | Required | Directory containing trained model |
| `--data-dir` | str | `.` | Data directory path |
| `--batch-size` | int | `1024` | Prediction batch size |
| `--output-dir` | str | './results' | Output directory (defaults to model-dir) |
```

## Output Files

### Training Output
- `best_multilevel_seed{SEED}.pth`: Best model weights
- `encoders_in.pkl`: Input label encoders
- `encoders_out.pkl`: Output label encoders
- `bidirectional_transformer_curves.png`: Training curves visualization
- `loss_curve.png`: Loss curve plot

### Prediction Output
- `cell_type_predictions.csv`: cell type prediction
- `prediction_summary_report.txt`: summarise cell type prediction
- 
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
