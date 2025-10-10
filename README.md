# MVCell

Automatic cell type annotation for single cell RNA data

## Bidirectional Hierarchical Transformer for Cell Type Annotation

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

### 2. Create the required conda environments

MVCell uses three conda environments:

- **MVCell** → Main environment containing both R and Python for running all celltype annotations (SingleR, sctype, scATOMIC, celltypist, scimilarity)
- **r-reticulate-mvcell** → Python environment for R's reticulate interface for scATOMIC

Create them from the provided YAML files:

```bash
# Create MVCell env
cd isntall
conda env create -f MVCell_env.yml   # or: mamba env create -f ...

# Create r-reticulate
conda env create -f r-reticulate-mvcell_env.yml    # or: mamba env create -f ...
```

**Why two environments?**

- Keeping r-reticulate-mvcell separate avoids Python version/package conflicts with the main MVCell environment.
- This separation also keeps r-reticulate alone for other R projects.

### 3. Install remaining R packages

```bash
#cd install
conda activate MVCell
Rscript MVCell_R_install_v2.R
```

## Quick Start

### 1) Full pipeline with run_all.sh (recommended)

`run_all.sh` runs the entire pipeline:

- **Step 1 (R)**: Convert .h5ad → *_sce.rds
- **Step 2 (R)**: Run SingleR / ScType / scATOMIC and export per-cell CSV
- **Step 3 (Python)**: Run CellTypist + SCimilarity and merge with Step2 CSV
- **Step 4 (Python)**: Transfer to HLCA ann_level_1..5
- **Step 5 (Python)**: MVCell prediction

```bash
# Move to MVCell main folder
cd ../

# Give execute permission to run_all.sh
chmod +x run_all.sh

# Run pipeline
./run_all.sh \
  --h5ad /path/to/your_dataset.h5ad \
  --cores 16 \
```

#### Options

- `--h5ad` (required): input .h5ad path
- `--cores`: CPU cores for R steps (default: available-2)
- `--outdir`: output directory (default: same as --h5ad dir)
- `--scimilarity-model`: SCimilarity model directory (default: ./scimilarity/model_v1.1)
  
  When the default path is used, the model tarball is auto-downloaded & extracted at first run.

- `--mapping-dir`: mapping DB for HLCA transfer (default: ./data/transform_db)
- `--ontology`: ontology CSV (default: ./data/HLCA_tree.csv)
- `--transformer-model`: transformer model dir (default: ./pretrained_model/)

#### Outputs

- `celltype_annotate/<BASENAME>_MVCell_annotate.csv` — R-based per-cell annotations
- `celltype_annotate/merged_with_pipeline.tsv` — merged with CellTypist & SCimilarity
- `celltype_annotate/merged_annlevel_transfered.tsv` — HLCA levels (ann_level_1..5)
- `mvcell_results/` — final transformer prediction artifacts



### 2) Standalone Transformer Prediction (predict.py only)

Run transformer inference directly on HLCA-transferred TSV (e.g., merged_annlevel_transfered.tsv):

```bash
python predict.py \
    --model-dir ./pretrained_model \
    --data-dir ./data \
    --test-file test_set.tsv \ # You can replace the data with your own files. Make sure to place them in the `data/` folder.
    --output-dir ./results
```

## Command Line Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--model-dir` | str | Required | Directory containing trained model |
| `--data-dir` | str | . | Data directory path |
| `--test-file` | str | test_data.tsv | Test data file (HLCA-transferred) |
| `--batch-size` | int | 1024 | Prediction batch size |
| `--seed` | int | 9401 | Random seed |
| `--output-dir` | str | ./results | Output directory |

## Data Layout (for training/inference)

```
data/
├── train_set.tsv            # (training)
├── test_set.tsv             # (inference)
└── HLCA_tree.csv            # ontology
```

**Expected columns (inference):**

- Tool scores: `{tool}_score` for tools in {singleR, sctype, scATOMIC, scimilarity, celltypist}
- (If present) Tool levels: `{tool}_ann{1..5}`
- (Optional) Ground truth levels for eval: `ann_level_{1..5}`

## GPU Requirements

- This project uses PyTorch built with CUDA 12.x.
- **Minimum NVIDIA GPU:** Compute Capability **≥ 7.0**.
- **Not supported:** Pascal (6.x), e.g. **GT 1030, GTX 1050/1060/1070/1080, Titan Xp**.
- **Supported families:** Volta/Turing/Ampere/Ada/Hopper and newer (7.x/8.x/9.x/…).
- If your GPU is older than CC 7.0, run on CPU (e.g., set `CUDA_VISIBLE_DEVICES=""`) or use a legacy PyTorch/CUDA stack.

## Citation If you use this code in your research, please cite:
```bibtex
@article{bidirectional_transformer_2025,
  title={Bidirectional Hierarchical Transformer for Cell Type Annotation},
  author={Your Name},
  journal={Journal Name},
  year={2025}
}
```

## License

MIT License — see LICENSE.

## Acknowledgments

- Thanks to the developers of the cell type annotation tools: singleR, sctype, scATOMIC, scimilarity, and celltypist 
- Human Lung Cell Atlas (HLCA) for providing the ontological framework
