# CellBRANCH

**Integrating heterogeneous cell-type annotators through masked bidirectional tree-structured hierarchical classification**

CellBRANCH integrates predictions from five independent cell-type annotation tools into a coherent path across five levels of the Human Lung Cell Atlas (HLCA) hierarchy. Tool-specific labels are first mapped to a shared naming system, after which a pretrained bidirectional hierarchical transformer combines their labels and confidence scores while enforcing valid parent-child transitions.

CellBRANCH is an integration model rather than a classifier that reads gene expression directly. The full pipeline in `run_all.sh` starts from an `.h5ad` count matrix, runs the five component annotators, harmonises their outputs, and then applies the CellBRANCH checkpoint.

## Installation

CellBRANCH is designed for Linux and uses two Conda environments. A CUDA-capable GPU is optional for transformer inference. The full component-annotation workflow also requires R and may take substantially longer on large datasets.

```bash
git clone https://github.com/NCC-sdklab/CellBRANCH.git
cd CellBRANCH

cd install
conda env create -f CellBRANCH_env.yml
conda env create -f r-reticulate-cellbranch_env.yml

conda activate CellBRANCH
Rscript CellBRANCH_R_install_v2.R
cd ..
```

Mamba can be used in place of Conda when creating the environments. The R installer should be run from the `install/` directory because it loads the accompanying dependency script by relative path.

The first full-pipeline run downloads the default SCimilarity model from Zenodo when `scimilarity/model_v1.1` is absent. Internet access is therefore required for initial environment setup and model download.

## Quick Start

### Steps 1-5: Full Pipeline from an h5ad File

The recommended entry point runs all five stages:

1. Convert the `.h5ad` object to a SingleCellExperiment object.
2. Run SingleR, ScType, and scATOMIC.
3. Run CellTypist and SCimilarity and merge all component outputs.
4. Map component labels to the five-level HLCA hierarchy.
5. Run CellBRANCH inference.

```bash
conda activate CellBRANCH

bash run_all.sh \
  --h5ad /path/to/dataset.h5ad \
  --cores 16 \
  --outdir results/dataset
```

The input must contain raw counts in the `counts` assay or in `X`. When only `X` is present, the pipeline treats it as raw counts. It creates `logcounts` when needed.

Useful options are:

| Option | Default | Purpose |
| --- | --- | --- |
| `--h5ad` | required | Input AnnData file |
| `--cores` | available cores minus two | CPU cores used by the R steps |
| `--outdir` | input file directory | Pipeline output directory |
| `--scimilarity-model` | `scimilarity/model_v1.1` | SCimilarity model directory |
| `--mapping-dir` | `data/transform_db` | Label-mapping tables |
| `--ontology` | `data/HLCA_tree.csv` | Five-level hierarchy |
| `--transformer-model` | `pretrained_model` | CellBRANCH checkpoint directory |

The full pipeline writes:

```text
<outdir>/
|-- celltype_annotate/
|   |-- <dataset>_CellBRANCH_annotate.csv
|   |-- merged_with_pipeline.tsv
|   `-- merged_annlevel_transfered.tsv
`-- cellbranch_results/
    |-- cell_type_predictions.csv
    `-- prediction_summary_report.txt
```

The spelling `merged_annlevel_transfered.tsv` is retained because it is the filename used by the current pipeline.

### Standalone CellBRANCH inference

Use `predict.py` when component predictions have already been mapped to the HLCA hierarchy. The bundled test set provides a reproducible example:

```bash
conda activate CellBRANCH

python predict.py \
  --model-dir pretrained_model \
  --data-dir data \
  --test-file test_set.tsv \
  --ont-file HLCA_tree.csv \
  --output-dir results/example
```

Key arguments are:

| Argument | Default | Purpose |
| --- | --- | --- |
| `--model-dir` | required | Directory containing the pretrained `.pt` checkpoint |
| `--data-dir` | `.` | Base directory for the input table and ontology |
| `--test-file` | `test_data.tsv` | Tab-separated mapped annotation table |
| `--ont-file` | `HLCA_tree.csv` | Hierarchy CSV relative to `--data-dir` |
| `--batch-size` | `1024` | Inference batch size |
| `--seed` | `9401` | Random seed |
| `--output-dir` | `results` | Prediction output directory |

### Steps 4-5: Harmonisation and CellBRANCH Inference

If raw outputs from the five component annotators are already available, Steps 4 and 5 can be run without repeating annotation:

```bash
conda activate CellBRANCH

bash tutorial/run_step45.sh \
  --merged-tsv data/step45_tutorial_data_set.tsv \
  --mapping-dir data/transform_db \
  --model-dir pretrained_model \
  --data-dir data \
  --outdir tutorial/step45_results
```

This writes the mapped table to `tutorial/step45_results/merged_annlevel_transfered.tsv` and final predictions to `tutorial/step45_results/step5_result/`.

## Workflow

```text
.h5ad count matrix
        |
        v
SingleR | ScType | scATOMIC | SCimilarity | CellTypist
        |
        v
Tool-specific labels and confidence scores
        |
        v
Mapping to the shared five-level HLCA hierarchy
        |
        v
CellBRANCH bidirectional hierarchical transformer
        |
        v
A valid cell-type path at levels 1-5
```

The model input contains five mapped labels from each annotator, giving 25 categorical features, plus one confidence score per annotator. During training, the forward and reverse decoders exchange information across hierarchy depth and are regularised for consistency. At inference, the masked forward decoder supplies the reported five-level path.

The bundled checkpoint uses a model dimension of 256, eight attention heads, three decoder layers, and dropout of 0.1007. The forward, reverse, and consistency loss weights are 2.0, 0.2, and 0.15, respectively.

## Input Contract

Standalone inference expects a tab-separated table with a `barcode` column. For each of the five annotators, it accepts the mapped columns `{tool}_ann1` through `{tool}_ann5` and a `{tool}_score` column, where `{tool}` is one of:

```text
singleR
sctype
scATOMIC
scimilarity
celltypist
```

Column names are case-sensitive. Missing annotation columns are represented internally as `Unknown`, and missing score columns are set to `0.0`. This permits inference with a missing component while retaining the same pretrained integration model. Optional `ann_level_1` through `ann_level_5` columns are used only when reference labels are available for evaluation.

For the Steps 4 and 5 workflow, `data/step45_tutorial_data_set.tsv` demonstrates the pre-mapping format: one `{tool}_annot` and one `{tool}_score` column per component, together with `barcode`.

## Repository Layout

```text
CellBRANCH/
|-- data/
|   |-- HLCA_tree.csv
|   |-- test_set.tsv
|   |-- step45_tutorial_data_set.tsv
|   `-- transform_db/
|-- install/
|   |-- CellBRANCH_env.yml
|   |-- r-reticulate-cellbranch_env.yml
|   `-- CellBRANCH_R_install_v2.R
|-- pretrained_model/
|   `-- model.pt
|-- tutorial/
|   `-- run_step45.sh
|-- run_all.sh
|-- transform_column.py
|-- predict.py
|-- model.py
`-- configs.py
```

## Reproducibility

The bundled checkpoint and hierarchy are the versions used for the current manuscript analyses:

```text
pretrained_model/model.pt
SHA-256: c102439710d11ab13efbe49163a7a30773dac928388ded7efe3aa3c6743d59fb

data/HLCA_tree.csv
SHA-256: a5e7fd0327a3ffa00af45e25d8974e117b082409a96d80bc49272c0f5cd4dbae
```

Verify them locally with:

```bash
sha256sum pretrained_model/model.pt data/HLCA_tree.csv
```

The included `data/test_set.tsv` contains 1,999 held-out cells. `data/step45_tutorial_data_set.tsv` contains 6,000 rows for demonstrating label harmonisation and inference.

## Citation

The manuscript citation will be added when it becomes available:

> CellBRANCH: Integrating heterogeneous cell-type annotators through masked bidirectional tree-structured hierarchical classification.

## Acknowledgements

We thank the developers of SingleR, ScType, scATOMIC, SCimilarity, and CellTypist, and the Human Lung Cell Atlas consortium for the reference atlas and cell-type hierarchy used in this work.
