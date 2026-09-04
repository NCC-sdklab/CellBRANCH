# CellBRANCH Pipeline Tutorial

This folder contains two concise tutorials:

- **A. Full Pipeline**: Run the entire CellBRANCH pipeline from an `.h5ad` file.
- **B. Step 4 + 5 only**: Run HLCA transfer (Step 4) and CellBRANCH prediction (Step 5) only.

> We use the demo file you created:
> - `tutorial_data_set.h5ad` — remove unknown label



## Prerequisites

- Conda envs are created per the main README:
  - `CellBRANCH` (main; R+Python)
  - `r-reticulate-cellbranch` (reticulate helper for scATOMIC)
- Run these from the repo root: `CellBRANCH/`

```bash
# From repo root
conda activate CellBRANCH

```

## A) Full Pipeline Tutorial (Steps 1 → 5)
This will run the entire pipeline using `run_all.sh`.

### 1. Prepare input
Place your `.h5ad` under a working directory, e.g.:
```
data/
└── tutorial_data_set.h5ad
```

### 2. Run
```bash
cd ./  # repo root (CellBRANCH)
chmod +x run_all.sh
./run_all.sh \
    --h5ad ./data/tutorial_data_set.h5ad \
    --cores 16
```

### 3. Outputs
- `celltype_annotate/_CellBRANCH_annotate.csv` — R tools (SingleR, ScType, scATOMIC) per-cell outputs
- `celltype_annotate/merged_with_pipeline.tsv` — merged with CellTypist & SCimilarity
- `celltype_annotate/merged_annlevel_transfered.tsv` — HLCA ann_level_1..5 (The result of Step 4)
- `cellbranch_results/` — CellBRANCH outputs (The result of Step 5)

## B) Step 4 + 5 Only
Use this when you already have outputs from Step 3 or Step 4.

**Input:** `celltype_annotate/merged_with_pipeline.tsv`

Run Step 4 (HLCA transfer) and Step 5 (Transformer) in one go:
```bash
# From repo root
chmod +x tutorial/run_step45.sh
tutorial/run_step45.sh \
    --merged-tsv ../data/step45_tutorial_data_set.tsv \
    --mapping-dir ../data/transform_db \
    --model-dir ../pretrained_model \
    --outdir ./step45_results
```

This will produce:
- `./merged_annlevel_transfered.tsv` — Step 4 result
- `./step45_results/prediction_summary_report.txt` — The summary of step 5 result
- `./step45_results/cell_type_predictions.csv` — Step 5 inference artifacts

## Notes & Tips


### Large data
Avoid committing `.h5ad` > 100 MiB to Git. Prefer Git LFS or release assets.

### Reproducibility
Fix `--seed` if your internal scripts support it.
