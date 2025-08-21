# MVCell Pipeline Tutorial

This folder contains two concise tutorials:

- **A. Full Pipeline**: Run the entire MVCell pipeline from an `.h5ad` file.
- **B. Step 4 + 5 only**: Run HLCA transfer (Step 4) and Transformer prediction (Step 5) only.

> We use the demo file you created:
> - `balanced_noUnknown.h5ad` — remove unknown label

---

## Prerequisites

- Conda envs are created per the main README:
  - `MVCell` (main; R+Python)
  - `r-reticulate-mvcell` (reticulate helper for scATOMIC)
- Run these from the repo root: `MVCell/`

```bash
# From repo root
conda activate MVCell

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
cd ./  # repo root (MVCell)
chmod +x run_all.sh
./run_all.sh \
    --h5ad ./data/tutorial_data_set.h5ad \
    --cores 16
```

### 3. Outputs
- `celltype_annotate/_MVCell_annotate.csv` — R tools (SingleR, ScType, scATOMIC) per-cell outputs
- `celltype_annotate/merged_with_pipeline.tsv` — merged with CellTypist & SCimilarity
- `celltype_annotate/merged_annlevel_transfered.tsv` — HLCA ann_level_1..5 (Step 4 결과)
- `mvcell_results/` — Transformer outputs (Step 5 결과)

If you want a lighter demo, you can also start with `data/balanced_annL1_1500perlabel.h5ad`.

## B) Step 4 + 5 Only
Use this when you already have outputs from Step 3 or Step 4.

### Option 1 — Start from Step 3 output
**Input:** `celltype_annotate/merged_with_pipeline.tsv`

Run Step 4 (HLCA transfer) and Step 5 (Transformer) in one go:
```bash
# From repo root
chmod +x tutorial/run_step45.sh
tutorial/run_step45.sh \
    --merged-tsv ./celltype_annotate/merged_with_pipeline.tsv \
    --mapping-dir ./data/transform_db \
    --ontology ./data/HLCA_tree.csv \
    --model-dir ./pretrained_model \
    --outdir ./step45_results
```

This will produce:
- `./step45_results/merged_annlevel_transfered.tsv` — Step 4
- `./step45_results/transformer/` — Step 5 inference artifacts

### Option 2 — You already have the HLCA-transferred TSV (Step 4 done)
**Input:** `celltype_annotate/merged_annlevel_transfered.tsv`

Run only Step 5 (Transformer):
```bash
python predict.py \
    --model-dir ./pretrained_model \
    --data-dir ./data \
    --test-file ./celltype_annotate/merged_annlevel_transfered.tsv \
    --output-dir ./step5_results
```

## Notes & Tips

### Where does Step 4 live?
We call the same script that `run_all.sh` uses for Step 4 internally. If your repository uses a different path/name, open `run_all.sh` and copy that path into `tutorial/run_step45.sh` → `STEP4_SCRIPT` variable.

### Large data
Avoid committing `.h5ad` > 100 MiB to Git. Prefer Git LFS or release assets.

### Reproducibility
Fix `--seed` if your internal scripts support it.
