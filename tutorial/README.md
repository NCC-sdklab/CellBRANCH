# MVCell Tutorials

This folder contains two concise tutorials:

- **A. Full Pipeline**: Run the entire MVCell pipeline from an `.h5ad` file.
- **B. Step 4 + 5 only**: Run HLCA transfer (Step 4) and Transformer prediction (Step 5) only.

> We use the demo files you created:
> - `balanced_noUnknown.h5ad` — Unknown 라벨 제외 버전
> - (Optional) `balanced_annL1_1500perlabel.h5ad` — ann_level_1 라벨별 1500개 추출 버전

---

## Prerequisites

- Conda envs are created per the main README:
  - `MVCell` (main; R+Python)
  - `r-reticulate-mvcell` (reticulate helper for scATOMIC)
- Run these from the repo root: `MVCell/`

```bash
# From repo root
conda activate MVCell
