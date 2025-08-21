# MVCell
Automatic Cell Type Annotation for Single-Cell RNA Data

## 1. Overview
MVCell provides a unified environment to run three R-based cell type annotation algorithms for single-cell RNA sequencing (scRNA-seq) data:

- **SingleR**
- **ScType**
- **scATOMIC**

These scripts are designed for reproducibility: package versions are fixed to ensure stability. 
The installation process has been tested on Linux-based systems. Other operating systems may require additional adjustments.


---
**🔹 Quick Facts**
- **Tested R version**: 4.4 (specified in `MVCell_R_env.yaml`).
- **Environment files**:  
  - `MVCell_R_env.yaml` – R + system dependencies for SingleR, ScType.  
  - `r-reticulate-mvcell_env.yaml` – Python dependencies for scATOMIC (via `reticulate`).
- **Input format**: `.h5ad` containing a **count matrix** (mandatory).  
  - Optionally: `logcounts` if you prefer your own log-normalization.  
  - Gene and cell metadata are **not used**.
- **References**: SingleR uses **Blueprint Encode Data**; ScType and scATOMIC use **default databases and parameters**.
- **Parallelization**:  
  - `[ncore]` argument controls cores for SingleR and ScType.  
  - scATOMIC is **always 1 core** for stability.  
  - Datasets >30,000 cells are automatically split for scATOMIC during imputation — no manual action required.
- **Reproducibility**: `sessionInfo()` is saved in the results directory after each run.

---
## 2. Requirements
- Linux-based system (tested on Ubuntu and CentOS)
- Conda or Mamba package manager
- R (version 4.4)
- Internet connection for package installation

## 3. Installation

1. **Clone the repository**  
```bash
git clone git@github.com:NCC-sdklab/MVCell.git
cd MVCell
```

2. **Create the conda environments**  
Two environments are prepared: `MVCell_R` for providing R environment, and `r-reticulate-mvcell` for imputation required for scATOMIC
```bash
cd r-install
conda env create -f MVCell_R_env.yaml 
conda env create -f r-reticulate-mvcell_env.yaml
conda activate MVCell_R
```

3. **Install R dependencies using `pak`**  
R-based script is prepared for consistent package management.
```r
Rscript MVCell_R_install_v2.R
```

## 4. Usage

### Step 1 — Preprocessing  
Run the first script to preprocess the `.h5ad` file into an R-compatible format:
```bash
Rscript MVCell_V2_R_init.R /path/to/input.h5ad [num_cores] [benchmark_id]
```

- `num_cores` *(optional)*: Number of CPU cores to use (default: all available minus 2)  
- `benchmark_id` *(optional)*: Identifier for benchmark tracking

### Step 2 — Annotation  
Run the second script to perform cell type annotation using all three algorithms:
```bash
Rscript MVCell_V2_Run_celltype.R /path/to/input.h5ad [num_cores] [benchmark_id]
```

### Step 3 — Output  
Results are saved as:
- Processed `.rds` file
- Per-method annotation results in `celltype_annotate/`
- Combined annotation CSV file

## 5. External intervention
When, whatever the reason, the celltype annotation does not work, or if you wish to perform it yourself, as long as you provide the correct name and the path, MVCell will integrate without issue.

### Step 1 - Preprocessing
`{sample_name}_sce.rds` is the file produced after `.h5ad` conversion. You may provide your version of `.rds` file. MVCell script checks for consistency; it needs to have `counts` and `logcounts` assays, and the dimensions must match the `.h5ad` file.

### Step 2 - Annotation
All celltype annotations are saved in the `celltype_annotate` folder located in the path where your `.h5ad` file is based.  

- For ScType, the outcome of `sctype_score()` is saved as `{sample_name}_sctype_result.rds`.  
- For SingleR, the outcome of `SingleR()` using `BlueprintEncodeData()` as a reference is saved as `{sample_name}_singleR_result.rds`.  
- For scATOMIC, the outcome of `create_summary_matrix()` is saved as `{sample_name}_scatomic_result.rds`.  


## 6. Notes
- Fixed package versions ensure stability but may not reflect upstream updates.
- `scATOMIC` requires both R and Python dependencies; the provided environment file includes these.
- The scripts automatically skip rerunning a method if previous results are found.

## 7. Troubleshooting
- If package installation fails, verify your network connection and that the CRAN and Bioconductor mirrors are accessible.
- For `scATOMIC` Python issues, ensure the correct conda environment is activated or set `RETICULATE_PYTHON`.

---
*Last updated: 2025-08-12*
