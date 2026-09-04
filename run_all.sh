#!/usr/bin/env bash
set -euo pipefail

#JHL: 제가 멍청해서.
# Ensure CellBRANCH conda environment is active
if [[ "$CONDA_DEFAULT_ENV" != "CellBRANCH" ]]; then
  echo "Error: Please activate the 'CellBRANCH' conda environment first."
  echo "Currently active: ${CONDA_DEFAULT_ENV:-none}"
  exit 1
fi


#############################################
# CellBRANCH full pipeline 
# - Step1: R preprocess (h5ad → *_sce.rds)
# - Step2: R annotate (SingleR, ScType, scATOMIC)  [single run]
# - Step3: python celltypist + scimilarity + merge
# - Step4: python transform → ann_level_1..5
# - Step5: python transformer prediction
#############################################

# -------- repo & script paths --------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# R scripts
R_STEP1="${REPO_ROOT}/CellBRANCH_V2_R_init.R"
R_STEP2="${REPO_ROOT}/CellBRANCH_V2_Run_celltype.R"

# Python scripts
PY_STEP3="${REPO_ROOT}/celltypist_scimilarity_merged.py"
PY_STEP4="${REPO_ROOT}/transform_column.py"
PY_STEP5="${REPO_ROOT}/predict.py"

# SCimilarity model bootstrap
SCIMIL_BASE="${REPO_ROOT}/scimilarity"
SCIMIL_TAR="${SCIMIL_BASE}/model_v1.1.tar.gz"
DEFAULT_SCIMIL_MODEL="${SCIMIL_BASE}/model_v1.1"
ZENO_URL="https://zenodo.org/records/10685499/files/model_v1.1.tar.gz?download=1"

# Other defaults
MAPPING_DIR="${REPO_ROOT}/data/transform_db"
ONTOLOGY_CSV="${REPO_ROOT}/data/HLCA_tree.csv"
TRANSF_MODEL_DIR="${REPO_ROOT}/pretrained_model/"
CELLTYP_MODEL="Cells_Lung_Airway.pkl"     # (Step3 내부에서도 강제 사용)
CORES="$(($(nproc) - 2))"
OUTDIR=""

usage() {
  cat <<EOF
Usage: $(basename "$0") --h5ad /path/to/file.h5ad [options]

Required:
  --h5ad PATH                 Input .h5ad file

Options:
  --cores N                   CPU cores for R steps (default: $(($(nproc)-2)))
  --outdir DIR                Output directory (default: same dir as --h5ad)
  --scimilarity-model DIR     SCimilarity model dir (default: ${DEFAULT_SCIMIL_MODEL})
  --mapping-dir DIR           Mapping DB dir for transform (default: ${MAPPING_DIR})
  --ontology CSV              Ontology CSV for prediction (default: ${ONTOLOGY_CSV})
  --transformer-model DIR     Transformer model dir (default: ${TRANSF_MODEL_DIR})
  -h, --help                  Show this help

Notes:
- 이 스크립트는 **청크 분할 없이** Step2를 그대로 실행합니다.
- SCimilarity 모델은 기본 경로 사용 시 자동 다운로드/해제합니다.
EOF
  exit 1
}

# -------- parse args --------
H5AD=""
SCIMIL_MODEL="${DEFAULT_SCIMIL_MODEL}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --h5ad) H5AD="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --scimilarity-model) SCIMIL_MODEL="$2"; shift 2 ;;
    --celltypist-model) CELLTYP_MODEL="$2"; shift 2 ;;
    --mapping-dir) MAPPING_DIR="$2"; shift 2 ;;
    --ontology) ONTOLOGY_CSV="$2"; shift 2 ;;
    --transformer-model) TRANSF_MODEL_DIR="$2"; shift 2 ;;
    --benchmark-num) BENCH_NUM="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

#adding benchmark sequences for testing. May be removed later
BENCH_NUM="${BENCH_NUM:-NA}"

[[ -z "${H5AD}" ]] && { echo "ERROR: --h5ad is required"; usage; }
[[ -f "${H5AD}" ]] || { echo "ERROR: h5ad not found: ${H5AD}"; exit 1; }
[[ -z "${OUTDIR}" ]] && OUTDIR="$(dirname "${H5AD}")"
BASENAME="$(basename "${H5AD}")"
BASENAME="${BASENAME%.h5ad}"

# -------- SCimilarity model bootstrap (only if default location is used) --------
if [[ "${SCIMIL_MODEL}" == "${DEFAULT_SCIMIL_MODEL}" ]]; then
  mkdir -p "${SCIMIL_BASE}"
  if [[ ! -d "${SCIMIL_MODEL}" || -z "$(ls -A "${SCIMIL_MODEL}" 2>/dev/null || true)" ]]; then
    echo "[Bootstrap] SCimilarity model not found at ${SCIMIL_MODEL}"
    if [[ ! -f "${SCIMIL_TAR}" ]]; then
      echo "[Bootstrap] Downloading model_v1.1 → ${SCIMIL_TAR}"
      command -v curl >/dev/null 2>&1 || { echo "ERROR: curl not found"; exit 1; }
      curl -fL -o "${SCIMIL_TAR}" "${ZENO_URL}"
    else
      echo "[Bootstrap] Found tarball: ${SCIMIL_TAR}"
    fi
    echo "[Bootstrap] Extracting tarball to ${SCIMIL_BASE}"
    tar -xzvf "${SCIMIL_TAR}" -C "${SCIMIL_BASE}"
    if [[ ! -d "${SCIMIL_MODEL}" ]]; then
      CAND="$(find "${SCIMIL_BASE}" -maxdepth 1 -type d -name 'model_*' -printf '%T@ %p\n' 2>/dev/null | sort -nr | head -n1 | awk '{print $2}')"
      if [[ -n "${CAND}" ]]; then
        SCIMIL_MODEL="${CAND}"
        echo "[Bootstrap] Using detected model dir: ${SCIMIL_MODEL}"
      fi
    fi
  fi
fi

echo "========== CellBRANCH Full Pipeline =========="
echo "Repo root          : ${REPO_ROOT}"
echo "H5AD               : ${H5AD}"
echo "Cores              : ${CORES}"
echo "Output dir         : ${OUTDIR}"
echo "SCimilarity model  : ${SCIMIL_MODEL}"
echo "CellTypist model   : ${CELLTYP_MODEL}"
echo "Mapping dir        : ${MAPPING_DIR}"
echo "Ontology CSV       : ${ONTOLOGY_CSV}"
echo "Transformer model  : ${TRANSF_MODEL_DIR}"
echo "====================================================="

SCE_RDS="${OUTDIR}/${BASENAME}_sce.rds"
CELLTYPE_DIR="${OUTDIR}/celltype_annotate"
mkdir -p "${CELLTYPE_DIR}"

# ------------------------------------------
# Step 1: R preprocess  (.h5ad -> *_sce.rds)
# ------------------------------------------
echo "[Step1] R preprocess → SCE RDS"
Rscript "${R_STEP1}" "${H5AD}" "${CORES}" "${BENCH_NUM}"
[[ -f "${SCE_RDS}" ]] || { echo "ERROR: Step1 output not found: ${SCE_RDS}"; exit 1; }

# ------------------------------------------
# Step 2: R annotate (SingleR/ScType/scATOMIC) — single run
# ------------------------------------------
echo "[Step2] R annotate"
Rscript "${R_STEP2}" "${H5AD}" "${CORES}" "${BENCH_NUM}"

MERGED_CSV="${CELLTYPE_DIR}/${BASENAME}_CellBRANCH_annotate.csv"
[[ -f "${MERGED_CSV}" ]] || {
  # 혹시 사용자가 full 같은 이름으로 바꿨을 때 대비(드물지만)
  ALT="${CELLTYPE_DIR}/${BASENAME}_full_CellBRANCH_annotate.csv"
  if [[ -f "${ALT}" ]]; then
    ln -sf "$(basename "${ALT}")" "${MERGED_CSV}"
  else
    echo "ERROR: Step2 output not found: ${MERGED_CSV}"
    exit 1
  fi
}

# ------------------------------------------
# Step 3: python celltypist + SCimilarity + merge
# ------------------------------------------
STEP3_OUT="${CELLTYPE_DIR}/merged_with_pipeline.tsv"
echo "[Step3] CellTypist + SCimilarity + merge"
python "${PY_STEP3}" "${MERGED_CSV}" "${STEP3_OUT}" \
  --scimilarity-model "${SCIMIL_MODEL}" \
  --h5ad "${H5AD}" \
  --celltypist-model "${CELLTYP_MODEL}"
[[ -f "${STEP3_OUT}" ]] || { echo "ERROR: Step3 output not found: ${STEP3_OUT}"; exit 1; }

# ------------------------------------------
# Step 4: transform → HLCA ann_level_1..5
# ------------------------------------------
STEP4_OUT="${CELLTYPE_DIR}/merged_annlevel_transfered.tsv"
echo "[Step4] Transform annotations → ann_level_1..5"
python "${PY_STEP4}" --input "${STEP3_OUT}" \
  --output "${STEP4_OUT}" \
  --mapping-dir "${MAPPING_DIR}"
[[ -f "${STEP4_OUT}" ]] || { echo "ERROR: Step4 output not found: ${STEP4_OUT}"; exit 1; }

# ------------------------------------------
# Step 5: transformer prediction
# ------------------------------------------
PRED_OUTDIR="${OUTDIR}/cellbranch_results"
mkdir -p "${PRED_OUTDIR}"

echo "[Step5] Transformer prediction"
python "${PY_STEP5}" \
  --model-dir "${TRANSF_MODEL_DIR}" \
  --data-dir "$(dirname "${STEP4_OUT}")" \
  --test-file "$(basename "${STEP4_OUT}")" \
  --ont-file "${ONTOLOGY_CSV}" \
  --output-dir "${PRED_OUTDIR}"

echo ""
echo "🎉 Done !"
echo " - Step2 CSV            : ${MERGED_CSV}"
echo " - Step3 merged TSV     : ${STEP3_OUT}"
echo " - Step4 annlevel TSV   : ${STEP4_OUT}"
echo " - Step5 results dir    : ${PRED_OUTDIR}"
