
---

# 2) `tutorial/run_step45.sh` (새 파일)

아래 스크립트는 **Step 4(Transfer to HLCA) → Step 5(Transformer)**를 한 번에 실행합니다.  
**중요:** `STEP4_SCRIPT` 위치는 실제 repo에서 `run_all.sh`가 호출하는 Step4 스크립트 경로로 수정하세요.

```bash
#!/usr/bin/env bash
# Run Step 4 (HLCA transfer) and Step 5 (Transformer) in one go.
# Usage:
#   tutorial/run_step45.sh \
#     --merged-tsv ./celltype_annotate/merged_with_pipeline.tsv \
#     --mapping-dir ./data/transform_db \
#     --ontology ./data/HLCA_tree.csv \
#     --model-dir ./pretrained_model \
#     --outdir ./step45_results
#
# Notes:
# - If you already have the HLCA-transferred TSV, skip this script and run predict.py directly.

set -euo pipefail

# ---------- EDIT HERE IF NEEDED ----------
# Path to the Step 4 Python script used by run_all.sh.
# Open run_all.sh in the repo root and copy the exact script path here.
STEP4_SCRIPT="./scripts/transfer_to_hlca.py"
# -----------------------------------------

MERGED_TSV=""
MAPPING_DIR=""
ONTOLOGY=""
MODEL_DIR=""
OUTDIR="./step45_results"

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --merged-tsv) MERGED_TSV="$2"; shift 2;;
    --mapping-dir) MAPPING_DIR="$2"; shift 2;;
    --ontology) ONTOLOGY="$2"; shift 2;;
    --model-dir) MODEL_DIR="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    *) echo "[Error] Unknown argument: $1"; exit 1;;
  esac
done

# Validate
[[ -z "$MERGED_TSV" ]] && { echo "[Error] --merged-tsv is required"; exit 1; }
[[ -z "$MAPPING_DIR" ]] && { echo "[Error] --mapping-dir is required"; exit 1; }
[[ -z "$ONTOLOGY" ]] && { echo "[Error] --ontology is required"; exit 1; }
[[ -z "$MODEL_DIR" ]] && { echo "[Error] --model-dir is required"; exit 1; }
[[ -f "$MERGED_TSV" ]] || { echo "[Error] merged TSV not found: $MERGED_TSV"; exit 1; }
[[ -d "$MAPPING_DIR" ]] || { echo "[Error] mapping dir not found: $MAPPING_DIR"; exit 1; }
[[ -f "$ONTOLOGY" ]] || { echo "[Error] ontology CSV not found: $ONTOLOGY"; exit 1; }
[[ -d "$MODEL_DIR" ]] || { echo "[Error] model dir not found: $MODEL_DIR"; exit 1; }

# Activate conda env (portable way)
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate MVCell || true
fi

mkdir -p "$OUTDIR"

# -------- Step 4: HLCA transfer --------
HLCA_TSV="${OUTDIR}/merged_annlevel_transfered.tsv"
echo "[Step 4] HLCA transfer → ${HLCA_TSV}"
python "${STEP4_SCRIPT}" \
  --input "${MERGED_TSV}" \
  --mapping-dir "${MAPPING_DIR}" \
  --ontology "${ONTOLOGY}" \
  --output "${HLCA_TSV}"

# -------- Step 5: Transformer inference --------
echo "[Step 5] Transformer prediction → ${OUTDIR}/transformer"
python predict.py \
  --model-dir "${MODEL_DIR}" \
  --data-dir ./data \
  --test-file "${HLCA_TSV}" \
  --output-dir "${OUTDIR}/transformer"

echo "[Done] Step 4 + 5 finished. Outputs under: ${OUTDIR}"
