#!/usr/bin/env bash
# Run Step 4 (HLCA transfer via transform_column.py) and Step 5 (Transformer) from any cwd.

set -euo pipefail

########## Resolve repo root & key paths ##########
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
REPO_ROOT="$( cd "${SCRIPT_DIR}/.." && pwd )"

STEP4_SCRIPT="${REPO_ROOT}/transform_column.py"   # <-- step4 파일명 반영
PREDICT_PY="${REPO_ROOT}/predict.py"

########## Defaults ##########
MERGED_TSV=""            # Step3 결과가 있으면 지정 (없으면 --config로 자동탐색)
CONFIG=""                # transform_column.py가 merged_with_pipeline.tsv를 자동 탐색하도록 할 때 사용
MAPPING_DIR=""           # 매핑 디렉토리(선택). 미지정 시 transform_column.py의 기본값 사용
TOOLS=""                 # 예: "singleR,sctype,scATOMIC,scimilarity,celltypist"
SEP="auto"               # 입력 구분자 (기본 auto)
OUT_SEP=""             # 출력 구분자 (기본 탭)
BARCODE_LAST=0           # 1이면 --barcode-last 전달

MODEL_DIR=""
OUTDIR="${SCRIPT_DIR}/step45_results"
DATA_DIR="${REPO_ROOT}/data"  # predict.py용; 필요 시 --data-dir로 교체

case "${OUT_SEP}" in
  "\\t") OUT_SEP=$'\t' ;;
  "\\n") OUT_SEP=$'\n' ;;
esac

########## Parse args ##########
while [[ $# -gt 0 ]]; do
  case "$1" in
    --merged-tsv) MERGED_TSV="$2"; shift 2;;
    --config) CONFIG="$2"; shift 2;;
    --mapping-dir) MAPPING_DIR="$2"; shift 2;;
    --tools) TOOLS="$2"; shift 2;;
    --sep) SEP="$2"; shift 2;;
    --out-sep) OUT_SEP="$2"; shift 2;;
    --barcode-last) BARCODE_LAST=1; shift 1;;
    --model-dir) MODEL_DIR="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --data-dir) DATA_DIR="$2"; shift 2;;
    *) echo "[Error] Unknown argument: $1"; exit 1;;
  esac
done

STEP4_OPTS=()
########## Validate minimal inputs ##########
[[ -z "${MODEL_DIR}" ]] && { echo "[Error] --model-dir is required"; exit 1; }
[[ -d "${MODEL_DIR}" ]] || { echo "[Error] model dir not found: ${MODEL_DIR}"; exit 1; }
[[ -d "${DATA_DIR}"  ]] || { echo "[Error] data dir not found: ${DATA_DIR}"; exit 1; }
[[ -f "${STEP4_SCRIPT}" ]] || { echo "[Error] Step4 script not found: ${STEP4_SCRIPT}"; exit 1; }
[[ -f "${PREDICT_PY}"  ]] || { echo "[Error] predict.py not found: ${PREDICT_PY}"; exit 1; }
[[ -n "${OUT_SEP}" ]] && STEP4_OPTS+=(--out-sep "${OUT_SEP}")

# merged-tsv 또는 config 둘 중 하나는 있어야 함
if [[ -z "${MERGED_TSV}" && -z "${CONFIG}" ]]; then
  echo "[Error] Provide either --merged-tsv or --config for Step 4 input resolution."
  exit 1
fi
# merged-tsv가 주어졌다면 실제 파일 확인
if [[ -n "${MERGED_TSV}" && ! -f "${MERGED_TSV}" ]]; then
  echo "[Error] merged TSV not found: ${MERGED_TSV}"
  exit 1
fi
# 매핑 디렉토리는 선택사항(주어졌다면만 확인)
if [[ -n "${MAPPING_DIR}" && ! -d "${MAPPING_DIR}" ]]; then
  echo "[Error] mapping dir not found: ${MAPPING_DIR}"
  exit 1
fi

########## Conda activate only if needed (avoid -u issues) ##########
if command -v conda >/dev/null 2>&1; then
  if [[ "${CONDA_DEFAULT_ENV:-}" != "CellBRANCH" && "${CONDA_PREFIX:-}" != *"/envs/CellBRANCH" ]]; then
    set +u
    eval "$(conda shell.bash hook)"
    conda activate CellBRANCH || true
    set -u
  fi
fi

mkdir -p "${OUTDIR}"

########## Step 4: HLCA transfer via transform_column.py ##########
HLCA_TSV="${OUTDIR}/merged_annlevel_transfered.tsv"
echo "[Step 4] HLCA transfer → ${HLCA_TSV}"

# 동적으로 옵션 구성
STEP4_OPTS=()
[[ -n "${MERGED_TSV}" ]] && STEP4_OPTS+=(--input "${MERGED_TSV}")
[[ -n "${CONFIG}"    ]] && STEP4_OPTS+=(--config "${CONFIG}")
[[ -n "${MAPPING_DIR}" ]] && STEP4_OPTS+=(--mapping-dir "${MAPPING_DIR}")
[[ -n "${TOOLS}"     ]] && STEP4_OPTS+=(--tools "${TOOLS}")
[[ "${SEP}" != "auto" ]] && STEP4_OPTS+=(--sep "${SEP}")
[[ -n "${OUT_SEP}"   ]] && STEP4_OPTS+=(--out-sep "${OUT_SEP}")
[[ "${BARCODE_LAST}" -eq 1 ]] && STEP4_OPTS+=(--barcode-last)

python "${STEP4_SCRIPT}" \
  "${STEP4_OPTS[@]}" \
  --output "${HLCA_TSV}"

########## Step 5: Transformer inference ##########

HLCA_TSV="${OUTDIR}/merged_annlevel_transfered.tsv"
HLCA_TSV_ABS="$(readlink -f "${HLCA_TSV}")"

MODEL_DIR_ABS="$(readlink -f "${MODEL_DIR}")"
DATA_DIR_ABS="$(readlink -f "${DATA_DIR}")"

# The checkpoint carries its own label encoders, so one file is the whole
# model. A directory prepared before that still works: predict.py falls back
# to best_multilevel_seed*.pth beside encoders_in.pkl and encoders_out.pkl.
if ! compgen -G "${MODEL_DIR_ABS}/*.pt" > /dev/null \
   && [[ ! -f "${MODEL_DIR_ABS}/encoders_in.pkl" ]]; then
  echo "[Error] no checkpoint found in: ${MODEL_DIR_ABS}"
  echo "        expected model.pt, or best_multilevel_seed*.pth with encoders_in.pkl"
  exit 1
fi

mkdir -p "${OUTDIR}/step5_result"

echo "[Step 5] Transformer prediction → ${OUTDIR}/step5_result"
python "${PREDICT_PY}" \
  --model-dir "${MODEL_DIR_ABS}" \
  --data-dir "${DATA_DIR_ABS}"  \
  --test-file "${HLCA_TSV_ABS}" \
  --output-dir "${OUTDIR}/step5_result"


echo "[Done] Step 4 + 5 finished. Outputs under: ${OUTDIR}"

