#!/usr/bin/env bash
set -euo pipefail

#bash batch_run.sh h5ad_list.txt

LIST_FILE=$1
CORES=8
SCIMIL_MODEL="/mnt/gluster_server/data/reference/SCimilarity/model_v1.1"
RUN_SCRIPT="./run_all.sh"   # adjust path if needed

if [[ ! -f "$LIST_FILE" ]]; then
  echo "List file not found: $LIST_FILE"
  exit 1
fi

while IFS= read -r H5AD; do
  [[ -z "$H5AD" ]] && continue   # skip empty lines

  echo ">>> Running $H5AD"
  bash "$RUN_SCRIPT" --h5ad "$H5AD" --cores "$CORES" --scimilarity-model "$SCIMIL_MODEL"

  echo "<<< Finished $H5AD"
  echo
done < "$LIST_FILE"