#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Transform tool-specific celltype annotations to HLCA ann_level_1..5 columns.

Usage examples:
  # 1) config 기반으로 merged_with_pipeline.tsv 자동 탐색 → 변환
  python transform_celltype_as_hlca.py --config ./cellbranch_parameter.csv

  # 2) 입력/출력 경로 명시
  python transform_celltype_as_hlca.py --input /path/to/merged_with_pipeline.tsv \
                                       --output /path/to/merged_annlevel_transfered.tsv

  # 3) 매핑 파일 디렉터리 통째로 지정(아래 파일명이 모두 있어야 함)
  python transform_celltype_as_hlca.py --mapping-dir /path/to/transform_each_sw_predict_to_hlca_celltype

Notes:
- 기본 입력 구분자는 자동 감지(, \t 등). --sep 로 강제 가능.
- 결과는 기본적으로 TSV로 저장(--out-sep 로 변경 가능).
"""

import os
import sys
import argparse
from pathlib import Path
import pandas as pd

# -------- 매핑 파일 기본(절대경로) --------
REPO_ROOT = Path(__file__).resolve().parent
# 환경변수로도 오버라이드 가능 (없으면 ./data)
DEFAULT_DATA_DIR = Path(os.environ.get("CELLBRANCH_DATA_DIR", REPO_ROOT / "data"))

# 매핑 파일 이름(디렉터리만 바꿔 받는 옵션에서 사용)
MAPPING_FILENAMES = {
    "singleR":    "SingleR_celltype_bt_hlca.csv",
    "sctype":     "sctype_celltype_bt_hlca.csv",
    "scATOMIC":   "scATOMIC_celltype_bt_hlca.csv",
    "scimilarity":"scimilarity_celltype_bt_hlca.csv",
    "celltypist": "celltypist_airway_bt_hlca.csv",
}

DEFAULT_MAPPING_DIR = DEFAULT_DATA_DIR / "transform_db"
DEFAULT_MAPPING_PATHS = {k: str(DEFAULT_MAPPING_DIR / v) for k, v in MAPPING_FILENAMES.items()}

# -------- Config 파서(경로 탐색에 사용) --------
class ConfigParser:
    @staticmethod
    def parse_config_file(config_path: str) -> dict:
        config_dict = {}
        df = pd.read_csv(config_path, sep=",", header=None, encoding="utf-8")
        for _, row in df.iterrows():
            if len(row) >= 2 and pd.notna(row.iloc[0]) and pd.notna(row.iloc[1]):
                key = str(row.iloc[0]).strip()
                value = str(row.iloc[1]).strip()
                if value.endswith(","):
                    value = value[:-1]
                config_dict[key] = value
        return config_dict

    @staticmethod
    def setup_directories(config: dict) -> dict:
        # 레포 상대 기본값으로 변경
        defaults = {
            "data_base": str(DEFAULT_DATA_DIR),
            "output_base": str(DEFAULT_DATA_DIR),
        }
        merged = {**defaults, **config}   # defaults로 채우고, config로 덮어씀
        return merged
# -------- 매핑 로더/도우미 --------
def build_mapping_paths(mapping_dir: str | None) -> dict[str, str]:
    if mapping_dir:
        base = Path(mapping_dir)
        return {k: str(base / fname) for k, fname in MAPPING_FILENAMES.items()}
    return DEFAULT_MAPPING_PATHS.copy()

def load_mapping_files(mapping_paths: dict[str, str]) -> dict[str, pd.DataFrame]:
    maps = {}
    for tool, path in mapping_paths.items():
        if not os.path.exists(path):
            raise FileNotFoundError(f"Mapping file not found: {path}")
        try:
            mp = pd.read_csv(path)
        except UnicodeDecodeError:
            mp = pd.read_csv(path, encoding="latin-1")
        if "original" not in mp.columns:
            raise KeyError(f"[{tool}] mapping file must have column 'original': {path}")
        for i in range(1, 6):
            col = f"ann_level_{i}"
            if col not in mp.columns:
                raise KeyError(f"[{tool}] mapping file must have column '{col}': {path}")
        maps[tool] = mp
    return maps

def map_to_levels(maps: dict[str, pd.DataFrame], tool: str, annot):
    if pd.isna(annot) or str(annot).strip() == "":
        return [None] * 5
    df_map = maps[tool]
    subset = df_map[df_map["original"] == annot]
    if subset.empty:
        return [None] * 5
    row = subset.iloc[0]
    return [row.get(f"ann_level_{i}") for i in range(1, 6)]

def resolve_col(df: pd.DataFrame, tool: str, kind: str = "annot") -> str | None:
    """
    kind: 'annot' or 'score'
    - singleR: singleR_main_annot/score 도 허용
    """
    if kind == "annot":
        candidates = [f"{tool}_annot"]
        if tool == "singleR":
            candidates = ["singleR_annot", "singleR_main_annot"]
    else:
        candidates = [f"{tool}_score"]
        if tool == "singleR":
            candidates = ["singleR_score", "singleR_main_score"]
    for c in candidates:
        if c in df.columns:
            return c
    return None

# -------- 입력 파일 자동 탐색 --------
def locate_input_file(config_path: str | None) -> Path:
    if not config_path:
        raise FileNotFoundError("No --input provided and no --config to auto-locate input.")
    cfg = ConfigParser.parse_config_file(config_path)
    cfg = ConfigParser.setup_directories(cfg)
    for key in ["output_base", "data_source", "tissue_source", "study_id"]:
        if key not in cfg:
            raise KeyError(f"config에 '{key}'가 없습니다. CONFIG={config_path}")

    base = Path(cfg["output_base"]) / cfg["data_source"] / cfg["tissue_source"] / cfg["study_id"]
    candidates = [
        base / "celltype_annotate" / "merged_with_pipeline.tsv",
        base / f"Processed/scRAP_{cfg.get('scrap_ver', '')}" / "merged_with_pipeline.tsv",
    ]
    for p in candidates:
        if p.is_file():
            return p
    hits = list(base.rglob("merged_with_pipeline*.tsv"))
    if hits:
        return hits[0]
    raise FileNotFoundError(f"merged_with_pipeline.tsv 를 찾을 수 없습니다. base={base}")

# -------- 본 변환 로직 --------
def apply_annlevel_transfer(
    in_df: pd.DataFrame,
    maps: dict[str, pd.DataFrame],
    tools: list[str],
    drop_src_annot: bool = True,
    move_barcode_last: bool = False
) -> pd.DataFrame:
    df = in_df.copy()

    # singleR_main_* → singleR_* 표준화
    df = df.rename(columns={
        "singleR_main_annot": "singleR_annot",
        "singleR_main_score": "singleR_score",
    })

    # 각 툴별 ann1~5 생성 및 원본 annot 삭제
    for tool in tools:
        src_col = resolve_col(df, tool, "annot")
        if not src_col:
            # 해당 툴 annot 없으면 스킵
            continue
        levels_df = df[src_col].map(lambda x: map_to_levels(maps, tool, x)).apply(pd.Series)
        levels_df.columns = [f"{tool}_ann{i}" for i in range(1, 6)]
        df = pd.concat([df, levels_df], axis=1)
        if drop_src_annot:
            df.drop(columns=[src_col], inplace=True)

    # 정렬: 툴별 ann1~5, score → 나머지 → (선택) barcode 맨 끝
    ordered = []
    for tool in tools:
        for i in range(1, 6):
            c = f"{tool}_ann{i}"
            if c in df.columns:
                ordered.append(c)
        sc = resolve_col(df, tool, "score")
        if sc and sc in df.columns:
            ordered.append(sc)

    remaining = [c for c in df.columns if c not in ordered]
    if move_barcode_last and "barcode" in remaining:
        remaining = [c for c in remaining if c != "barcode"] + ["barcode"]

    df = df[ordered + remaining]
    return df

# -------- 입출력(read/write) --------
def read_table(path: Path, sep: str) -> pd.DataFrame:
    if sep == "auto":
        return pd.read_csv(path, sep=None, engine="python", dtype=str)
    return pd.read_csv(path, sep=sep, dtype=str)

def write_table(df: pd.DataFrame, path: Path, sep: str):
    df.to_csv(path, sep=sep, index=False)

# -------- CLI --------
def parse_args():
    default_config = os.environ.get("PIPELINE_CONFIG", "/home/seokwon/nas1_deep/filled_celltype/cellbranch_parameter.csv")
    p = argparse.ArgumentParser(
        description="Transform tool-specific celltype annotations to HLCA ann levels (ann1..ann5)."
    )
    p.add_argument("--input", "-i", type=str, default=None,
                   help="입력 테이블 경로(미지정 시 --config로 merged_with_pipeline.tsv 자동 탐색).")
    p.add_argument("--config", "-c", type=str, default=default_config,
                   help=f"Config CSV 경로(기본: {default_config}). --input 미지정 시 자동 탐색에 사용.")
    p.add_argument("--output", "-o", type=str, default=None,
                   help="출력 경로(기본: 입력과 같은 폴더의 merged_annlevel_transfered.tsv).")
    p.add_argument("--sep", type=str, default="auto",
                   help="입력 구분자(기본: auto). 예: '\\t' 또는 ','.")
    p.add_argument("--out-sep", type=str, default="\t",
                   help="출력 구분자(기본: '\\t').")
    p.add_argument("--mapping-dir", type=str, default=None,
                   help="매핑 파일들이 들어있는 디렉터리(파일명은 고정). 지정 안 하면 내장 기본 경로 사용.")
    p.add_argument("--tools", type=str, default="singleR,sctype,scATOMIC,scimilarity,celltypist",
                   help="처리할 툴 목록(콤마 구분).")
    p.add_argument("--barcode-last", action="store_true",
                   help="출력에서 barcode 컬럼을 맨 끝으로 이동.")
    return p.parse_args()

def main():
    args = parse_args()

    # 입력 경로 결정
    if args.input:
        in_path = Path(args.input)
        if not in_path.is_file():
            raise FileNotFoundError(f"Input not found: {in_path}")
    else:
        in_path = locate_input_file(args.config)

    # 출력 경로 결정
    if args.output:
        out_path = Path(args.output)
    else:
        out_path = in_path.with_name("merged_annlevel_transfered.tsv")

    # 매핑 파일 로드
    mapping_paths = build_mapping_paths(args.mapping_dir)
    maps = load_mapping_files(mapping_paths)

    # 입력 읽기
    df_in = read_table(in_path, args.sep)

    # 변환 수행
    tools = [t.strip() for t in args.tools.split(",") if t.strip()]
    df_out = apply_annlevel_transfer(
        df_in,
        maps=maps,
        tools=tools,
        drop_src_annot=True,
        move_barcode_last=args.barcode_last
    )

    # 저장
    write_table(df_out, out_path, args.out_sep)  # <-- typo 방지 위해 아래로 수정

if __name__ == "__main__":
    try:
        # 작은 오타 방지(위에서 한 줄 교체)
        args = parse_args()
        if args.input:
            in_path = Path(args.input)
            if not in_path.is_file():
                raise FileNotFoundError(f"Input not found: {in_path}")
        else:
            in_path = locate_input_file(args.config)

        out_path = Path(args.output) if args.output else in_path.with_name("merged_annlevel_transfered.tsv")
        mapping_paths = build_mapping_paths(args.mapping_dir)
        maps = load_mapping_files(mapping_paths)

        df_in = read_table(in_path, args.sep)
        tools = [t.strip() for t in args.tools.split(",") if t.strip()]
        df_out = apply_annlevel_transfer(
            df_in,
            maps=maps,
            tools=tools,
            drop_src_annot=True,
            move_barcode_last=args.barcode_last
        )

        write_table(df_out, out_path, args.out_sep)
        print(f"[OK] Saved: {out_path}")
    except Exception as e:
        print(f"[ERROR] {e}")
        sys.exit(1)
