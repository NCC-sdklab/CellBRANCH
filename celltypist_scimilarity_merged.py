#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
from scimilarity.utils import lognorm_counts, align_dataset
from scimilarity import CellAnnotation

# ============== 공용 유틸 ==============

def _pick_first(cols: pd.Index, candidates: list[str]) -> str | None:
    return next((c for c in candidates if c in cols), None)

def _get_barcode_series(df: pd.DataFrame) -> pd.Series:
    for col in ["barcode", "cell_name", "CellID", "cell", "cell_id"]:
        if col in df.columns:
            return df[col].astype(str).rename("barcode")
    return pd.Series(df.index.astype(str), index=df.index, name="barcode")

def _resolve_h5ad(args) -> Path:
    if not args.h5ad:
        raise ValueError("Please provide --h5ad <file.h5ad>")
    p = Path(args.h5ad)
    if not p.is_file():
        raise FileNotFoundError(f"--h5ad not found: {p}")
    return p.resolve()

# ============== 병합 저장 ==============

def merge_and_save_results(
    tsv_path: str,
    celltypist_results: pd.DataFrame,
    scimilarity_results: pd.DataFrame,
    output_path: str
):
    db_df = pd.read_csv(tsv_path, sep=None, engine="python", dtype=str)
    if "barcode" not in db_df.columns:
        b = _get_barcode_series(db_df)
        db_df = db_df.copy()
        db_df.insert(0, "barcode", b.values)

    # CellTypist
    ct = celltypist_results.copy()
    if "barcode" not in ct.columns:
        raise KeyError("There are not 'barcode' in celltypist_results")
    ct_label_col = _pick_first(ct.columns, ["predicted_labels", "majority_voting", "cell_type"])
    ct_score_col = _pick_first(
        ct.columns,
        ["conf_score", "confidence", "predicted_labels_score", "majority_voting_score", "predicted_proba"]
    )
    if ct_label_col is None:
        raise KeyError(f"MVCell does not find the label column of CellTypist. Exsting column : {list(ct.columns)}")
    if ct_score_col is None:
        ct["celltypist_score"] = pd.NA
        ct_df = ct[["barcode", ct_label_col, "celltypist_score"]].rename(
            columns={ct_label_col: "celltypist_annot"}
        )
    else:
        ct_df = ct[["barcode", ct_label_col, ct_score_col]].rename(
            columns={ct_label_col: "celltypist_annot", ct_score_col: "celltypist_score"}
        )

    # SCimilarity
    scdf = scimilarity_results.copy()
    if "barcode" not in scdf.columns:
        raise KeyError("There are not 'barcode' in scimilarity_results")
    sc_label_col = _pick_first(scdf.columns, ["predictions_unconstrained", "prediction", "annot"])
    sc_score_col = _pick_first(scdf.columns, ["vsAll_weighted", "score", "conf"])
    if sc_label_col is None:
        raise KeyError(f"MVCell does not find the label column of SCimilarity. Exsting column : {list(scdf.columns)}")
    if sc_score_col is None:
        scdf["scimilarity_score"] = pd.NA
        sc_df = scdf[["barcode", sc_label_col, "scimilarity_score"]].rename(
            columns={sc_label_col: "scimilarity_annot"}
        )
    else:
        sc_df = scdf[["barcode", sc_label_col, sc_score_col]].rename(
            columns={sc_label_col: "scimilarity_annot", sc_score_col: "scimilarity_score"}
        )

    merged = db_df.merge(ct_df, on="barcode", how="left").merge(sc_df, on="barcode", how="left")
    merged = merged.rename(columns={
        "singleR_main_annot": "singleR_annot",
        "singleR_main_score": "singleR_score",
    })
    merged = merged[[c for c in merged.columns if c != "barcode"] + ["barcode"]]
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, sep="\t", index=False)
    print(f"Merged results saved to: {output_path}")

# ============== pipelines ==============

class CellTypistPipeline:
    def __init__(self, forced_model: str | None = None):
        self.model = None
        self.model_name = forced_model or "Cells_Lung_Airway.pkl"

    def load_model(self):
        forced = self.model_name
        try:
            print(f"Loading CellTypist model: {forced}")
            self.model = models.Model.load(model=forced)
        except Exception as e:
            print(f"[warn] Failed to load {forced} ({e}), falling back to Immune_All_Low.pkl")
            self.model = models.Model.load(model="Immune_All_Low.pkl")
            self.model_name = "Immune_All_Low.pkl"

    def _process(self, adata):
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)
        return adata

    def run(self, h5ad_file: Path) -> pd.DataFrame:
        if self.model is None:
            self.load_model()
        print(f"[CellTypist] {h5ad_file}")
        adata = sc.read(h5ad_file)
        adata = self._process(adata)
        pred = celltypist.annotate(adata, model=self.model, majority_voting=True).to_adata()
        df = pred.obs.copy()
        barcode = _get_barcode_series(df)
        label_col = _pick_first(df.columns, ["predicted_labels", "majority_voting", "cell_type"])
        score_col = _pick_first(
            df.columns,
            ["conf_score", "confidence", "predicted_labels_score", "majority_voting_score", "predicted_proba"]
        )
        if label_col is None:
            raise KeyError(f"CellTypist 결과 컬럼 누락: {list(df.columns)}")
        out = pd.DataFrame({
            "barcode": barcode.values,
            "predicted_labels": df[label_col].values,
            "celltypist_score": (df[score_col].values if score_col else pd.NA),
        })
        return out[["barcode", "predicted_labels", "celltypist_score"]]

class SCimilarityPipeline:
    def __init__(self, model_path: str):
        self.model_path = model_path
        self.ca = None

    def load_model(self):
        print(f"Loading SCimilarity model: {self.model_path}")
        self.ca = CellAnnotation(model_path=self.model_path)

    def _process(self, adata):
        adata = align_dataset(adata, self.ca.gene_order)
        if "counts" in adata.layers:
            adata = lognorm_counts(adata)
        return adata

    def _predict(self, adata):
        adata.obsm["X_scimilarity"] = self.ca.get_embeddings(adata.X)
        predictions, nn_idxs, nn_dists, nn_stats = self.ca.get_predictions_knn(adata.obsm["X_scimilarity"])
        adata.obs["predictions_unconstrained"] = predictions.values
        adata.obs["vsAll_weighted"] = nn_stats["vsAll_weighted"].values
        return adata

    def run(self, h5ad_file: Path) -> pd.DataFrame:
        if self.ca is None:
            self.load_model()
        print(f"[SCimilarity] {h5ad_file}")
        adata = sc.read(h5ad_file)
        adata = self._process(adata)
        adata = self._predict(adata)
        df = adata.obs.copy()
        barcode = _get_barcode_series(df)
        label_col = _pick_first(df.columns, ["predictions_unconstrained", "prediction", "annot"])
        score_col = _pick_first(df.columns, ["vsAll_weighted", "score", "conf"])
        if label_col is None:
            raise KeyError(f"SCimilarity 결과 컬럼 누락: {list(df.columns)}")
        out = pd.DataFrame({
            "barcode": barcode.values,
            "predictions_unconstrained": df[label_col].values,
            "vsAll_weighted": (df[score_col].values if score_col else pd.NA),
        })
        return out[["barcode", "predictions_unconstrained", "vsAll_weighted"]]

# ============== 메인 ==============

def parse_args():
    repo = Path(__file__).resolve().parent
    default_scim_path = repo / "scimilarity" / "model_v1.1"

    p = argparse.ArgumentParser(
        description="Run CellTypist + SCimilarity and merge with an existing per-cell TSV (no config & no auto-discovery)."
    )
    p.add_argument("input_tsv_path", help="Step2 output CSV/TSV (e.g., .../celltype_annotate/..._MVCell_annotate.csv)")
    p.add_argument("output_tsv_path", nargs="?", default=None,
                   help="Output TSV (default: merged_with_pipeline.tsv next to input)")
    p.add_argument("-m", "--scimilarity-model", default=default_scim_path,
                   help=f"SCimilarity 모델 디렉터리 (기본: {default_scim_path})")

    # ✅ 반드시 하나 지정
    p.add_argument("--h5ad", required=True, help="분석에 사용할 .h5ad 파일 경로(단일)")

    # ✅ CellTypist 모델 선택(기본은 Lung Airway)
    p.add_argument("--celltypist-model", default="Cells_Lung_Airway.pkl",
                   help="CellTypist .pkl 모델명 (기본: Cells_Lung_Airway.pkl)")

    return p.parse_args()

def main():
    args = parse_args()

    tsv_input = Path(args.input_tsv_path)
    if not tsv_input.is_file():
        raise FileNotFoundError(f"Input file not found: {tsv_input}")

    merged_output = str(Path(args.output_tsv_path) if args.output_tsv_path
                        else tsv_input.parent / "merged_with_pipeline.tsv")

    # h5ad 파일 (단일)
    h5ad_file = _resolve_h5ad(args)
    print(f"🔎 Using h5ad: {h5ad_file}")

    # 파이프라인 실행
    ct_pipe = CellTypistPipeline(forced_model=args.celltypist_model)
    sc_pipe = SCimilarityPipeline(model_path=str(args.scimilarity_model))

    print("\n" + "="*60 + "\n1) CellTypist running\n" + "="*60)
    ct_df = ct_pipe.run(h5ad_file)

    print("\n" + "="*60 + "\n2) SCimilarity running\n" + "="*60)
    sc_df = sc_pipe.run(h5ad_file)

    print("\n" + "="*60 + "\n3) Merging\n" + "="*60)
    merge_and_save_results(
        tsv_path=str(tsv_input),
        celltypist_results=ct_df,
        scimilarity_results=sc_df,
        output_path=merged_output
    )

    print("\n🎉 Done!")
    print(f"Input : {tsv_input}")
    print(f"Output: {merged_output}")

if __name__ == "__main__":
    main()
