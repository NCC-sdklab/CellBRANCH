#!/usr/bin/env python3
# %%
import os
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
import matplotlib.pyplot as plt
from pathlib import Path
from scimilarity.utils import lognorm_counts, align_dataset
from scimilarity import CellAnnotation
import openpyxl
import sys
import argparse

# =========================
# 공용 유틸
# =========================
def _resolve_data_dir(config: dict) -> Path:
    """
    data_dir 후보 둘 중 실제로 존재하는 경로를 반환.
    우선순위: 1) Processed/scRAP_{scrap_ver} 2) celltype_annotate
    둘 다 없으면 FileNotFoundError 발생.
    """
    base = Path(config["output_base"]) / config["data_source"] / config["tissue_source"] / config["study_id"]

    scrap_ver = config.get("scrap_ver")
    primary = base / f"Processed/scRAP_{scrap_ver}" if scrap_ver else base / "Processed"
    secondary = base / "celltype_annotate"

    if primary.is_dir():
        return primary
    if secondary.is_dir():
        return secondary

    raise FileNotFoundError(
        "data_dir 후보 경로가 존재하지 않습니다.\n"
        f"- 1차 후보: {primary}\n"
        f"- 2차 후보: {secondary}\n"
        "config 값(data_source, tissue_source, study_id, scrap_ver, output_base)을 확인하세요."
    )


def _find_h5ad_files(root: Path, prefer_name: str = "sce_qc.h5ad") -> list[Path]:
    """
    root 하위에서 재귀적으로 *.h5ad 탐색.
    없으면 바로 상위 폴더에서 한 번 더 재귀 탐색.
    prefer_name이 있으면 그 파일들만 우선 사용, 없으면 전체 *.h5ad 사용.
    """
    def _search(dirpath: Path) -> list[Path]:
        if not dirpath.exists():
            return []
        all_h5ad = list(dirpath.rglob("*.h5ad"))
        if not all_h5ad:
            return []
        preferred = [p for p in all_h5ad if p.name == prefer_name]
        return preferred if preferred else all_h5ad

    # 1) root에서 검색
    results = _search(root)
    if results:
        return results

    # 2) 상위 폴더에서 한 번 더 검색
    parent = root.parent
    if parent and parent != root:
        return _search(parent)

    return []

def _get_barcode_series(df: pd.DataFrame) -> pd.Series:
    """
    obs DataFrame에서 바코드 시리즈를 안전하게 추출.
    우선순위: 'barcode' -> 'cell_name' -> 'CellID' -> 'cell' -> 'cell_id' -> index
    """
    for col in ["barcode", "cell_name", "CellID", "cell", "cell_id"]:
        if col in df.columns:
            return df[col].astype(str).rename("barcode")
    # 컬럼에 없으면 index에서 가져오기
    return pd.Series(df.index.astype(str), index=df.index, name="barcode")

def _pick_first(cols: pd.Index, candidates: list[str]) -> str | None:
    """cols 안에서 candidates 리스트의 첫 번째로 존재하는 컬럼명을 반환"""
    return next((c for c in candidates if c in cols), None)

def merge_and_save_results(
    tsv_path: str,
    celltypist_results: pd.DataFrame,
    scimilarity_results: pd.DataFrame,
    output_path: str
):
    """
    barcode 컬럼을 기준으로 메모리 내 DataFrame과 기존 TSV 데이터를 병합해 저장합니다.
    """
    # db_df = pd.read_csv(tsv_path, sep='\t', dtype={'barcode': str})
    db_df = pd.read_csv(tsv_path, sep=None, engine="python", dtype=str)
    if "barcode" not in db_df.columns:
        b = _get_barcode_series(db_df)  # 'barcode' / 'cell_name' / 'CellID' / 'cell' / 'cell_id' / index 순으로 추출
        db_df = db_df.copy()
        db_df.insert(0, "barcode", b.values)

    # ---- CellTypist 결과 컬럼 안전 매핑 ----
    ct = celltypist_results.copy()
    if "barcode" not in ct.columns:
        raise KeyError("celltypist_results에 'barcode' 컬럼이 없습니다.")
    ct_label_col = _pick_first(ct.columns, ["predicted_labels", "majority_voting", "cell_type"])
    ct_score_col = _pick_first(ct.columns, ["confidence", "predicted_labels_score", "majority_voting_score", "predicted_proba"])
    if ct_label_col is None:
        raise KeyError(f"CellTypist 라벨 컬럼을 찾지 못했습니다. 존재 컬럼: {list(ct.columns)}")
    if ct_score_col is None:
        # 점수 없으면 NaN으로 채움
        ct["celltypist_score"] = pd.NA
        ct_df = ct[["barcode", ct_label_col, "celltypist_score"]].rename(
            columns={ct_label_col: "celltypist_annot"}
        )
    else:
        ct_df = ct[["barcode", ct_label_col, ct_score_col]].rename(
            columns={ct_label_col: "celltypist_annot", ct_score_col: "celltypist_score"}
        )

    # ---- SCimilarity 결과 컬럼 안전 매핑 ----
    scdf = scimilarity_results.copy()
    if "barcode" not in scdf.columns:
        raise KeyError("scimilarity_results에 'barcode' 컬럼이 없습니다.")
    sc_label_col = _pick_first(scdf.columns, ["predictions_unconstrained", "prediction", "annot"])
    sc_score_col = _pick_first(scdf.columns, ["vsAll_weighted", "score", "conf"])
    if sc_label_col is None:
        raise KeyError(f"SCimilarity 라벨 컬럼을 찾지 못했습니다. 존재 컬럼: {list(scdf.columns)}")
    if sc_score_col is None:
        scdf["scimilarity_score"] = pd.NA
        sc_df = scdf[["barcode", sc_label_col, "scimilarity_score"]].rename(
            columns={sc_label_col: "scimilarity_annot"}
        )
    else:
        sc_df = scdf[["barcode", sc_label_col, sc_score_col]].rename(
            columns={sc_label_col: "scimilarity_annot", sc_score_col: "scimilarity_score"}
        )

    merged = db_df.merge(ct_df, on='barcode', how='left').merge(sc_df, on='barcode', how='left')
    merged = merged.rename(columns={
        "singleR_main_annot": "singleR_annot",
        "singleR_main_score": "singleR_score",
        })

    merged = merged[[c for c in merged.columns if c != 'barcode'] + ['barcode']]
    
    merged.to_csv(output_path, sep='\t', index=False)
    print(f"Merged results saved to: {output_path}")

# =========================
# Config
# =========================
class ConfigParser:
    """Config 파일 파싱을 위한 클래스"""
    @staticmethod
    def parse_config_file(config_path):
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
    def get_tissue_and_organism(config_path):
        config = ConfigParser.parse_config_file(config_path)
        tissue_source = config.get("tissue_source", None)
        org_type = config.get("org_type", None)
        return tissue_source, org_type

# =========================
# CellTypist 모델 선택기
# =========================
class CellTypistModelSelector:
    def __init__(self):
        self.organism_mapping = {
            "hs": ["Human", "Adult_Human", "Fetal_Human", "Healthy_Human", "Cells_Human"],
            "ms": ["Mouse", "Adult_Mouse", "Healthy_Mouse", "Developing_Mouse"],
        }
        self.tissue_mapping = {
            "colon": ["Intestinal", "Gut", "Colorectal", "Colon"],
            "lung": ["Lung", "Airway", "Pulmonary"],
            "breast": ["Breast", "Mammary"],
            "brain": ["Brain", "Neural", "Hippocampus", "Cortex"],
            "liver": ["Liver", "Hepatic"],
            "kidney": ["Kidney", "Renal"],
            "heart": ["Heart", "Cardiac"],
            "skin": ["Skin", "Dermal"],
            "pancreas": ["Pancreas", "Pancreatic"],
            "blood": ["Blood", "PBMC", "Immune"],
        }

    def find_best_model(self, tissue_name, organism_type, log_file_path="celltypist_model_selection_log.csv"):
        available_models = pd.DataFrame(models.models_description())
        tissue_lower = tissue_name.lower()
        organism_keywords = self.organism_mapping.get(organism_type, [])
        tissue_keywords = self.tissue_mapping.get(tissue_lower, [tissue_name])

        log_data = []
        model_scores = {}
        for model in available_models["model"]:
            score = 0
            model_lower = model.lower()
            log_entry = self._initialize_log_entry(model)
            score += self._calculate_organism_score(model_lower, organism_keywords, log_entry)
            score += self._calculate_tissue_score(model_lower, tissue_keywords, log_entry)
            score += self._calculate_special_score(model_lower, tissue_lower, log_entry)
            score += self._calculate_age_score(model_lower, log_entry)
            score += self._calculate_penalty_score(model_lower, log_entry)
            model_scores[model] = score
            log_entry["total_score"] = score
            log_data.append(log_entry)

        self._save_log(log_data, log_file_path)
        return self._select_best_model(model_scores)

    def _initialize_log_entry(self, model):
        return {
            "model_name": model,
            "model_parts": model.replace(".pkl", ""),
            "organism_match": False,
            "tissue_match": False,
            "special_match": False,
            "age_type": "none",
            "organism_score": 0,
            "tissue_score": 0,
            "special_score": 0,
            "age_score": 0,
            "penalty_score": 0,
            "total_score": 0,
            "matched_organism_keyword": "",
            "matched_tissue_keyword": "",
        }

    def _calculate_organism_score(self, model_lower, organism_keywords, log_entry):
        for org_keyword in organism_keywords:
            if org_keyword.lower() in model_lower:
                log_entry["organism_match"] = True
                log_entry["organism_score"] = 10
                log_entry["matched_organism_keyword"] = org_keyword
                return 10
        return 0

    def _calculate_tissue_score(self, model_lower, tissue_keywords, log_entry):
        for tissue_keyword in tissue_keywords:
            if tissue_keyword.lower() in model_lower:
                log_entry["tissue_match"] = True
                log_entry["tissue_score"] = 8
                log_entry["matched_tissue_keyword"] = tissue_keyword
                return 8
        return 0

    def _calculate_special_score(self, model_lower, tissue_lower, log_entry):
        if tissue_lower == "colon":
            if any(keyword in model_lower for keyword in ["intestinal", "gut", "colorectal"]):
                log_entry["special_match"] = True
                log_entry["special_score"] = 5
                return 5
        return 0

    def _calculate_age_score(self, model_lower, log_entry):
        if "adult" in model_lower:
            log_entry["age_type"] = "adult"
            log_entry["age_score"] = 2
            return 2
        elif any(keyword in model_lower for keyword in ["developing", "fetal"]):
            log_entry["age_type"] = "developing/fetal"
            log_entry["age_score"] = 1
            return 1
        return 0

    def _calculate_penalty_score(self, model_lower, log_entry):
        if "immune_all" in model_lower:
            log_entry["penalty_score"] = -5
            return -5
        return 0

    def _save_log(self, log_data, log_file_path):
        log_df = pd.DataFrame(log_data).sort_values("total_score", ascending=False)
        log_df.to_csv(log_file_path, index=False)

    def _select_best_model(self, model_scores):
        if model_scores:
            best_model = max(model_scores, key=model_scores.get)
            if model_scores[best_model] > 0:
                return best_model
        return "Immune_All_Low.pkl"

# =========================
# CellTypist 파이프라인
# =========================
class CellTypistPipeline:
    """CellTypist 파이프라인 메인 클래스"""

    def __init__(self, config_path):
        self.config_path = config_path
        self.model_selector = CellTypistModelSelector()
        self.model = None
        self.model_name = None
        self.all_results = None

    def load_model_from_config(self, log_dir=None):
        tissue_source, org_type = ConfigParser.get_tissue_and_organism(self.config_path)
        if not tissue_source:
            print("Warning: tissue_source not found, assuming immune tissue")
            tissue_source = "immune"
        if not org_type:
            print("Warning: org_type not found, assuming human (hs)")
            org_type = "hs"

        log_file_path = Path(log_dir)/"celltypist_model_selection_log.csv" if log_dir else Path("celltypist_model_selection_log.csv")
        _ = self.model_selector.find_best_model(tissue_source, org_type, str(log_file_path))

        print("⚠️" * 20)
        print("WARNING: CellTypist 모델이 강제로 'Cells_Lung_Airway.pkl'로 고정되었습니다!")
        print("자동 모델 선택이 비활성화되었습니다.")
        print("⚠️" * 20)

        forced_model = "Cells_Lung_Airway.pkl"
        try:
            model = models.Model.load(model=forced_model)
            print(f"Successfully loaded CellTypist model: {forced_model}")
        except Exception as e:
            print(f"Error loading model {forced_model}: {e}")
            model = models.Model.load(model="Immune_All_Low.pkl")
            print("Loaded fallback model: Immune_All_Low.pkl")
        self.model = model
        self.model_name = forced_model
        return model, forced_model

    def setup_directories(self):
        defaults = {"data_base": "/mnt/dks_nas2/data/", "output_base": "/mnt/dks_nas2/data/"}
        config = ConfigParser.parse_config_file(self.config_path)
        config.update(defaults)
        config["data_dir"] = _resolve_data_dir(config)
        return config

    def process_data(self, adata):
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)
        return adata

    def predict_cell_types(self, adata):
        if self.model is None:
            raise ValueError("Model not loaded. Call load_model_from_config() first.")
        predictions = celltypist.annotate(adata, model=self.model, majority_voting=True)
        return predictions.to_adata()

    def run_pipeline(self):
        self.load_model_from_config()
        print("=" * 60)
        print(f"🎯 Final selected model: {self.model_name}")
        print("=" * 60)

        config = self.setup_directories()
        data_dir = Path(config['data_dir'])
        print(f"📁 CellTypist data_dir: {data_dir}")

        h5ad_files = _find_h5ad_files(data_dir)
        if not h5ad_files:
            raise ValueError(
                f"No .h5ad files found under {data_dir}\n"
                "Looked for 'sce_qc.h5ad' or '*.h5ad' recursively."
            )
        print(f"🔎 Found {len(h5ad_files)} h5ad files (showing up to 5):")
        for p in h5ad_files[:5]:
            print("  -", p)

        aggregated = []
        for h5ad_path in h5ad_files:
            print(f"Processing {h5ad_path} ...")
            adata = sc.read(h5ad_path)
            adata = self.process_data(adata)
            adata_predict = self.predict_cell_types(adata)
            df = adata_predict.obs.copy()
            barcode = _get_barcode_series(df)
            # 컬럼 매핑(라벨/스코어 안전 처리)
            label_col = _pick_first(df.columns, ["predicted_labels", "majority_voting", "cell_type"])
            score_col = _pick_first(df.columns, ["confidence", "predicted_labels_score", "majority_voting_score", "predicted_proba"])
            if label_col is None:
                raise KeyError(f"CellTypist 결과 컬럼을 찾을 수 없습니다. 존재 컬럼: {list(df.columns)}")

            out = pd.DataFrame({
                "barcode": barcode.values,
                "predicted_labels": df[label_col].values,
                "confidence": (df[score_col].values if score_col else pd.NA),
            })
            
            aggregated.append(out[["barcode", "predicted_labels", "confidence"]])

        if not aggregated:
            raise ValueError("CellTypist: 수집된 결과가 없습니다(aggregated empty). h5ad 파일과 obs 컬럼을 확인하세요.")
        self.all_results = pd.concat(aggregated, ignore_index=True)
        return self.all_results

# =========================
# SCimilarity 파이프라인
# =========================
class SCimilarityPipeline:
    """SCimilarity를 사용한 세포 타입 예측 파이프라인"""

    def __init__(self, config_path, model_path):
        self.config_path = config_path
        self.model_path = model_path
        self.ca = None
        self.all_results = None

    def setup_directories(self):
        defaults = {"data_base": "/mnt/dks_nas2/data/", "output_base": "/mnt/dks_nas2/data/"}
        config = ConfigParser.parse_config_file(self.config_path)
        config.update(defaults)
        config["data_dir"] = _resolve_data_dir(config)
        print(f"📁 SCimilarity Data directory: {config['data_dir']}")
        return config

    def load_model(self):
        self.ca = CellAnnotation(model_path=self.model_path)
        print(f"SCimilarity 모델 로드 준비: {self.model_path}")
        return True

    def process_data(self, adata):
        adata = align_dataset(adata, self.ca.gene_order)
        if "counts" in adata.layers:
            adata = lognorm_counts(adata)
        return adata

    def get_embeddings_and_predictions(self, adata):
        adata.obsm["X_scimilarity"] = self.ca.get_embeddings(adata.X)
        predictions, nn_idxs, nn_dists, nn_stats = self.ca.get_predictions_knn(adata.obsm["X_scimilarity"])
        adata.obs["predictions_unconstrained"] = predictions.values
        adata.obs["vsAll_weighted"] = nn_stats["vsAll_weighted"].values
        return adata

    def run_pipeline(self):
        print(" SCimilarity 파이프라인 시작...")
        if not self.load_model():
            print("모델 로드 실패로 파이프라인 중단")
            return

        config = self.setup_directories()
        data_dir = Path(config["data_dir"])
        print(f"현재 디렉토리 : {data_dir}")

        h5ad_files = _find_h5ad_files(data_dir)
        if not h5ad_files:
            raise ValueError(
                f"No .h5ad files found under {data_dir}\n"
                "Looked for 'sce_qc.h5ad' or '*.h5ad' recursively."
            )
        print(f"🔎 Found {len(h5ad_files)} h5ad files (showing up to 5):")
        for p in h5ad_files[:5]:
            print("  -", p)

        aggregated = []
        for h5ad_path in h5ad_files:
            print(f"Processing {h5ad_path} ...")
            adata = sc.read(h5ad_path)
            adata = self.process_data(adata)
            adata = self.get_embeddings_and_predictions(adata)
            df = adata.obs.copy()
            barcode = _get_barcode_series(df)

            label_col = _pick_first(df.columns, ["predictions_unconstrained", "prediction", "annot"])
            score_col = _pick_first(df.columns, ["vsAll_weighted", "score", "conf"])
            if label_col is None:
                raise KeyError(f"SCimilarity 결과 컬럼을 찾을 수 없습니다. 존재 컬럼: {list(df.columns)}")

            out = pd.DataFrame({
                "barcode": barcode.values,
                "predictions_unconstrained": df[label_col].values,
                "vsAll_weighted": (df[score_col].values if score_col else pd.NA),
            })
            aggregated.append(out[["barcode", "predictions_unconstrained", "vsAll_weighted"]])

        if not aggregated:
            raise ValueError("SCimilarity: 수집된 결과가 없습니다(aggregated empty). h5ad 파일과 obs 컬럼을 확인하세요.")
        self.all_results = pd.concat(aggregated, ignore_index=True)
        return self.all_results

# =========================
# 통합 파이프라인
# =========================
class CombinedPipeline:
    """CellTypist와 SCimilarity를 함께 실행하는 통합 파이프라인"""

    def __init__(self, config_path, scimilarity_model_path):
        self.celltypist_pipeline = CellTypistPipeline(config_path)
        self.scimilarity_pipeline = SCimilarityPipeline(config_path, scimilarity_model_path)

    def run_celltypist_only(self):
        print("🔬 CellTypist 파이프라인만 실행...")
        return self.celltypist_pipeline.run_pipeline()

    def run_scimilarity_only(self):
        print("🧬 SCimilarity 파이프라인만 실행...")
        return self.scimilarity_pipeline.run_pipeline()

    def run_both_and_merge(self, tsv_input: str, merged_output: str):
        print("🔬🧬 CellTypist + SCimilarity 통합 파이프라인 실행...")
        print("\n" + "=" * 60)
        print("1  CellTypist 파이프라인 시작")
        print("=" * 60)
        ct_df = self.celltypist_pipeline.run_pipeline()

        print("\n" + "=" * 60)
        print("2  SCimilarity 파이프라인 시작")
        print("=" * 60)
        sc_df = self.scimilarity_pipeline.run_pipeline()
        
        print("\n" + "=" * 60)
        print("3  starting merge")
        print("=" * 60)
        merge_and_save_results(
            tsv_path=tsv_input,
            celltypist_results=ct_df,
            scimilarity_results=sc_df,
            output_path=merged_output
        )
        print("\n🎉 통합 파이프라인 완료!")

# =========================
# main
# =========================
def main():
    parser = argparse.ArgumentParser(
        description="Run CellTypist + SCimilarity pipelines and merge results with an existing TSV."
    )
    parser.add_argument(
        "input_tsv_path",
        help="기존 TSV 경로 (예: celltype_annotate_per_cell_extended.tsv)"
    )
    parser.add_argument(
        "output_tsv_path",
        nargs="?",
        default=None,
        help="출력 TSV 경로 (생략 시 입력 파일과 같은 폴더에 'merged_with_pipeline.tsv'로 저장)"
    )
    parser.add_argument(
        "-c", "--config",
        default="/home/seokwon/nas1_deep/filled_celltype/mvcell_parameter.csv",
        help="Config CSV 경로 (기본: /home/seokwon/nas1_deep/filled_celltype/mvcell_parameter.csv)"
    )
    parser.add_argument(
        "-m", "--scimilarity-model",
        default="/home/seokwon/nas1_deep/SCimilarity/model_v1.1",
        help="SCimilarity 모델 디렉토리 경로 (기본: /home/seokwon/nas1_deep/SCimilarity/model_v1.1)"
    )
    args = parser.parse_args()

    tsv_input = args.input_tsv_path
    if not Path(tsv_input).is_file():
        raise FileNotFoundError(f"Input file not found: {tsv_input}")

    # 출력 경로: 인자로 주어지면 사용, 없으면 기본값
    if args.output_tsv_path is not None:
        merged_output = args.output_tsv_path
    else:
        merged_output = str(Path(tsv_input).with_name("merged_with_pipeline.tsv"))

    config_file_path = args.config
    scimilarity_model_path = args.scimilarity_model

    print(f"[Config] config_file_path      : {config_file_path}")
    print(f"[Config] scimilarity_model_path: {scimilarity_model_path}")

    combined_pipeline = CombinedPipeline(
        config_path=config_file_path,
        scimilarity_model_path=scimilarity_model_path
    )

    combined_pipeline.run_both_and_merge(
        tsv_input=tsv_input,
        merged_output=merged_output
    )

    print("Pipeline completed successfully!")
    print(f"Input : {tsv_input}")
    print(f"Output: {merged_output}")


if __name__ == "__main__":
    main()

# %%
