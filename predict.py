#!/usr/bin/env python3
"""
Bidirectional Hierarchical Transformer Prediction Script
"""

import os
import argparse
import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, 
    f1_score, matthews_corrcoef
)

# Import custom modules
from utils import (
    set_seed, get_all_possible_values, identify_terminal_nodes,
    load_encoders, get_terminal_nodes_idx, preprocess_data
)
from model import MultiLevelDataset, SimpleBidirectionalTransformer

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Predict with Bidirectional Transformer')
    parser.add_argument('--model-dir', type=str, required=True, 
                       help='Directory containing trained model and encoders')
    parser.add_argument('--data-dir', type=str, default='.', 
                       help='Data directory')
    parser.add_argument('--test-file', type=str, 
                       default='test_data.tsv',
                       help='Test data file')
    parser.add_argument('--ont-file', type=str,
                       default='HLCA_tree.csv',
                       help='Ontology file')
    parser.add_argument('--batch-size', type=int, default=1024,
                       help='Batch size for prediction')
    parser.add_argument('--seed', type=int, default=9401,
                       help='Random seed (should match training)')
    parser.add_argument('--output-dir', type=str, default='./results',
                       help='Output directory (default: ./results)')
    
    return parser.parse_args()

def load_model_and_encoders(model_dir, device, ont_df):
    """Load trained model and encoders"""
    print(f"Loading model and encoders from {model_dir}")
    
    # Load encoders
    enc_in, enc_out = load_encoders(model_dir)
    
    # Constants (should match training)
    TOOLS = ["singleR", "sctype", "scATOMIC", "scimilarity", "celltypist"]
    SCORE_COLS = [f"{t}_score" for t in TOOLS]
    
    # Load model configuration (you might want to save this during training)
    # For now, using default values that match training
    d_model = 256
    nhead = 8
    num_layers = 2
    dropout = 0.5
    
    # Get terminal nodes
    terminal_nodes = identify_terminal_nodes(ont_df)
    terminal_nodes_idx = get_terminal_nodes_idx(terminal_nodes, enc_out)
    
    # Initialize model
    model = SimpleBidirectionalTransformer(
        enc_in, enc_out, ont_df, len(SCORE_COLS), terminal_nodes_idx,
        d_model=d_model, nhead=nhead, num_layers=num_layers,
        dropout=dropout, device=device
    ).to(device)
    
    # Load model weights
    model_path = None
    for f in os.listdir(model_dir):
        if f.startswith("best_multilevel_seed") and f.endswith(".pth"):
            model_path = os.path.join(model_dir, f)
            break
    
    if model_path is None:
        raise FileNotFoundError("No model file found in model directory")
    
    print(f"Loading model weights from {model_path}")
    checkpoint = torch.load(model_path, map_location=device)
    if any(key.startswith('_orig_mod.') for key in checkpoint.keys()):
        # strip torch.compile prefix
        checkpoint = { (k[10:] if k.startswith('_orig_mod.') else k): v for k, v in checkpoint.items() }

    model.load_state_dict(checkpoint)
    model.eval()
    
    return model, enc_in, enc_out, terminal_nodes_idx

def get_model_predictions(model, loader, device, enc_in, target_cols):
    """Get model predictions"""
    all_preds = {lvl: [] for lvl in target_cols}
    
    with torch.no_grad():
        for cat, scores, targets in loader:
            cat = cat.to(device)
            scores = scores.to(device)
            
            forward_outputs, _ = model(cat, scores, enc_in=enc_in)
            
            for i, lvl in enumerate(target_cols):
                preds = forward_outputs[lvl].argmax(1).cpu().numpy()
                all_preds[lvl].extend(preds)
    
    return all_preds

def save_prediction_results(model_preds, output_dir,barcodes=None):
    """Save prediction results to CSV (predictions only)"""
    # Create predictions dataframe
    predictions_df = pd.DataFrame({
        f"{lvl}_prediction": model_preds[lvl] 
        for lvl in model_preds.keys()
    })
    if barcodes is not None:
        if len(barcodes) != len(predictions_df):
            raise ValueError(f"Barcode length ({len(barcodes)}) does not match predictions ({len(predictions_df)})")
        predictions_df.insert(0, "barcode", barcodes)
    
    # Save results 
    pred_path = os.path.join(output_dir, "cell_type_predictions.csv")
    predictions_df.to_csv(pred_path, index=False)
    
    return pred_path, predictions_df

def print_prediction_summary(model_preds, total_cells, output_dir):
    """Print summary of predictions and save report"""
    # Terminal output
    print(f"\n🎯 Prediction Results:")
    print(f"📊 Total cells predicted: {total_cells:,}")
    
    # Prepare report content
    report_lines = []
    report_lines.append("Cell Type Annotation Prediction Report")
    report_lines.append("=" * 50)
    report_lines.append(f"Total cells predicted: {total_cells:,}")
    report_lines.append("")
    
    for lvl in model_preds.keys():
        pred_counts = pd.Series(model_preds[lvl]).value_counts()
        
        # Terminal output
        print(f"\n{lvl} predictions (Top 5):")
        
        # Report content
        report_lines.append(f"{lvl} predictions (Top 5):")
        report_lines.append("-" * 30)
        
        for cell_type, count in pred_counts.head(5).items():
            percentage = (count / total_cells) * 100
            line = f"  {cell_type}: {count:,} cells ({percentage:.1f}%)"
            print(line)  # Terminal
            report_lines.append(line)  # Report
        
        report_lines.append("")
    
    # Save report to file
    report_path = os.path.join(output_dir, "prediction_summary_report.txt")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"\n📄 Summary report saved to: {report_path}")
    
    return report_path

def _read_table_auto(path: str) -> pd.DataFrame:            # NEW
    # 탭/콤마 자동 판별
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception:
        # 마지막 폴백: 탭
        return pd.read_csv(path, sep="\t")

def _ensure_inputs(df: pd.DataFrame, tools: list[str]) -> pd.DataFrame:   # NEW
    """
    모델 입력에 필요한 컬럼(툴별 ann1~5, score)을 보장.
    없으면 ann은 'Unknown', score는 0.0 으로 생성.
    """
    out = df.copy()
    for t in tools:
        for i in range(1, 6):
            c = f"{t}_ann{i}"
            if c not in out.columns:
                out[c] = "Unknown"
        s = f"{t}_score"
        if s not in out.columns:
            out[s] = 0.0
    return out

def _encode_inputs(df: pd.DataFrame, enc_in, tools: list[str]) -> pd.DataFrame:  # NEW
    """
    enc_in으로 툴별 ann1~5를 인덱스로 변환. Unknown 처리.
    """
    out = df.copy()
    for t in tools:
        for i in range(1, 6):
            col = f"{t}_ann{i}"
            out[col] = out[col].fillna("Unknown").astype(str)
            try:
                out[f"{col}_idx"] = enc_in[col].transform(out[col])
            except ValueError:
                # enc_in에 없는 범주는 Unknown으로 대체
                classes = set(enc_in[col].classes_)
                out[col] = out[col].apply(lambda x: x if x in classes else "Unknown")
                out[f"{col}_idx"] = enc_in[col].transform(out[col])
    # score 숫자화
    for t in tools:
        s = f"{t}_score"
        out[s] = pd.to_numeric(out[s], errors="coerce").fillna(0.0)
    return out

    
def main():
    args = parse_args()
    
    # Setup
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    set_seed(args.seed)
    
    # Constants
    TOOLS = ["singleR", "sctype", "scATOMIC", "scimilarity", "celltypist"]
    SCORE_COLS = [f"{t}_score" for t in TOOLS]
    ANN_COLS = [f"{t}_ann{i}" for t in TOOLS for i in range(1, 6)]
    ANN_LEVEL = [f"ann_level_{i}" for i in range(1, 6)]
    TARGET_COLS = [f"ann_level_{i}" for i in range(1, 6)]
    level_cols = TARGET_COLS
    
    # Output directory
    output_dir = args.output_dir if args.output_dir else "./results"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    print("Loading data...")
    ont_df = pd.read_csv(os.path.join(args.data_dir, args.ont_file))
    df_eval = pd.read_csv(os.path.join(args.data_dir, args.test_file), sep='\t')
    
    print(f"Input data: {df_eval.shape}")
    
    barcodes = df_eval["barcode"].astype(str).tolist()
        
    has_gt = all(c in df_eval.columns for c in TARGET_COLS)
    
    # Preprocess evaluation data
    # df_eval = preprocess_data(df_eval, SCORE_COLS, ANN_COLS, ANN_LEVEL)
    
    # Load model and encoders
    model, enc_in, enc_out, terminal_nodes_idx = load_model_and_encoders(
        args.model_dir, device, ont_df
    )
    
    df_eval = _ensure_inputs(df_eval, TOOLS)
    df_eval = _encode_inputs(df_eval, enc_in, TOOLS)

    if has_gt:
        for lvl in TARGET_COLS:
            df_eval[lvl] = df_eval[lvl].fillna("Unknown").astype(str)
            try:
                df_eval[f"{lvl}_idx"] = enc_out[lvl].transform(df_eval[lvl])
            except ValueError:
                classes = set(enc_out[lvl].classes_)
                df_eval[lvl] = df_eval[lvl].apply(lambda x: x if x in classes else "Unknown")
                df_eval[f"{lvl}_idx"] = enc_out[lvl].transform(df_eval[lvl])
        y_dummy = df_eval[[f"{lvl}_idx" for lvl in TARGET_COLS]].to_numpy()
    else:
        y_dummy = np.zeros((len(df_eval), 5), dtype=int)

    # Encode evaluation data
    print("Encoding input data...")

    if has_gt:
        for lvl in TARGET_COLS:
            df_eval[lvl] = df_eval[lvl].fillna("Unknown").astype(str)
            try:
                df_eval[f"{lvl}_idx"] = enc_out[lvl].transform(df_eval[lvl])
            except ValueError:
                classes = set(enc_out[lvl].classes_)
                df_eval[lvl] = df_eval[lvl].apply(lambda x: x if x in classes else "Unknown")
                df_eval[f"{lvl}_idx"] = enc_out[lvl].transform(df_eval[lvl])
        y_dummy = df_eval[[f"{lvl}_idx" for lvl in TARGET_COLS]].to_numpy()
    else:
        y_dummy = np.zeros((len(df_eval), 5), dtype=int)

    # Prepare data
    feature_cols = [f"{c}_idx" for c in ANN_COLS] + SCORE_COLS
    eval_X = df_eval[feature_cols]
    
    # Create dataset and dataloader
    eval_ds = MultiLevelDataset(eval_X, y_dummy, ANN_COLS, SCORE_COLS)
    eval_loader = DataLoader(eval_ds, batch_size=args.batch_size, shuffle=False, 
                            num_workers=4, pin_memory=True if device.type == 'cuda' else False)
    
    # Get model predictions
    print("Running predictions...")
    eval_preds_idx = get_model_predictions(model, eval_loader, device, enc_in, TARGET_COLS)
    
    # Convert to string labels
    model_preds = {
        lvl: enc_out[lvl].inverse_transform(eval_preds_idx[lvl])
        for lvl in level_cols
    }
    
    # Save prediction results
    pred_path, result_df = save_prediction_results(model_preds, output_dir, barcodes=barcodes)
    print(f"\n💾 Predictions saved to: {pred_path}")
    
    # Print prediction summary
    args.output_dir = os.path.abspath(args.output_dir)
    print_prediction_summary(model_preds, len(df_eval), args.output_dir)
    
    print(f"\n✅ Prediction completed. Results saved to {output_dir}")

if __name__ == "__main__":
    main()
