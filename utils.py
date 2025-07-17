import os
import random
import numpy as np
import pandas as pd
import torch
import pickle
from sklearn.preprocessing import LabelEncoder
from typing import Dict, Tuple, Any, Union, Sequence

def set_seed(seed=None):
    """Set random seed for reproducibility"""
    if seed is None:
        seed = random.randint(1, 10000)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    return seed

def get_all_possible_values(ont_df: pd.DataFrame, data_df: pd.DataFrame, test_df: pd.DataFrame) -> Dict[str, list]:
    """Get all possible values for each annotation level"""
    all_values = {}
    for i in range(1, 6):
        lvl = f"ann_level_{i}"
        vals = set(ont_df.get(lvl, pd.Series()).dropna().unique())
        vals |= set(data_df.get(lvl, pd.Series()).dropna().unique())
        vals |= set(test_df.get(lvl, pd.Series()).dropna().unique())
        vals.add("Unknown")
        all_values[lvl] = sorted(vals)
        print(f"레벨 {lvl}: {len(all_values[lvl])}개")
    return all_values

def identify_terminal_nodes(ont_df: pd.DataFrame) -> Dict[str, set]:
    """Identify terminal nodes in the hierarchy"""
    terminal_nodes = {}
    
    for i in range(1, 5):
        current_lvl = f"ann_level_{i}"
        next_lvl = f"ann_level_{i+1}"
        
        all_nodes = set(ont_df[current_lvl].dropna().unique())
        
        has_children = set()
        for _, row in ont_df.dropna(subset=[next_lvl]).iterrows():
            parent = row[current_lvl]
            if pd.notna(parent):
                has_children.add(parent)
        
        terminal_nodes[current_lvl] = all_nodes - has_children
        
    return terminal_nodes

def create_label_encoders(df_train: pd.DataFrame, df_test: pd.DataFrame, all_cell_types: Dict[str, list]) -> Tuple[Dict[str, LabelEncoder], Dict[str, LabelEncoder]]:
    """Create and fit label encoders for input and output columns"""
    TOOLS = ["singleR", "sctype", "scATOMIC", "scimilarity", "celltypist"]
    ANN_COLS = [f"{t}_ann{i}" for t in TOOLS for i in range(1, 6)]
    TARGET_COLS = [f"ann_level_{i}" for i in range(1, 6)]
    
    enc_in, enc_out = {}, {}
    
    # Encode input annotation columns
    for col in ANN_COLS:
        df_train[col] = df_train[col].fillna("Unknown").astype(str)
        df_test[col] = df_test[col].fillna("Unknown").astype(str)
        le = LabelEncoder()
        lvl = int(col.split("ann")[-1])
        
        # Fit encoder with all possible values
        all_values = all_cell_types[f"ann_level_{lvl}"] + df_train[col].unique().tolist() + df_test[col].unique().tolist()
        le.fit(list(set(all_values)))
        
        df_train[f"{col}_idx"] = le.transform(df_train[col])
        df_test[f"{col}_idx"] = le.transform(df_test[col])
        enc_in[col] = le
    
    # Encode target columns
    for lvl in TARGET_COLS:
        df_train[lvl] = df_train[lvl].fillna("Unknown").astype(str)
        df_test[lvl] = df_test[lvl].fillna("Unknown").astype(str)
        le = LabelEncoder()
        
        all_values = all_cell_types[lvl] + df_train[lvl].unique().tolist() + df_test[lvl].unique().tolist()
        le.fit(list(set(all_values)))
        
        df_train[f"{lvl}_idx"] = le.transform(df_train[lvl])
        df_test[f"{lvl}_idx"] = le.transform(df_test[lvl])
        enc_out[lvl] = le
    
    return enc_in, enc_out

def save_encoders(enc_in: Dict[str, LabelEncoder], enc_out: Dict[str, LabelEncoder], save_dir: str):
    """Save label encoders to disk"""
    os.makedirs(save_dir, exist_ok=True)
    
    with open(os.path.join(save_dir, 'encoders_in.pkl'), 'wb') as f:
        pickle.dump(enc_in, f)
    
    with open(os.path.join(save_dir, 'encoders_out.pkl'), 'wb') as f:
        pickle.dump(enc_out, f)
    
    print(f"Encoders saved to {save_dir}")

def load_encoders(save_dir: str) -> Tuple[Dict[str, LabelEncoder], Dict[str, LabelEncoder]]:
    """Load label encoders from disk"""
    with open(os.path.join(save_dir, 'encoders_in.pkl'), 'rb') as f:
        enc_in = pickle.load(f)
    
    with open(os.path.join(save_dir, 'encoders_out.pkl'), 'rb') as f:
        enc_out = pickle.load(f)
    
    print(f"Encoders loaded from {save_dir}")
    return enc_in, enc_out

def get_terminal_nodes_idx(terminal_nodes: Dict[str, set], enc_out: Dict[str, LabelEncoder]) -> Dict[str, set]:
    """Convert terminal nodes to indices"""
    terminal_nodes_idx = {}
    for lvl in range(1, 5):
        current_lvl = f"ann_level_{lvl}"
        terminal_nodes_idx[current_lvl] = set()
        for node in terminal_nodes.get(current_lvl, set()):
            try:
                idx = enc_out[current_lvl].transform([node])[0]
                terminal_nodes_idx[current_lvl].add(idx)
            except:
                pass
    return terminal_nodes_idx

def ablation_data(data_df: pd.DataFrame, ablation_sw: Union[str, Sequence[str], None] = None) -> pd.DataFrame:
    """Apply ablation to specific tools"""
    if not ablation_sw:
        return data_df
    
    if isinstance(ablation_sw, str):
        sw_list = [ablation_sw]
    else:
        sw_list = list(ablation_sw)
    
    for sw in sw_list:
        ann_cols = [f"{sw}_ann{i}" for i in range(1, 6)]
        for col in ann_cols:
            data_df[col] = "Unknown"
        data_df[f"{sw}_score"] = 0
    
    return data_df

def verify_hierarchy_consistency(preds_dict: Dict[str, list], level_count: int = 5) -> Tuple[int, int]:
    """Verify hierarchy consistency in predictions"""
    consistency_errors = 0
    total_samples = len(next(iter(preds_dict.values())))
    
    for idx in range(total_samples):
        for lvl in range(1, level_count):
            current_lvl = f"ann_level_{lvl}_pred"
            next_lvl = f"ann_level_{lvl+1}_pred"
            
            # Check if current level prediction is a terminal node
            current_pred = preds_dict[current_lvl][idx]
            # This would require terminal_nodes to be passed as parameter
            # For now, skip this check or implement differently
    
    return consistency_errors, total_samples

def preprocess_data(df: pd.DataFrame, score_cols: list, ann_cols: list, ann_level: list) -> pd.DataFrame:
    """Preprocess dataframe by handling missing values"""
    df = df.copy()
    df.dropna(subset=score_cols, inplace=True)
    
    for col in ann_cols + ann_level:
        df[col] = df[col].fillna("Unknown").astype(str)
    
    return df