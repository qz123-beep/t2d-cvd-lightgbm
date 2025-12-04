# utils.py
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import config

def savefig_pub(path_no_ext: str):
    """不加后缀的路径，自动存 pdf 与 png（期刊友好）"""
    Path(os.path.dirname(path_no_ext)).mkdir(parents=True, exist_ok=True)
    for ext in [".pdf", ".png"]:
        plt.savefig(path_no_ext + ext, bbox_inches="tight", dpi=600)
    plt.close()

def read_csv_any(p, encodings=("utf-8-sig","utf-8","gb18030","cp936","latin1")):
    p = str(p)
    for enc in encodings:
        try:
            return pd.read_csv(p, encoding=enc)
        except Exception:
            continue
    return pd.read_csv(p)

def enforce_unique_columns(df: pd.DataFrame) -> pd.DataFrame:
    """防止重复列名引发各种坑"""
    names = list(df.columns)
    used = set()
    for i, name in enumerate(names):
        if name not in used:
            used.add(name)
            continue
        j = 2
        new = f"{name}__{j}"
        while new in used:
            j += 1
            new = f"{name}__{j}"
        names[i] = new
        used.add(new)
    out = df.copy()
    out.columns = names
    return out

def _dedup_keep_order(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            out.append(x)
            seen.add(x)
    return out

def split_num_cat(X: pd.DataFrame, thresh_numeric: float = 0.7):
    """自动区分数值型和分类型列（与原脚本逻辑一致）"""
    Xc = X.copy()
    num_cols, cat_cols = [], []
    cols = list(Xc.columns)

    for i, col in enumerate(cols):
        s = Xc.iloc[:, i]
        if pd.api.types.is_numeric_dtype(s):
            num_cols.append(col)
        else:
            coerced = pd.to_numeric(s, errors="coerce")
            ratio = float(np.isfinite(coerced).mean())
            if ratio >= thresh_numeric:
                Xc.iloc[:, i] = coerced
                num_cols.append(col)
            else:
                cat_cols.append(col)

    num_cols = _dedup_keep_order(num_cols)
    cat_cols = [c for c in _dedup_keep_order(cat_cols) if c not in set(num_cols)]
    return Xc, num_cols, cat_cols

from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline

def make_ohe(drop_first=None):
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False, drop=drop_first)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=False, drop=drop_first)

def make_preprocessor(NUM_COLS, CAT_COLS, scale_numeric=True, drop_first_cat=False):
    num_pipe = Pipeline(steps=[
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler() if scale_numeric else "passthrough")
    ])
    cat_pipe = Pipeline(steps=[
        ("imputer", SimpleImputer(strategy="most_frequent")),
        ("ohe", make_ohe(drop_first="first" if drop_first_cat else None))
    ])
    pre = ColumnTransformer(
        transformers=[
            ("num", num_pipe, list(NUM_COLS)),
            ("cat", cat_pipe, list(CAT_COLS)),
        ],
        remainder="drop"
    )
    return pre
