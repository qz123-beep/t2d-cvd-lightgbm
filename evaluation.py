# evaluation.py
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import config
import utils

from sklearn.metrics import (
    roc_auc_score, average_precision_score, brier_score_loss,
    roc_curve, precision_recall_curve
)
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline

# ========= 校准相关 =========
def safe_logit(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))

def calib_slope_intercept(y, p):
    z = safe_logit(p).reshape(-1, 1)
    lr = LogisticRegression(solver="liblinear", C=1e6, max_iter=1000)
    lr.fit(z, y.astype(int))
    b = float(lr.coef_[0, 0])
    a = float(lr.intercept_[0])
    return a, b

def platt_recalibrate(p, a, b):
    z = safe_logit(p)
    return 1 / (1 + np.exp(-(a + b * z)))

def isotonic_recalibrate(y, p):
    ir = IsotonicRegression(out_of_bounds="clip")
    return ir.fit_transform(p, y)

# ========= ROC / PR / Calibration / DCA =========
def plot_roc_pr(y_true, y_proba, prefix):
    fpr, tpr, _ = roc_curve(y_true, y_proba)
    prec, rec, _ = precision_recall_curve(y_true, y_proba)

    # ROC
    plt.figure(figsize=(config.FIG_W, config.FIG_H))
    plt.plot(fpr, tpr, label=f"AUROC={roc_auc_score(y_true, y_proba):.3f}")
    plt.plot([0, 1], [0, 1], "--", color="gray", lw=1)
    plt.xlabel("1 - Specificity")
    plt.ylabel("Sensitivity")
    plt.title("ROC")
    plt.legend(loc="lower right")
    utils.savefig_pub(os.path.join(config.OUT_DIR, f"{prefix}_ROC"))

    # PR
    plt.figure(figsize=(config.FIG_W, config.FIG_H))
    ap = average_precision_score(y_true, y_proba)
    plt.plot(rec, prec, label=f"AUPRC={ap:.3f}")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Precision–Recall")
    base = y_true.mean()
    plt.hlines(base, 0, 1, linestyles="--", lw=1)
    plt.legend(loc="lower left")
    utils.savefig_pub(os.path.join(config.OUT_DIR, f"{prefix}_PR"))

def plot_calibration(y_true, y_proba, prefix, bins=10):
    prob_true, prob_pred = calibration_curve(
        y_true, y_proba, n_bins=bins, strategy="quantile"
    )
    plt.figure(figsize=(config.FIG_W, config.FIG_H))
    plt.plot([0, 1], [0, 1], "--", color="gray", lw=1, label="Perfect")
    plt.plot(prob_pred, prob_true, marker="o", lw=1.2, label="Model")
    plt.xlabel("Predicted probability")
    plt.ylabel("Observed probability")
    plt.title("Calibration plot")
    plt.legend(loc="upper left")
    utils.savefig_pub(os.path.join(config.OUT_DIR, f"{prefix}_Calibration"))

def youden_threshold(y_true, y_proba):
    fpr, tpr, thr = roc_curve(y_true, y_proba)
    j = tpr - fpr
    idx = int(np.argmax(j))
    return float(thr[idx])

def confusion_at(y_true, y_proba, thr):
    y_pred = (y_proba >= thr).astype(int)
    y_true = np.asarray(y_true).astype(int)

    TP = int(((y_true == 1) & (y_pred == 1)).sum())
    TN = int(((y_true == 0) & (y_pred == 0)).sum())
    FP = int(((y_true == 0) & (y_pred == 1)).sum())
    FN = int(((y_true == 1) & (y_pred == 0)).sum())

    sens = TP / (TP + FN) if (TP + FN) > 0 else np.nan
    spec = TN / (TN + FP) if (TN + FP) > 0 else np.nan
    ppv  = TP / (TP + FP) if (TP + FP) > 0 else np.nan
    npv  = TN / (TN + FN) if (TN + FN) > 0 else np.nan
    return dict(TP=TP, TN=TN, FP=FP, FN=FN,
                Sensitivity=sens, Specificity=spec,
                PPV=ppv, NPV=npv)

def threshold_table(y, p, name_prefix="Internal test"):
    rows = []
    youden = youden_threshold(y, p)
    for thr in config.PRESET_THRESH + [youden]:
        cm = confusion_at(y, p, thr)
        rows.append({"Set": name_prefix, "Threshold": thr, **cm})
    return pd.DataFrame(rows)

def plot_dca(y_true, y_proba, prefix, label="Model",
             pts=np.linspace(0.01, 0.99, 99),
             ann_pts=None):
    if ann_pts is None:
        ann_pts = config.PRESET_THRESH
    N = len(y_true)
    y_true = np.asarray(y_true).astype(int)
    base = y_true.mean()

    nb_model, nb_all, nb_none = [], [], []
    for pt in pts:
        thr = pt
        y_pred = (y_proba >= thr).astype(int)
        TP = ((y_true == 1) & (y_pred == 1)).sum()
        FP = ((y_true == 0) & (y_pred == 1)).sum()
        nb = (TP / N) - (FP / N) * (pt / (1 - pt))
        nb_model.append(nb)
        nb_all.append(base - (1 - base) * (pt / (1 - pt)))
        nb_none.append(0.0)

    plt.figure(figsize=(config.FIG_W, config.FIG_H))
    plt.plot(pts, nb_model, label=label)
    plt.plot(pts, nb_all, "--", label="Treat all")
    plt.plot(pts, nb_none, "--", label="Treat none")

    for a in ann_pts:
        ia = np.argmin(np.abs(pts - a))
        plt.scatter([pts[ia]], [nb_model[ia]], s=18)
        plt.text(pts[ia], nb_model[ia], f"  {a:.0%}", fontsize=8, va="bottom")

    plt.xlabel("Threshold probability")
    plt.ylabel("Net benefit")
    plt.title("Decision curve analysis")
    plt.legend(loc="best")
    utils.savefig_pub(os.path.join(config.OUT_DIR, f"{prefix}_DCA"))
