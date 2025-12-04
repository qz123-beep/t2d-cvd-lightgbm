# run_pipeline.py
import os
from pathlib import Path
import numpy as np
import pandas as pd

import config
import utils
import models
import evaluation
import interpretability

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline

def main():
    print(f"[info] 输出目录: {config.OUT_DIR}")
    # ===== 1. 读 NHANES 开发集 =====
    df = utils.read_csv_any(config.CSV_PATH)
    if config.TARGET not in df.columns:
        raise ValueError(f"Target '{config.TARGET}' not in columns.")
    df = utils.enforce_unique_columns(df)

    present = utils._dedup_keep_order(
        [c for c in config.RAW_CANDIDATES if c in df.columns]
    )
    y = pd.to_numeric(df[config.TARGET], errors="coerce").astype(int)
    X = df[present].copy()

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=config.TEST_SIZE,
        stratify=y, random_state=config.SEED
    )
    print(f"[info] Train/Test shape: {X_train.shape} / {X_test.shape}")
    print(f"[info] Train prevalence: {y_train.mean():.3f} | Test: {y_test.mean():.3f}")

    # ===== 2. 预处理 + 6 模型 CV 比较 =====
    X_train_c, NUM, CAT = utils.split_num_cat(X_train)
    pre = utils.make_preprocessor(NUM, CAT, scale_numeric=True)

    cv_df = models.cv_compare_models(X_train_c, y_train, pre)
    cv_path = config.OUT_DIR / "S1_CV_performance.csv"
    cv_df.to_csv(cv_path, index=False, encoding="utf-8-sig")
    print("[saved]", cv_path)

    primary_name = models.pick_primary(cv_df)
    best_params = None
    if primary_name == "LightGBM" and models.HAS_LGBM:
        if models.HAS_OPTUNA:
            best_params, _ = models.tune_lgbm_optuna(X_train_c, y_train, pre)
        else:
            best_params, _ = models.tune_lgbm_random(X_train_c, y_train, pre)

    # ===== 3. 训练基线模型，内部测试评估 =====
    baseline = models.fit_primary_model(X_train_c, y_train, pre,
                                        primary_name, best_params)
    y_test_proba = baseline.predict_proba(X_test)[:, 1]
    auroc = roc_auc_score_safe(y_test, y_test_proba)
    auprc = evaluation.average_precision_score(y_test, y_test_proba) \
        if hasattr(evaluation, "average_precision_score") else None
    brier = brier_score_loss_safe(y_test, y_test_proba)
    a, b = evaluation.calib_slope_intercept(y_test, y_test_proba)
    print(f"[Test] AUROC={auroc:.3f}  Brier={brier:.3f}  Calib=({a:.3f},{b:.3f})")

    evaluation.plot_roc_pr(y_test, y_test_proba, prefix="test_baseline")
    evaluation.plot_calibration(y_test, y_test_proba, prefix="test_baseline")
    evaluation.plot_dca(y_test, y_test_proba, prefix="test_baseline",
                        label=primary_name)

    thr_df_internal = evaluation.threshold_table(
        y_test.values, y_test_proba, name_prefix="Internal test"
    )

    # ===== 4. SHAP – 基线模型 =====
    try:
        interpretability.shap_for_port(baseline, X_train, "train_full")
        interpretability.shap_for_port(baseline, X_test,  "test_full")
    except Exception as e:
        print(f"[提示] 基线模型 SHAP 失败：{e}")

    # ===== 5. 置换重要性 & 紧凑 K=6/7/8 模型 =====
    pipe_full = Pipeline([("pre", pre), ("clf", baseline.clf)])
    pipe_full.fit(X_train_c, y_train)
    stab = models.permutation_stability(
        pipe_full, X_train_c, y_train,
        n_runs=40, random_state=config.SEED
    )
    stab_path = config.OUT_DIR / "S2_permutation_stability.csv"
    stab.to_csv(stab_path, index=False, encoding="utf-8-sig")
    print("[saved]", stab_path)

    ordered = list(stab["feature"])
    compact_rows = []
    bestK, best_feats, best_tuple = None, None, None
    for K in [6, 7, 8]:
        feats_k = models.apply_constraints_ordered(
            ordered, k=K, exclusions=models.EXCLUSIONS
        )
        au, ap, br = models.cv_score_subset(
            X_train, y_train, feats_k,
            chosen_name=primary_name, base_params=best_params
        )
        compact_rows.append({
            "K": K,
            "CV_AUROC": au,
            "CV_AUPRC": ap,
            "CV_Brier": br,
            "Features": "; ".join(feats_k),
        })
        if (best_tuple is None) or (au > best_tuple[0]):
            best_tuple = (au, ap, br)
            bestK = K
            best_feats = feats_k

    pd.DataFrame(compact_rows).to_csv(
        config.OUT_DIR / "S2b_compactK_CV.csv",
        index=False, encoding="utf-8-sig"
    )

    # 紧凑模型训练 + 测试
    X_train_comp = X_train[best_feats]
    X_train_comp_c, NUMc, CATc = utils.split_num_cat(X_train_comp)
    pre_c = utils.make_preprocessor(NUMc, CATc, scale_numeric=True)
    compact = models.fit_primary_model(X_train_comp_c, y_train,
                                       pre_c, primary_name, best_params)
    y_test_comp = compact.predict_proba(X_test[best_feats])[:, 1]
    evaluation.plot_roc_pr(y_test, y_test_comp,
                           prefix=f"test_compact_K{bestK}")
    evaluation.plot_calibration(y_test, y_test_comp,
                                prefix=f"test_compact_K{bestK}")
    evaluation.plot_dca(y_test, y_test_comp,
                        prefix=f"test_compact_K{bestK}",
                        label=f"Compact K={bestK}")

    try:
        interpretability.shap_for_port(
            compact, X_train[best_feats],
            f"train_compact_K{bestK}"
        )
        interpretability.shap_for_port(
            compact, X_test[best_feats],
            f"test_compact_K{bestK}"
        )
    except Exception as e:
        print(f"[提示] 紧凑模型 SHAP 失败：{e}")

    # ===== 6. 外部可迁移模型（portable） =====
    ext_thr_df = None
    if config.EXT_CSV_PATH and Path(config.EXT_CSV_PATH).exists():
        ext = utils.read_csv_any(config.EXT_CSV_PATH).rename(columns=config.EXTERNAL_RENAME)
        ext = utils.enforce_unique_columns(ext)
        if config.TARGET in ext.columns:
            y_ext = pd.to_numeric(ext[config.TARGET], errors="coerce").astype(int)
            inter = [c for c in X_train.columns if c in ext.columns]
            if len(inter) >= 3:
                Xtr_port = X_train[inter]
                Xtr_port_c, NUMp, CATp = utils.split_num_cat(Xtr_port)
                pre_p = utils.make_preprocessor(NUMp, CATp, scale_numeric=True)
                port = models.fit_primary_model(Xtr_port_c, y_train,
                                                pre_p, primary_name, best_params)
                X_ext_use = ext[inter].copy()
                y_ext_proba = port.predict_proba(X_ext_use)[:, 1]

                evaluation.plot_roc_pr(y_ext, y_ext_proba,
                                       prefix="external_portable")
                evaluation.plot_calibration(y_ext, y_ext_proba,
                                            prefix="external_portable")
                evaluation.plot_dca(y_ext, y_ext_proba,
                                    prefix="external_portable",
                                    label="Portable")
                ext_thr_df = evaluation.threshold_table(
                    y_ext.values, y_ext_proba,
                    name_prefix="External portable"
                )

                try:
                    interpretability.shap_for_port(port, Xtr_port,
                                                   "portable_train")
                    interpretability.shap_for_port(port, X_ext_use,
                                                   "portable_external")
                except Exception as e:
                    print(f"[提示] 外部 SHAP 失败：{e}")

    # ===== 7. 阈值表合并输出 =====
    if ext_thr_df is not None:
        thr_all = pd.concat([thr_df_internal, ext_thr_df],
                            axis=0, ignore_index=True)
    else:
        thr_all = thr_df_internal
    thr_all.to_csv(
        config.OUT_DIR / "S3_thresholds.csv",
        index=False, encoding="utf-8-sig"
    )

    print("\n[done] 所有图和表已输出至：", config.OUT_DIR)


# 小工具：安全算 AUROC / Brier（防止偶尔崩）
from sklearn.metrics import roc_auc_score, brier_score_loss
def roc_auc_score_safe(y, p):
    try:
        return roc_auc_score(y, p)
    except Exception:
        return np.nan

def brier_score_loss_safe(y, p):
    try:
        return brier_score_loss(y, p)
    except Exception:
        return np.nan

if __name__ == "__main__":
    main()
