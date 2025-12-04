# interpretability.py
import os
import numpy as np
import matplotlib.pyplot as plt

import config
import utils

try:
    import shap
    HAS_SHAP = True
except Exception:
    HAS_SHAP = False


def ensure_dense_float(X):
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X, dtype=float)

def get_feature_names_from_pre(pre, input_cols=None):
    try:
        names = pre.get_feature_names_out()
        return [str(n).replace("num__","").replace("cat__","") for n in names]
    except Exception:
        pass

    names = []
    for name, trans, cols in pre.transformers_:
        if name == "remainder" and trans == "drop":
            continue
        col_list = list(cols) if isinstance(cols, (list, tuple, np.ndarray)) else [cols]
        if name == "num":
            names.extend([str(c) for c in col_list])
        elif name == "cat":
            ohe = None
            if hasattr(trans, "named_steps"):
                for key in trans.named_steps:
                    try:
                        from sklearn.preprocessing import OneHotEncoder
                        if isinstance(trans.named_steps[key], OneHotEncoder):
                            ohe = trans.named_steps[key]
                            break
                    except Exception:
                        pass
            if ohe is not None and hasattr(ohe, "get_feature_names_out"):
                out = ohe.get_feature_names_out(col_list)
                fixed = []
                for x in out:
                    s = str(x).replace("x0_", "")
                    if "_" in s:
                        var, lvl = s.split("_", 1)
                        fixed.append(f"{var}: {lvl.replace('_',' ')}")
                    else:
                        fixed.append(s)
                names.extend(fixed)
            else:
                names.extend([str(c) for c in col_list])
        else:
            names.extend([str(c) for c in col_list])
    return names

def shap_bar_group_by_parent(values_all, out_names, out_prefix, top_k=None):
    parents = []
    for s in out_names:
        p = s
        if s.startswith("num__"):
            p = s.split("num__", 1)[1]
        elif s.startswith("cat__"):
            p = s.split("cat__", 1)[1]
            if "_" in p:
                p = p[:p.rfind("_")]
        parents.append(p)

    df = (
        shap.utils.pd.DataFrame({
            "parent": parents,
            "abs_mean": np.nanmean(np.abs(values_all), axis=0),
        })
        if hasattr(shap, "utils") else
        __import__("pandas").DataFrame({
            "parent": parents,
            "abs_mean": np.nanmean(np.abs(values_all), axis=0),
        })
    )

    g = (
        df.groupby("parent", as_index=False)["abs_mean"]
        .sum()
        .sort_values("abs_mean", ascending=False)
    )
    if top_k is not None:
        g = g.head(top_k)

    plt.figure(figsize=(6, max(3, 0.3 * len(g) + 1)))
    plt.barh(g["parent"][::-1], g["abs_mean"][::-1])
    plt.xlabel("mean(|SHAP|) grouped by original feature")
    plt.title("SHAP – grouped bar")
    utils.savefig_pub(os.path.join(config.OUT_DIR, f"{out_prefix}_SHAP_bar_grouped"))

def shap_explain_any(clf, X_t, names, out_prefix,
                     max_display=20, sample_beeswarm=2000):
    if not HAS_SHAP:
        print("[提示] shap 未安装，跳过 SHAP 分析。")
        return

    import shap as _shap
    from sklearn.linear_model import LogisticRegression

    X_t = ensure_dense_float(X_t)
    n, m = X_t.shape

    names = list(names) if names is not None else [f"F{i}" for i in range(m)]
    if len(names) > m:
        names = names[:m]
    if len(names) < m:
        names = names + [f"F{i}" for i in range(len(names), m)]

    rng = np.random.RandomState(0)
    bg_idx = rng.choice(n, size=min(200, n), replace=False)
    background = X_t[bg_idx]

    def _build_explainer(model, background):
        try:
            return _shap.Explainer(model, background)
        except Exception:
            pass
        try:
            import lightgbm as _lgb
            if isinstance(model, _lgb.LGBMClassifier):
                try:
                    return _shap.TreeExplainer(
                        model,
                        feature_perturbation="interventional",
                        model_output="probability",
                    )
                except TypeError:
                    return _shap.TreeExplainer(model)
        except Exception:
            pass
        if isinstance(model, LogisticRegression):
            try:
                return _shap.LinearExplainer(model, background)
            except Exception:
                return _shap.explainers.Linear(model, background)
        if hasattr(model, "predict_proba"):
            f = lambda data: model.predict_proba(data)[:, 1]
        else:
            f = lambda data: model.decision_function(data)
        return _shap.KernelExplainer(f, background)

    def _to_nm(vals, X):
        vals = np.asarray(vals)
        if vals.ndim == 3:
            vals = vals[-1]
        nX, mX = X.shape
        if vals.shape == (mX, nX):
            vals = vals.T
        if vals.shape[0] == mX and vals.shape[1] == nX:
            vals = vals.T
        if vals.shape[1] > mX:
            vals = vals[:, :mX]
        elif vals.shape[1] < mX:
            pad = np.zeros((vals.shape[0], mX - vals.shape[1]))
            vals = np.concatenate([vals, pad], axis=1)
        return vals

    def _get_values(exp, X):
        try:
            sv = exp(X, check_additivity=False)
        except TypeError:
            sv = exp(X)
        vals = sv.values if hasattr(sv, "values") else sv
        return _to_nm(vals, X)

    explainer = _build_explainer(clf, background)

    # bar: 全数据
    values_all = _get_values(explainer, X_t)
    mean_abs = np.nanmean(np.abs(values_all), axis=0)
    order = np.argsort(-mean_abs)[:max_display]
    plt.figure(figsize=(6, max(3, 0.3 * len(order) + 1)))
    plt.barh([names[i] for i in order][::-1], mean_abs[order][::-1])
    plt.xlabel("mean(|SHAP value|)")
    plt.title("SHAP – bar (top features)")
    utils.savefig_pub(os.path.join(config.OUT_DIR, f"{out_prefix}_SHAP_bar"))

    # 额外：按原始变量聚合
    shap_bar_group_by_parent(values_all, names, out_prefix + "_main", top_k=10)

    # beeswarm: 子样本
    if n > sample_beeswarm:
        idx = rng.choice(n, size=sample_beeswarm, replace=False)
        X_bee = X_t[idx]
        vals_bee = _get_values(explainer, X_bee)
    else:
        X_bee = X_t
        vals_bee = values_all

    try:
        exp = shap.Explanation(vals_bee, data=X_bee, feature_names=names)
        shap.plots.beeswarm(exp, show=False, max_display=max_display)
    except Exception:
        shap.summary_plot(vals_bee, X_bee, feature_names=names,
                          show=False, max_display=max_display)
    plt.title("SHAP – beeswarm")
    utils.savefig_pub(os.path.join(config.OUT_DIR, f"{out_prefix}_SHAP_beeswarm"))

def shap_for_port(baseline_port, X_df, label_prefix):
    """给 BaselinePort 模型做 SHAP"""
    if not HAS_SHAP:
        print("[提示] shap 未安装，跳过 SHAP。")
        return
    pre = baseline_port.pre
    clf = baseline_port.clf
    X_t = pre.transform(X_df)
    X_t = ensure_dense_float(X_t)
    names = get_feature_names_from_pre(pre, input_cols=list(X_df.columns))
    shap_explain_any(clf, X_t, names, out_prefix=label_prefix,
                     max_display=30, sample_beeswarm=2000)
