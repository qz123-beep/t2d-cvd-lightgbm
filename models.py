# models.py
import numpy as np
import pandas as pd

import config
import utils

from sklearn.model_selection import (
    train_test_split, RepeatedStratifiedKFold, cross_validate
)
from sklearn.pipeline import Pipeline
from sklearn.metrics import (
    roc_auc_score, average_precision_score, brier_score_loss,
)
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.inspection import permutation_importance

# LightGBM / Optuna / shap 标志位
try:
    import lightgbm as lgb
    HAS_LGBM = True
except Exception:
    HAS_LGBM = False

try:
    import optuna
    HAS_OPTUNA = True
except Exception:
    HAS_OPTUNA = False


# ========= 模型集合 =========
def build_models():
    models = {}
    models["EN-Logistic"] = LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        l1_ratio=0.5,
        max_iter=4000,
        class_weight="balanced"
    )
    models["SVM-RBF"] = SVC(
        kernel="rbf",
        probability=True,
        class_weight="balanced",
        random_state=config.SEED
    )
    models["LDA"] = LinearDiscriminantAnalysis()
    models["KNN"] = KNeighborsClassifier(n_neighbors=31)
    models["GaussianNB"] = GaussianNB()
    if HAS_LGBM:
        models["LightGBM"] = lgb.LGBMClassifier(
            n_estimators=1000,
            learning_rate=0.03,
            num_leaves=31,
            subsample=0.8,
            colsample_bytree=0.8,
            min_child_samples=20,
            reg_alpha=0.0,
            reg_lambda=0.0,
            class_weight="balanced",
            random_state=config.SEED,
            deterministic=True,
            force_row_wise=True,
            bagging_seed=config.SEED,
            feature_fraction_seed=config.SEED,
            data_random_seed=config.SEED,
            verbosity=-1,
            n_jobs=1,
        )
    return models

def neg_brier_scorer(est, X, y):
    if hasattr(est, "predict_proba"):
        p = est.predict_proba(X)[:, 1]
    elif hasattr(est, "decision_function"):
        z = est.decision_function(X)
        z = (z - z.min()) / (z.max() - z.min() + 1e-9)
        p = z
    else:
        p = est.predict(X)
    return -brier_score_loss(y, p)

def cv_compare_models(X_train, y_train, preprocessor):
    """5×2 重复分层 CV 比较 6 个候选模型"""
    cv = RepeatedStratifiedKFold(
        n_splits=5, n_repeats=2, random_state=config.SEED
    )
    scorers = {
        "AUROC": "roc_auc",
        "AUPRC": "average_precision",
        "negBrier": neg_brier_scorer,
    }
    rows = []
    for name, base in build_models().items():
        pipe = Pipeline([("pre", preprocessor), ("clf", base)])
        print(f"[CV] {name} ...")
        res = cross_validate(
            pipe,
            X_train,
            y_train,
            cv=cv,
            scoring=scorers,
            n_jobs=-1,
            error_score="raise",
        )
        rows.append({
            "Model":       name,
            "AUROC_mean":  float(np.mean(res["test_AUROC"])),
            "AUROC_sd":    float(np.std(res["test_AUROC"])),
            "AUPRC_mean":  float(np.mean(res["test_AUPRC"])),
            "AUPRC_sd":    float(np.std(res["test_AUPRC"])),
            "Brier_mean":  float(np.mean(-res["test_negBrier"])),
            "Brier_sd":    float(np.std(-res["test_negBrier"])),
        })
    df = (
        pd.DataFrame(rows)
        .sort_values(["AUROC_mean", "AUPRC_mean"], ascending=False)
        .reset_index(drop=True)
    )
    return df

# ========= LightGBM 调参（Optuna / 随机搜索） =========
def fit_lgbm_with_params(Xtr, ytr, Xva, yva, params):
    clf = lgb.LGBMClassifier(**params)
    callbacks = []
    try:
        callbacks.append(lgb.early_stopping(100, verbose=False))
    except Exception:
        pass
    clf.fit(
        Xtr,
        ytr,
        eval_set=[(Xva, yva)],
        eval_metric="auc",
        callbacks=callbacks,
    )
    return clf

def tune_lgbm_optuna(X, y, pre, n_trials=config.OPTUNA_TRIALS, seed=config.SEED):
    if not (HAS_LGBM and HAS_OPTUNA):
        return None, None

    def objective(trial):
        params = dict(
            n_estimators   = trial.suggest_int("n_estimators", 800, 3000),
            learning_rate  = trial.suggest_float("learning_rate", 0.01, 0.08),
            num_leaves     = trial.suggest_int("num_leaves", 15, 63),
            min_child_samples = trial.suggest_int("min_child_samples", 10, 100),
            subsample      = trial.suggest_float("subsample", 0.6, 1.0),
            colsample_bytree = trial.suggest_float("colsample_bytree", 0.6, 1.0),
            reg_alpha      = trial.suggest_float("reg_alpha", 0.0, 5.0),
            reg_lambda     = trial.suggest_float("reg_lambda", 0.0, 5.0),
            class_weight   = "balanced",
            random_state   = seed + trial.number,
            deterministic  = True,
            force_row_wise = True,
            bagging_seed   = seed + trial.number,
            feature_fraction_seed = seed + trial.number,
            data_random_seed      = seed + trial.number,
            verbosity      = -1,
            n_jobs         = 1,
        )
        Xtr, Xva, ytr, yva = train_test_split(
            X, y, test_size=0.1, stratify=y, random_state=seed + trial.number
        )
        pre.fit(Xtr, ytr)
        Xtr_t = pre.transform(Xtr)
        Xva_t = pre.transform(Xva)
        if hasattr(Xtr_t, "toarray"):
            Xtr_t = Xtr_t.toarray()
        if hasattr(Xva_t, "toarray"):
            Xva_t = Xva_t.toarray()
        Xtr_t = np.asarray(Xtr_t, dtype=float)
        Xva_t = np.asarray(Xva_t, dtype=float)
        clf = fit_lgbm_with_params(Xtr_t, ytr, Xva_t, yva, params)
        proba = clf.predict_proba(Xva_t)[:, 1]
        return roc_auc_score(yva, proba)

    sampler = optuna.samplers.TPESampler(seed=seed)
    study = optuna.create_study(direction="maximize", sampler=sampler)
    study.optimize(objective, n_trials=n_trials, show_progress_bar=False)
    return study.best_params, study.best_value

def tune_lgbm_random(X, y, pre, n_iter=config.RANDSEARCH_ITERS, seed=config.SEED):
    if not HAS_LGBM:
        return None, None
    rng = np.random.RandomState(seed)
    best = -1.0
    best_params = None
    for i in range(n_iter):
        s = seed + i
        params = dict(
            n_estimators   = int(rng.randint(800, 3001)),
            learning_rate  = float(rng.uniform(0.01, 0.08)),
            num_leaves     = int(rng.randint(15, 64)),
            min_child_samples = int(rng.randint(10, 101)),
            subsample      = float(rng.uniform(0.6, 1.0)),
            colsample_bytree = float(rng.uniform(0.6, 1.0)),
            reg_alpha      = float(rng.uniform(0.0, 5.0)),
            reg_lambda     = float(rng.uniform(0.0, 5.0)),
            class_weight   = "balanced",
            random_state   = s,
            deterministic  = True,
            force_row_wise = True,
            bagging_seed   = s,
            feature_fraction_seed = s,
            data_random_seed      = s,
            verbosity      = -1,
            n_jobs         = 1,
        )
        Xtr, Xva, ytr, yva = train_test_split(
            X, y, test_size=0.1, stratify=y, random_state=s
        )
        pre.fit(Xtr, ytr)
        Xtr_t = pre.transform(Xtr)
        Xva_t = pre.transform(Xva)
        if hasattr(Xtr_t, "toarray"):
            Xtr_t = Xtr_t.toarray()
        if hasattr(Xva_t, "toarray"):
            Xva_t = Xva_t.toarray()
        Xtr_t = np.asarray(Xtr_t, dtype=float)
        Xva_t = np.asarray(Xva_t, dtype=float)
        clf = fit_lgbm_with_params(Xtr_t, ytr, Xva_t, yva, params)
        proba = clf.predict_proba(Xva_t)[:, 1]
        auc = roc_auc_score(yva, proba)
        if auc > best:
            best = auc
            best_params = params
    return best_params, best


# ========= 一个小包装：把预处理和 clf 打包，便于 SHAP / portable =========
class BaselinePort:
    def __init__(self, pre, clf):
        self.pre = pre
        self.clf = clf

    def predict_proba(self, X):
        Xt = self.pre.transform(X)
        if hasattr(Xt, "toarray"):
            Xt = Xt.toarray()
        Xt = np.asarray(Xt, dtype=float)
        if hasattr(self.clf, "predict_proba"):
            return self.clf.predict_proba(Xt)
        if hasattr(self.clf, "decision_function"):
            z = self.clf.decision_function(Xt)
            z = (z - z.min())/(z.max() - z.min() + 1e-9)
            return np.c_[1 - z, z]
        y = self.clf.predict(Xt)
        return np.c_[1 - y, y]

def fit_primary_model(Xtrain, ytrain, pre, chosen_name, best_params=None):
    """在训练集内部再划 0.1 做 early stopping 验证，最终返回 BaselinePort"""
    Xtr, Xva, ytr, yva = train_test_split(
        Xtrain, ytrain, test_size=0.1,
        stratify=ytrain, random_state=config.SEED
    )
    pre.fit(Xtr, ytr)
    Xtr_t = pre.transform(Xtr)
    Xva_t = pre.transform(Xva)
    if hasattr(Xtr_t, "toarray"):
        Xtr_t = Xtr_t.toarray()
    if hasattr(Xva_t, "toarray"):
        Xva_t = Xva_t.toarray()
    Xtr_t = np.asarray(Xtr_t, dtype=float)
    Xva_t = np.asarray(Xva_t, dtype=float)

    models_all = build_models()

    if chosen_name == "LightGBM" and HAS_LGBM:
        params = best_params or dict(
            n_estimators=1500, learning_rate=0.03, num_leaves=31,
            subsample=0.8, colsample_bytree=0.8, min_child_samples=20,
            reg_alpha=0.0, reg_lambda=0.0,
            class_weight="balanced",
            random_state=config.SEED, deterministic=True,
            force_row_wise=True,
            bagging_seed=config.SEED,
            feature_fraction_seed=config.SEED,
            data_random_seed=config.SEED,
            verbosity=-1, n_jobs=1
        )
        clf = fit_lgbm_with_params(Xtr_t, ytr, Xva_t, yva, params)
        return BaselinePort(pre, clf)

    clf = models_all[chosen_name]
    clf.fit(Xtr_t, ytr)
    return BaselinePort(pre, clf)

def pick_primary(cv_df: pd.DataFrame):
    """按你预设的非劣效规则挑最终算法"""
    df = cv_df.copy().reset_index(drop=True)
    top = df.iloc[0]
    if top["Model"] == config.PREFERRED_MODEL:
        print(f"[选择] {config.PREFERRED_MODEL} 已是最高模型。")
        return config.PREFERRED_MODEL
    if (df["Model"] == config.PREFERRED_MODEL).any():
        cand = df[df["Model"] == config.PREFERRED_MODEL].iloc[0]
        auroc_gap = float(top["AUROC_mean"] - cand["AUROC_mean"])
        auprc_gap = float(top["AUPRC_mean"] - cand["AUPRC_mean"])
        brier_better = bool(cand["Brier_mean"] <= top["Brier_mean"] + 1e-12)
        if (auroc_gap <= config.DELTA_AUROC) and (brier_better or (auprc_gap <= config.DELTA_AUROC)):
            print(f"[选择] 将 {config.PREFERRED_MODEL} 设为最终模型："
                  f"ΔAUROC {auroc_gap:.3f} ≤ {config.DELTA_AUROC} 且 Brier/AUPRC 不劣。")
            return config.PREFERRED_MODEL
    print(f"[选择] 采用最高模型：{top['Model']}")
    return top["Model"]

# ========= 置换重要性 & 互斥约束选紧凑模型 =========
def build_exclusion_map(groups):
    ex = {}
    for comp, parts in groups.items():
        ex.setdefault(comp, set()).update(parts)
        for p in parts:
            ex.setdefault(p, set()).add(comp)
    return ex

EXCLUSIONS = build_exclusion_map(config.EXCLUSIVE_GROUPS)

def permutation_stability(pipe, Xtrain_c, ytrain, n_runs=40, random_state=config.SEED):
    rng = np.random.RandomState(random_state)
    feats = list(Xtrain_c.columns)
    counts = {f: 0 for f in feats}
    for i in range(n_runs):
        res = permutation_importance(
            pipe,
            Xtrain_c,
            ytrain,
            scoring="roc_auc",
            n_repeats=1,
            random_state=int(rng.randint(0, 1_000_000_000)),
            n_jobs=1,
        )
        for j, f in enumerate(feats):
            if res.importances_mean[j] > 0:
                counts[f] += 1
    total = float(n_runs)
    stab = (
        pd.DataFrame({
            "feature": feats,
            "stability": [counts[f] / total for f in feats],
        })
        .sort_values("stability", ascending=False)
        .reset_index(drop=True)
    )
    return stab

def apply_constraints_ordered(sorted_feats, k, exclusions=EXCLUSIONS):
    chosen = []
    blocked = set()
    for f in sorted_feats:
        if f in blocked:
            continue
        chosen.append(f)
        for g in exclusions.get(f, []):
            blocked.add(g)
        if len(chosen) >= k:
            break
    return chosen

def cv_score_subset(Xtrain, ytrain, feat_list, chosen_name, base_params=None):
    Xs = Xtrain[feat_list]
    Xs_c, NUM, CAT = utils.split_num_cat(Xs)
    pre = utils.make_preprocessor(NUM, CAT, scale_numeric=True)
    base = fit_primary_model(Xs_c, ytrain, pre, chosen_name, best_params=base_params)

    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=config.SEED)
    au, ap, br = [], [], []
    for tr, va in cv.split(Xs_c, ytrain):
        Xtr, Xva = Xs_c.iloc[tr], Xs_c.iloc[va]
        ytr, yva = ytrain.iloc[tr], ytrain.iloc[va]
        pre.fit(Xtr, ytr)
        model = fit_primary_model(Xtr, ytr, pre, chosen_name, best_params=base_params)
        proba = model.predict_proba(Xva)[:, 1]
        au.append(roc_auc_score(yva, proba))
        ap.append(average_precision_score(yva, proba))
        br.append(brier_score_loss(yva, proba))
    return float(np.mean(au)), float(np.mean(ap)), float(np.mean(br))
