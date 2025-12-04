# config.py
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

# ========= 基本路径（按你本地情况改） =========
CSV_PATH     = r"C:\zq\数据_SCI列名.csv"              # NHANES 开发集
EXT_CSV_PATH = r"C:\Users\Admin\Desktop\CVD统计表格.csv"  # 外部验证集（可留空）

OUT_DIR = Path(r"C:\zq\CVD\图4")   # 输出目录
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ========= 建模参数 =========
SEED      = 251
TEST_SIZE = 0.25
TARGET    = "CVD broad"

PREFERRED_MODEL = "LightGBM"
DELTA_AUROC     = 0.01

OPTUNA_TRIALS    = 120
RANDSEARCH_ITERS = 100

# 预设诊断阈值（5%、10%、20%）
PRESET_THRESH = [0.05, 0.10, 0.20]

# 图像尺寸 & 字体
FIG_W, FIG_H = 3.35, 3.35

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['axes.labelsize'] = 9
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.titlesize'] = 9

# ========= 候选特征 & 互斥约束 =========
RAW_CANDIDATES = [
    "Sex","Age","Height","Weight","BMI","Waist circumference","Hip circumference",
    "Smoking","Physical","Antihypertensive medication","Hypertension",
    "Fasting plasma glucose","HbA1c","Fasting insulin",
    "Total cholesterol","LDL-C","HDL-C","Triglycerides",
    "C-reactive protein","ALT","AST","GGT","Albumin","Serum creatinine","Blood cadmium",
    "Platelets","Neutrophils","Lymphocytes","SII","LSM","CAP",
    "HOMA-IR","TyG index","METS-IR","Triglyceride-to-HDL ratio","AIP",
    "FIB-4","AAR","HSI","NFS","APRI","eGFR","eGDR"
]

# 外部字段对齐
EXTERNAL_RENAME = {
    "Gender": "Sex",
    "HTN": "Hypertension",
    "Smoking status": "Smoking",
    "WC": "Waist circumference",
    "CVD": "CVD broad",
}

# 复合指标 ↔ 组成变量互斥（和文章里写的一致）
EXCLUSIVE_GROUPS = {
    "FIB-4": ["Age","AST","ALT","Platelets"],
    "AAR":   ["AST","ALT"],
    "eGDR":  ["Waist circumference","HbA1c","Hypertension"],
    "SII":   ["Neutrophils","Lymphocytes","Platelets"],
    "APRI":  ["AST","Platelets"],
}

# 临床基线模型变量
CLIN_BASE_VARS = ["Age","Sex","Hypertension","HbA1c"]
