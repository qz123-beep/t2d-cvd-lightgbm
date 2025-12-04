# =========================================================
# 03_make_S1_variable_mapping_from_raw.R
#
# 功能：
#   直接从原始 NHANES XPT 文件构建“最终分析队列（DM 且 30–75 岁，
#   剔除肝病/丙肝/肾衰/恶性肿瘤、重度饮酒等）”，
#   在 *未插补* 的前提下：
#     - 规范化变量命名（与建模脚本一致）
#     - 计算各变量的分布、范围、原始缺失率
#     - 生成 Supplementary Table S1：
#       变量映射 + 取值范围 + NHANES 缺失率（原始）
#   可选：如果已经准备好了外部验证集，也可以同时计算外部队列分布和缺失率。
#
# 注意：
#   这里不会做 na.aggregate 或 mice 插补，
#   因此 Missingness_NHANES_percent 代表的是“原始缺失”。
# =========================================================

rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(foreign)  # read.xport
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
})

## 全局开关：是否把“被建议阿司匹林 (RXQ510==1)”纳入 broad？
incl_aspirin_in_broad <- FALSE  # 如需纳入，改为 TRUE

## 工具函数 ----
chk <- function(df, tag) { cat(sprintf("[%s] n = %d\n", tag, nrow(df))); df }
to01 <- function(v) as.integer(v == 1)
miss_rate <- function(x) round(mean(is.na(x)) * 100, 1)

# ---------- 0. 外部验证集路径（可先留空） -----------------

# 若暂时没有外部验证集，就留空字符串 "" 即可
path_external <- ""   # 例如 "C:/zq/外部验证/DM_CVD_external.csv"

# 若外部验证集列名与 NHANES 不同，可以在后面 external_name_map 中指定映射


# =========================================================
# 1) 按照原始主脚本路径：逐块读取 NHANES 各模块
#    —— 完全不做插补，只做重编码
# =========================================================

## ========= 1) 空腹血糖 =========
P_GLU <- read.xport("C:/zq/Nhanes/空腹血糖/P_GLU.xpt")
GLU_J <- read.xport("C:/zq/Nhanes/空腹血糖/GLU_J.xpt")
glucose <- rbind(
  P_GLU[, c("SEQN","LBDGLUSI","LBXGLU")],
  GLU_J[, c("SEQN","LBDGLUSI","LBXGLU")]
) %>% rename(FPG1 = LBDGLUSI,   # mmol/L
             FPG2 = LBXGLU)     # mg/dL

## ========= 2) 糖尿病问卷 =========
DIQ_L <- read.xport("C:/zq/Nhanes/糖尿病问卷/DIQ_L.xpt")
P_DIQ <- read.xport("C:/zq/Nhanes/糖尿病问卷/P_DIQ.xpt")
DIQ_J <- read.xport("C:/zq/Nhanes/糖尿病问卷/DIQ_J.xpt")
diabete <- rbind(
  P_DIQ[, c("SEQN","DIQ010","DIQ050","DIQ070")],
  DIQ_J[, c("SEQN","DIQ010","DIQ050","DIQ070")]
) %>% rename(
  diabetes = DIQ010,   # 1是, 2否, 7拒, 9不详
  Insulin  = DIQ050,   # 1是, 2否
  sugar    = DIQ070    # 1是, 2否
)
data1 <- full_join(glucose, diabete, by = "SEQN")

## ========= 3) 糖化 & 胰岛素 =========
GHB_L <- read.xport("C:/zq/Nhanes/糖化血红蛋白/GHB_L.xpt")
P_GHB <- read.xport("C:/zq/Nhanes/糖化血红蛋白/P_GHB.xpt")
GHB_J <- read.xport("C:/zq/Nhanes/糖化血红蛋白/GHB_J.xpt")
tanghua <- rbind(
  GHB_L[, c("SEQN","LBXGH")],
  P_GHB[, c("SEQN","LBXGH")],
  GHB_J[, c("SEQN","LBXGH")]
) %>% rename(HbA1c = LBXGH)   # %

INS_L <- read.xport("C:/zq/Nhanes/胰岛素/INS_L.xpt")
P_INS <- read.xport("C:/zq/Nhanes/胰岛素/P_INS.xpt")
INS_J <- read.xport("C:/zq/Nhanes/胰岛素/INS_J.xpt")
yidaosu <- rbind(
  INS_L[, c("SEQN","LBXIN")],
  P_INS[, c("SEQN","LBXIN")],
  INS_J[, c("SEQN","LBXIN")]
) %>% rename(FINS = LBXIN)    # μU/mL
data2 <- full_join(tanghua, yidaosu, by = "SEQN")

## ========= 4) 脂质 & CRP =========
TCHOL_L <- read.xport("C:/zq/Nhanes/胆固醇总计/TCHOL_L.xpt")
P_TCHOL <- read.xport("C:/zq/Nhanes/胆固醇总计/P_TCHOL.xpt")
TCHOL_J <- read.xport("C:/zq/Nhanes/胆固醇总计/TCHOL_J.xpt")
zongdgc <- rbind(
  TCHOL_L[, c("SEQN","LBDTCSI","LBXTC")],
  P_TCHOL[, c("SEQN","LBDTCSI","LBXTC")],
  TCHOL_J[, c("SEQN","LBDTCSI","LBXTC")]
) %>% rename(TC1 = LBDTCSI,  # mmol/L
             TC2 = LBXTC)    # mg/dL

P_TRIGLY <- read.xport("C:/zq/Nhanes/低密度脂蛋白/P_TRIGLY.xpt")
TRIGLY_J <- read.xport("C:/zq/Nhanes/低密度脂蛋白/TRIGLY_J.xpt")
dizhidanbai <- rbind(
  P_TRIGLY[, c("SEQN","LBDTRSI","LBXTR","LBDLDNSI","LBDLDLN")],
  TRIGLY_J[, c("SEQN","LBDTRSI","LBXTR","LBDLDNSI","LBDLDLN")]
) %>% rename(
  TG1    = LBDTRSI,   # mmol/L
  TG2    = LBXTR,     # mg/dL
  LDL_c1 = LBDLDNSI,  # mmol/L
  LDL_c2 = LBDLDLN    # mg/dL
)
data3 <- full_join(zongdgc, dizhidanbai, by = "SEQN")

HDL_L <- read.xport("C:/zq/Nhanes/高密度脂蛋白/HDL_L.xpt")
P_HDL <- read.xport("C:/zq/Nhanes/高密度脂蛋白/P_HDL.xpt")
HDL_J <- read.xport("C:/zq/Nhanes/高密度脂蛋白/HDL_J.xpt")
gaomdzhidanbai <- rbind(
  HDL_L[, c("SEQN","LBDHDDSI","LBDHDD")],
  P_HDL[, c("SEQN","LBDHDDSI","LBDHDD")],
  HDL_J[, c("SEQN","LBDHDDSI","LBDHDD")]
) %>% rename(HDL_c1 = LBDHDDSI,  # mmol/L
             HDL_c2 = LBDHDD)    # mg/dL

HSCRP_L <- read.xport("C:/zq/Nhanes/C反应蛋白/HSCRP_L.xpt")
P_HSCRP <- read.xport("C:/zq/Nhanes/C反应蛋白/P_HSCRP.xpt")
HSCRP_J <- read.xport("C:/zq/Nhanes/C反应蛋白/HSCRP_J.xpt")
cfan <- rbind(
  HSCRP_L[, c("SEQN","LBXHSCRP")],
  P_HSCRP[, c("SEQN","LBXHSCRP")],
  HSCRP_J[, c("SEQN","LBXHSCRP")]
) %>% rename(CRP = LBXHSCRP)
data4 <- full_join(gaomdzhidanbai, cfan, by = "SEQN")

## ========= 5) 金属 & 生化 =========
P_PBCD <- read.xport("C:/zq/Nhanes/铅、镉、总汞、硒和锰-血液/P_PBCD.xpt")
PBCD_L <- read.xport("C:/zq/Nhanes/铅、镉、总汞、硒和锰-血液/PBCD_L.xpt")
PBCD_J <- read.xport("C:/zq/Nhanes/铅、镉、总汞、硒和锰-血液/PBCD_J.xpt")
jinshu <- rbind(
  PBCD_L[, c("SEQN","LBXBCD")],
  P_PBCD[, c("SEQN","LBXBCD")],
  PBCD_J[, c("SEQN","LBXBCD")]
) %>% rename(Cd = LBXBCD)

P_BIOPRO <- read.xport("C:/zq/Nhanes/标准生化/P_BIOPRO.xpt")
BIOPRO_J <- read.xport("C:/zq/Nhanes/标准生化/BIOPRO_J.xpt")
shenhua <- rbind(
  P_BIOPRO[, c("SEQN","LBXSATSI","LBXSASSI","LBXSGTSI","LBXSAL","LBXSCR")],
  BIOPRO_J[, c("SEQN","LBXSATSI","LBXSASSI","LBXSGTSI","LBXSAL","LBXSCR")]
) %>% rename(
  ALT     = LBXSATSI,  # U/L
  AST     = LBXSASSI,  # U/L
  GGT     = LBXSGTSI,  # U/L
  Albumin = LBXSAL,    # g/dL
  Scr     = LBXSCR     # mg/dL
)
data5 <- full_join(jinshu, shenhua, by = "SEQN")

## ========= 6) 肝成像 & CBC =========
LUX_L <- read.xport("C:/zq/Nhanes/肝成像/LUX_L.xpt")
P_LUX <- read.xport("C:/zq/Nhanes/肝成像/P_LUX.xpt")
LUX_J <- read.xport("C:/zq/Nhanes/肝成像/LUX_J.xpt")
ganchx <- rbind(
  LUX_L[, c("SEQN","LUXSMED","LUXCAPM")],
  P_LUX[, c("SEQN","LUXSMED","LUXCAPM")],
  LUX_J[, c("SEQN","LUXSMED","LUXCAPM")]
) %>% rename(LSM = LUXSMED,   # kPa
             CAP = LUXCAPM)   # dB/m

CBC_L <- read.xport("C:/zq/Nhanes/全血/CBC_L.xpt")
P_CBC <- read.xport("C:/zq/Nhanes/全血/P_CBC.xpt")
CBC_J <- read.xport("C:/zq/Nhanes/全血/CBC_J.xpt")
quanxue <- rbind(
  CBC_L[, c("SEQN","LBXPLTSI","LBDNENO","LBDLYMNO")],
  P_CBC[, c("SEQN","LBXPLTSI","LBDNENO","LBDLYMNO")],
  CBC_J[, c("SEQN","LBXPLTSI","LBDNENO","LBDLYMNO")]
) %>% rename(
  platelet    = LBXPLTSI,  # 10^3/μL
  neutrophils = LBDNENO,
  lymphocytes = LBDLYMNO
) %>% mutate(SII = platelet * neutrophils / lymphocytes)
data6 <- full_join(ganchx, quanxue, by = "SEQN")

## ========= 7) 身体测量 & 是否服用降压药 =========
BMX_L <- read.xport("C:/zq/Nhanes/身体测量/BMX_L.xpt")
P_BMX <- read.xport("C:/zq/Nhanes/身体测量/P_BMX.xpt")
BMX_J <- read.xport("C:/zq/Nhanes/身体测量/BMX_J.xpt")
shenti <- rbind(
  BMX_L[, c("SEQN","BMXHT","BMXWT","BMXBMI","BMXWAIST","BMXHIP")],
  P_BMX[, c("SEQN","BMXHT","BMXWT","BMXBMI","BMXWAIST","BMXHIP")],
  BMX_J[, c("SEQN","BMXHT","BMXWT","BMXBMI","BMXWAIST","BMXHIP")]
) %>% rename(
  height    = BMXHT,
  weiht     = BMXWT,
  BMI       = BMXBMI,
  waistline = BMXWAIST,
  hip       = BMXHIP
)

BPQ_L <- read.xport("C:/zq/Nhanes/血压/BPQ_L.xpt")
P_BPQ <- read.xport("C:/zq/Nhanes/血压/P_BPQ.xpt")
BPQ_J <- read.xport("C:/zq/Nhanes/血压/BPQ_J.xpt")
.bp_med_candidates <- c("BPQ040A","BPQ040","BPQ050A","BPQ050","BPQ150")
extract_bp_med <- function(df) {
  stopifnot("SEQN" %in% names(df))
  vars <- intersect(names(df), .bp_med_candidates)
  out <- data.frame(SEQN = df$SEQN)
  if (length(vars) == 0) { out$HMedcine <- NA_integer_; return(out) }
  m_rec <- as.data.frame(lapply(df[, vars, drop = FALSE], function(x) {
    ifelse(x == 1, 1L, ifelse(x == 2, 0L, NA_integer_))
  }))
  out$HMedcine <- apply(m_rec, 1, function(r) {
    if (all(is.na(r))) NA_integer_ else as.integer(any(r == 1, na.rm = TRUE))
  })
  out
}
bp_med_L <- extract_bp_med(BPQ_L)
bp_med_P <- extract_bp_med(P_BPQ)
bp_med_J <- extract_bp_med(BPQ_J)
xueya1 <- bind_rows(bp_med_L, bp_med_P, bp_med_J) %>%
  group_by(SEQN) %>%
  summarise(HMedcine = if (all(is.na(HMedcine))) NA_integer_ else max(HMedcine, na.rm = TRUE),
            .groups = "drop")
data7 <- full_join(shenti, xueya1, by = "SEQN")

## ========= 8) 血压（被告知/测量） & 高血压判定 =========
xueyaq <- rbind(
  BPQ_L[, c("SEQN","BPQ020","BPQ030")],
  P_BPQ[, c("SEQN","BPQ020","BPQ030")],
  BPQ_J[, c("SEQN","BPQ020","BPQ030")]
) %>% rename(xy1 = BPQ020, xy2 = BPQ030)   # 1=Yes, 2=No

bpxo_path1 <- "C:/zq/Nhanes/血压/P_BPXO.xpt"
bpxo_path2 <- "C:/zq/Nhanes/血压/BPXO_P.xpt"
if (file.exists(bpxo_path1)) {
  P_BPXO <- read.xport(bpxo_path1)
} else if (file.exists(bpxo_path2)) {
  P_BPXO <- read.xport(bpxo_path2)
} else {
  stop("未找到 P_BPXO 文件，请检查路径：", bpxo_path1, " 或 ", bpxo_path2)
}
BPXO_L <- read.xport("C:/zq/Nhanes/血压/BPXO_L.xpt")
BPXO_J <- read.xport("C:/zq/Nhanes/血压/BPXO_J.xpt")

xueyac <- rbind(
  BPXO_L[, c("SEQN","BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3")],
  P_BPXO[, c("SEQN","BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3")],
  BPXO_J[, c("SEQN","BPXOSY1","BPXODI1","BPXOSY2","BPXODI2","BPXOSY3","BPXODI3")]
)

xueyaall <- merge(xueyaq, xueyac, by = "SEQN", all.x = TRUE)
xueyaall$HNT <- ifelse(
  (xueyaall$xy1 == 1 | xueyaall$xy2 == 1) |
    (xueyaall$BPXOSY1 >= 140 | xueyaall$BPXODI1 >= 90) |
    (xueyaall$BPXOSY2 >= 140 | xueyaall$BPXODI2 >= 90) |
    (xueyaall$BPXOSY3 >= 140 | xueyaall$BPXODI3 >= 90),
  1, 0
)
xueyaall <- xueyaall[, c("SEQN","HNT")]
data8 <- full_join(data7, xueyaall, by = "SEQN")

## ========= 9) 心脑血管问卷 =========
MCQ_L <- read.xport("C:/zq/Nhanes/其他健康问卷/MCQ_L.xpt")
P_MCQ <- read.xport("C:/zq/Nhanes/其他健康问卷/P_MCQ.xpt")
MCQ_J <- read.xport("C:/zq/Nhanes/其他健康问卷/MCQ_J.xpt")
disease <- rbind(
  MCQ_L[, c("MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F","MCQ160L","MCQ220","SEQN")],
  P_MCQ[, c("MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160L","MCQ160F","MCQ220","SEQN")],
  MCQ_J[, c("MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160L","MCQ160F","MCQ220","SEQN")]
) %>% rename(
  CHF    = MCQ160B,
  CHD    = MCQ160C,
  AP     = MCQ160D,
  heart  = MCQ160E,
  stroke = MCQ160F,
  liver  = MCQ160L,
  cancer = MCQ220
) %>% mutate(
  CHF    = ifelse(CHF==1,1, ifelse(CHF %in% c(2,7,9),0,NA)),
  CHD    = ifelse(CHD==1,1, ifelse(CHD %in% c(2,7,9),0,NA)),
  AP     = ifelse(AP==1,1, ifelse(AP %in% c(2,7,9),0,NA)),
  heart  = ifelse(heart==1,1, ifelse(heart %in% c(2,7,9),0,NA)),
  stroke = ifelse(stroke==1,1, ifelse(stroke %in% c(2,7,9),0,NA)),
  liver  = ifelse(liver==1,1, ifelse(liver %in% c(2,7,9),0,NA)),
  cancer = ifelse(cancer==1,1, ifelse(cancer %in% c(2,7,9),0,NA))
)
data8 <- full_join(data8, disease, by = "SEQN")

## ========= 10) 阿司匹林 & 运动 =========
P_RXQASA <- read.xport("C:/zq/Nhanes/阿司匹林/P_RXQASA.xpt")
RXQASA_J <- read.xport("C:/zq/Nhanes/阿司匹林/RXQASA_J.xpt")
as <- rbind(
  P_RXQASA[, c("SEQN","RXQ510")],
  RXQASA_J[, c("SEQN","RXQ510")]
) %>% mutate(ax = ifelse(RXQ510 == 1, 1L, 0L)) %>% select(SEQN, ax)

PAQ_L <- read.xport("C:/zq/Nhanes/运动/PAQ_L.xpt")
P_PAQ <- read.xport("C:/zq/Nhanes/运动/P_PAQ.xpt")
PAQ_J <- read.xport("C:/zq/Nhanes/运动/PAQ_J.xpt")
PAQ_L <- PAQ_L[, c("SEQN","PAD790Q","PAD790U","PAD800")]
PAQ_L <- PAQ_L %>%
  mutate(across(
    c("PAD790Q","PAD790U","PAD800"),
    ~ suppressWarnings(as.numeric(
      ifelse(is.na(.), 0,
             ifelse(.=="W",1,
                    ifelse(.=="D",7,
                           ifelse(.=="M",1/4,
                                  ifelse(.=="Y",1/52, .)))))))
  )) %>%
  mutate(
    Product = PAD790Q * PAD790U * PAD800,
    Level   = dplyr::case_when(
      is.na(Product) ~ 1,
      Product < 150  ~ 1,
      Product <= 300 ~ 2,
      Product > 300  ~ 3,
      TRUE ~ NA_real_
    )
  )

PAQ2 <- rbind(
  P_PAQ[, c("SEQN","PAQ605","PAQ620")],
  PAQ_J[, c("SEQN","PAQ605","PAQ620")]
) %>% mutate(
  Level = ifelse(PAQ620==1 & PAQ605==1, 4,
                 ifelse(PAQ620==1, 2, ifelse(PAQ605==1, 3, 1))),
  Level = ifelse(Level==4,3, ifelse(Level %in% c(2,3),2,1))
) %>% select(SEQN, Level)

sport <- rbind(PAQ_L[, c("SEQN","Level")], PAQ2[, c("SEQN","Level")])
data9 <- full_join(as, sport, by = "SEQN")

## ========= 11) 丙肝 & 饮酒 =========
HEPC_J <- read.xport("C:/zq/Nhanes/丙肝/HEPC_J.xpt")
P_HEPC <- read.xport("C:/zq/Nhanes/丙肝/P_HEPC.xpt")
binggan <- rbind(
  HEPC_J[, c("SEQN","LBDHCI")],
  P_HEPC[, c("SEQN","LBDHCI")]
) %>% rename(BING = LBDHCI)

ALQ_J <- read.xport("C:/zq/Nhanes/饮酒/ALQ_J.xpt")
P_ALQ <- read.xport("C:/zq/Nhanes/饮酒/P_ALQ.xpt")
ALQ <- rbind(
  ALQ_J[, c("SEQN","ALQ121","ALQ130")],
  P_ALQ[, c("SEQN","ALQ121","ALQ130")]
)
# ALQ121=过去一年饮酒天数；ALQ130=饮酒日标准杯数；每杯约 14g 纯酒精
ALQ <- ALQ %>%
  mutate(
    weekly_alc_g = suppressWarnings(as.numeric(ALQ130)) * 14 *
      (suppressWarnings(as.numeric(ALQ121))/52)
  ) %>%
  transmute(SEQN, weekly_alc = weekly_alc_g)
data10 <- full_join(binggan, ALQ, by = "SEQN")

## ========= 12) 吸烟 & 肾脏问卷 =========
P_SMQ <- read.xport("C:/zq/Nhanes/吸烟/P_SMQ.xpt")
SMQ_J <- read.xport("C:/zq/Nhanes/吸烟/SMQ_J.xpt")
SMQ <- rbind(
  P_SMQ[, c("SEQN","SMQ020","SMQ040")],
  SMQ_J[, c("SEQN","SMQ020","SMQ040")]
) %>% mutate(
  Smoke = ifelse(SMQ020==2, 0,
                 ifelse(SMQ040==3, 1,
                        ifelse(SMQ040 %in% c(1,2), 2, NA)))
) %>% select(SEQN, Smoke)

KIQ_U_J <- read.xport("C:/zq/Nhanes/肾脏/KIQ_U_J.xpt")
P_KIQ_U  <- read.xport("C:/zq/Nhanes/肾脏/P_KIQ_U.xpt")
kidney <- rbind(
  KIQ_U_J[, c("SEQN","KIQ022","KIQ025")],
  P_KIQ_U[, c("SEQN","KIQ022","KIQ025")]
)  # 1=Yes
data11 <- full_join(SMQ, kidney, by = "SEQN")

## ========= 13) 人口学 & 处方药 =========
DEMO_L <- read.xport("C:/zq/Nhanes/人口学特征/DEMO_L.xpt")
P_DEMO <- read.xport("C:/zq/Nhanes/人口学特征/P_DEMO.xpt")
DEMO_J <- read.xport("C:/zq/Nhanes/人口学特征/DEMO_J.xpt")
DEMO_L$Year <- "2021-2023"
P_DEMO$Year <- "2017-2020"
DEMO_J$Year <- "2017-2018"

P_DEMO$Marital_status <- ifelse(P_DEMO$DMDMARTZ==3,1,
                                ifelse(P_DEMO$DMDMARTZ==1,2,
                                       ifelse(P_DEMO$DMDMARTZ %in% c(2,77,99),3,NA)))
DEMO_J$Marital_status <- ifelse(DEMO_J$DMDMARTL==5,1,
                                ifelse(DEMO_J$DMDMARTL %in% c(1,6),2,
                                       ifelse(DEMO_J$DMDMARTL %in% c(2,3,4),3,NA)))

P_DEMO <- P_DEMO %>% mutate(Poverty_level = case_when(
  INDFMPIR <= 1 ~ 1,
  INDFMPIR > 1 & INDFMPIR <= 4 ~ 2,
  INDFMPIR > 4 ~ 3,
  TRUE ~ NA_real_
))
DEMO_J <- DEMO_J %>% mutate(Poverty_level = case_when(
  INDFMPIR <= 1 ~ 1,
  INDFMPIR > 1 & INDFMPIR <= 4 ~ 2,
  INDFMPIR > 4 ~ 3,
  TRUE ~ NA_real_
))

P_DEMO$Eduction <- ifelse(P_DEMO$DMDEDUC2 %in% c(1,2),1,
                          ifelse(P_DEMO$DMDEDUC2==3,2,
                                 ifelse(P_DEMO$DMDEDUC2 %in% c(4,5),3,NA)))
DEMO_J$Eduction <- ifelse(DEMO_J$DMDEDUC3 %in% c(0:11,55,66),1,
                          ifelse(DEMO_J$DMDEDUC3 %in% c(12,13,14),2,
                                 ifelse(DEMO_J$DMDEDUC3==15,3,
                                        ifelse(DEMO_J$DMDEDUC2 %in% c(1,2),1,
                                               ifelse(DEMO_J$DMDEDUC2==3,2,
                                                      ifelse(DEMO_J$DMDEDUC2 %in% c(4,5),3,NA))))))

DEMO <- rbind(
  P_DEMO[, c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1",
             "Marital_status","Eduction","Poverty_level","Year")],
  DEMO_J[, c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH1",
             "Marital_status","Eduction","Poverty_level","Year")]
) %>% rename(Gender = RIAGENDR,
             Age    = RIDAGEYR,
             Race   = RIDRETH1)

P_RXQ_RX <- read.xport("C:/zq/Nhanes/处方药物/P_RXQ_RX.xpt")
RXQ_RX_J <- read.xport("C:/zq/Nhanes/处方药物/RXQ_RX_J.xpt")
yaowu <- rbind(
  P_RXQ_RX[, c("SEQN","RXDUSE","RXDRSC1","RXDRSC2","RXDRSC3")],
  RXQ_RX_J[, c("SEQN","RXDUSE","RXDRSC1","RXDRSC2","RXDRSC3")]
) %>% mutate(
  across(starts_with("RXDRSC"),
         ~ ifelse(is.na(.), "", toupper(trimws(as.character(.))))),
  cvd_code = grepl("^I2[0-5]|^I6[0-9]|^Z95", RXDRSC1) |
    grepl("^I2[0-5]|^I6[0-9]|^Z95", RXDRSC2) |
    grepl("^I2[0-5]|^I6[0-9]|^Z95", RXDRSC3)
) %>%
  group_by(SEQN) %>%
  summarise(yw = as.integer(any(cvd_code, na.rm = TRUE)),
            .groups = "drop")

data12 <- full_join(DEMO, yaowu, by = "SEQN")

## ========= 14) 大合并 & CVD/DM 标签 =========
data_list <- list(data1, data2, data3, data4, data5, data6,
                  data8, data9, data10, data11, data12)
data_list_unique <- data_list %>%
  map(~ dplyr::distinct(., SEQN, .keep_all = TRUE))
Data <- Reduce(function(x, y) dplyr::full_join(x, y, by = "SEQN"),
               data_list_unique) %>%
  chk("合并后")

# 删除旧 CVD 列
Data <- Data %>% dplyr::select(-dplyr::any_of(c("CVD_strict","CVD_broad")))

# DM 定义
Data <- Data %>%
  mutate(
    DM = ifelse(
      diabetes == 1 |
        (!is.na(FPG1)  & FPG1  >= 7) |
        (!is.na(HbA1c) & HbA1c >= 6.5) |
        Insulin == 1 |
        sugar   == 1,
      1L, 0L)
  )

# strict：MCQ160C/D/E/F 任一 = 1
strict_any <- (to01(Data$AP) + to01(Data$CHD) +
                 to01(Data$heart) + to01(Data$stroke)) > 0
Data$CVD_strict <- ifelse(strict_any, 1L, 0L)

CHF01 <- if ("CHF" %in% names(Data)) to01(Data$CHF) else rep(0L, nrow(Data))
Data$yw <- as.integer(!is.na(Data$yw) & Data$yw == 1L)
Data$ax <- as.integer(!is.na(Data$ax) & Data$ax == 1L)

Data$CVD_broad <- as.integer(
  strict_any | (CHF01 == 1L) | (Data$yw == 1L) |
    (incl_aspirin_in_broad & Data$ax == 1L)
)

Data$CVD_strict <- as.integer(Data$CVD_strict %in% 1)
Data$CVD_broad  <- as.integer(Data$CVD_broad  %in% 1)
stopifnot(all(Data$CVD_broad %in% 0:1),
          all(Data$CVD_strict %in% 0:1))
bad <- sum(Data$CVD_strict == 1 & Data$CVD_broad == 0, na.rm = TRUE)
if (bad > 0) stop(sprintf("逻辑错误：有 %d 行 strict=1 但 broad=0", bad))

## ========= 15) 年龄 + 基本排除（肝病/丙肝/肾病/恶性肿瘤/重度饮酒） =========
Data <- Data %>%
  filter(DM == 1, Age >= 30, Age <= 75) %>%
  chk("保留 DM==1 且年龄 30–75")

Data <- Data %>%
  mutate(
    liver_pos  = dplyr::case_when(
      liver == 1 ~ TRUE,
      liver %in% c(0) ~ FALSE,
      TRUE ~ NA
    ),
    BING_pos   = dplyr::case_when(
      BING == 1 ~ TRUE,
      BING %in% c(0,2) ~ FALSE,
      TRUE ~ NA
    ),
    kidney_pos = dplyr::case_when(
      KIQ022 == 1 | KIQ025 == 1 ~ TRUE,
      KIQ022 %in% c(2, NA) &
        (is.na(KIQ025) | KIQ025 == 2) ~ FALSE,
      TRUE ~ NA
    ),
    cancer_pos = dplyr::case_when(
      cancer == 1 ~ TRUE,
      cancer %in% c(0) ~ FALSE,
      TRUE ~ NA
    )
  ) %>%
  filter(is.na(liver_pos)  | liver_pos  == FALSE) %>% chk("剔除肝病阳性") %>%
  filter(is.na(BING_pos)   | BING_pos   == FALSE) %>% chk("剔除丙肝阳性") %>%
  filter(is.na(kidney_pos) | kidney_pos == FALSE) %>% chk("剔除肾衰/透析阳性") %>%
  filter(is.na(cancer_pos) | cancer_pos == FALSE) %>% chk("剔除恶性肿瘤阳性")

# 重度饮酒（仅在有数时剔除；单位：克/周）
Data <- Data %>%
  filter(!(Gender == 1 & !is.na(weekly_alc) & weekly_alc > 210)) %>%  # 男 >210 g/周
  filter(!(Gender == 2 & !is.na(weekly_alc) & weekly_alc > 140)) %>%  # 女 >140 g/周
  chk("剔除重度饮酒")


# =========================================================
# 2) 变量规范化命名 & 各种指数（全部基于原始缺失的 Data）
#    —— 与建模脚本 Data_out 一致，但不插补
# =========================================================

calculate_eGFR_2021 <- function(Scr, Age, sex){
  kappa <- ifelse(sex == "Female", 0.7, 0.9)
  alpha <- ifelse(sex == "Female", -0.241, -0.302)
  142 * pmin(Scr/kappa, 1)^alpha *
    pmax(Scr/kappa, 1)^(-1.200) *
    (0.9938^Age) *
    ifelse(sex == "Female", 1.012, 1)
}

Data_std <- Data %>%
  mutate(
    ## 基础特征
    sex               = factor(Gender, levels = c(1,2),
                               labels = c("Male","Female")),
    age               = Age,
    race              = Race,
    height_cm         = height,
    weight_kg         = weiht,
    bmi               = BMI,
    waist_cm          = waistline,
    hip_cm            = hip,
    ## 生化（统一 mg/dL / 标准单位）
    fpg_mgdl          = FPG2,
    hba1c_pct         = HbA1c,
    insulin_uuml      = FINS,
    tc_mgdl           = TC2,
    hdl_mgdl          = HDL_c2,
    ldl_mgdl          = LDL_c2,
    tg_mgdl           = TG2,
    crp_mgL           = CRP,
    alt_uL            = ALT,
    ast_uL            = AST,
    ggt_uL            = GGT,
    albumin_gdl       = Albumin,
    scr_mgdl          = Scr,
    blood_cadmium_ugL = Cd,
    platelets_10e3_uL = platelet,
    neutrophils_10e3_uL = neutrophils,
    lymphocytes_10e3_uL = lymphocytes,
    sii               = SII,
    lsm_kpa           = LSM,
    cap_dbm           = CAP,
    ## 行为/血压/疾病标签
    pa_level          = Level,
    smoking_status    = factor(Smoke,
                               levels = c(0,1,2),
                               labels = c("Never","Former","Current")),
    antihypertensive_med   = HMedcine,
    hypertension_dx_or_meas = HNT,
    dm                = DM
  ) %>%
  ## 各类指数（注意：这里仍然是原始缺失；只要有一个原始缺失，该指数也会 NA）
  mutate(
    homa_ir        = (fpg_mgdl * 0.0555) * insulin_uuml / 22.5,
    tyg_index      = log((tg_mgdl * fpg_mgdl) / 2),
    mets_ir        = log(2*fpg_mgdl + tg_mgdl) * bmi / log(hdl_mgdl),
    tg_hdl_ratio   = tg_mgdl / hdl_mgdl,
    aip            = log10((tg_mgdl*0.01129) / (hdl_mgdl*0.02586)),
    fib4           = age * ast_uL / (platelets_10e3_uL * sqrt(alt_uL)),
    aar            = ast_uL / alt_uL,
    hsi            = 8*(alt_uL/ast_uL) + bmi +
      2*(sex=="Female") + 2*(dm==1),
    nfs            = -1.675 + 0.037*age + 0.094*bmi + 1.13*1 +
      0.99*(ast_uL/alt_uL) -
      0.013*platelets_10e3_uL - 0.66*albumin_gdl,
    apri           = ifelse(
      sex=="Male",
      (ast_uL/37)/platelets_10e3_uL*100,
      ifelse(
        sex=="Female",
        (ast_uL/31)/platelets_10e3_uL*100,
        NA
      )
    )
  )

Data_std$egfr <- calculate_eGFR_2021(
  Data_std$scr_mgdl,
  Data_std$age,
  as.character(Data_std$sex)
)
Data_std$egdr <- 21.158 -
  0.09*Data_std$waist_cm -
  3.407*Data_std$hypertension_dx_or_meas -
  0.551*Data_std$hba1c_pct


# =========================================================
# 3) 精简为“建模同款变量集” Data_out_raw（未插补）
# =========================================================

Data_out_raw <- Data_std %>%
  select(
    ## 标签
    dm, cvd_strict = CVD_strict, cvd_broad = CVD_broad,
    ## 人口学
    sex, age, race,
    ## 体测
    height_cm, weight_kg, bmi, waist_cm, hip_cm,
    ## 行为
    smoking_status, pa_level,
    ## 血压与用药
    antihypertensive_med, hypertension_dx_or_meas,
    ## 生化
    fpg_mgdl, hba1c_pct, insulin_uuml,
    tc_mgdl, ldl_mgdl, hdl_mgdl, tg_mgdl,
    crp_mgL, alt_uL, ast_uL, ggt_uL, albumin_gdl, scr_mgdl,
    blood_cadmium_ugL,
    platelets_10e3_uL, neutrophils_10e3_uL, lymphocytes_10e3_uL, sii,
    ## 肝超声
    lsm_kpa, cap_dbm,
    ## 派生指数
    homa_ir, tyg_index, mets_ir, tg_hdl_ratio, aip,
    fib4, aar, hsi, nfs, apri, egfr, egdr
  ) %>%
  chk("最终可用数据（未插补版本，用于 S1）")


# =========================================================
# 4) label_map：与主脚本完全一致的 SCI 列名映射
# =========================================================

label_map <- c(
  dm  = "Diabetes (DM), 0/1",
  cvd_strict = "CVD (strict), 0/1",
  cvd_broad  = "CVD (broad), 0/1",
  sex = "Sex",
  age = "Age, years",
  race = "Race/ethnicity",
  height_cm = "Height, cm",
  weight_kg = "Weight, kg",
  bmi = "BMI, kg/m^2",
  waist_cm = "Waist circumference, cm",
  hip_cm = "Hip circumference, cm",
  smoking_status = "Smoking status",
  pa_level = "Physical activity level",
  antihypertensive_med = "Antihypertensive medication, 0/1",
  hypertension_dx_or_meas = "Hypertension (diagnosed or measured), 0/1",
  fpg_mgdl = "Fasting plasma glucose, mg/dL",
  hba1c_pct = "HbA1c, %",
  insulin_uuml = "Fasting insulin, μU/mL",
  tc_mgdl = "Total cholesterol, mg/dL",
  ldl_mgdl = "LDL-C, mg/dL",
  hdl_mgdl = "HDL-C, mg/dL",
  tg_mgdl = "Triglycerides, mg/dL",
  crp_mgL = "C-reactive protein, mg/L",
  alt_uL = "ALT, U/L",
  ast_uL = "AST, U/L",
  ggt_uL = "GGT, U/L",
  albumin_gdl = "Albumin, g/dL",
  scr_mgdl = "Serum creatinine, mg/dL",
  blood_cadmium_ugL = "Blood cadmium, μg/L",
  platelets_10e3_uL = "Platelets, ×10^3/μL",
  neutrophils_10e3_uL = "Neutrophils, ×10^3/μL",
  lymphocytes_10e3_uL = "Lymphocytes, ×10^3/μL",
  sii = "Systemic immune-inflammation index (SII)",
  lsm_kpa = "Liver stiffness (LSM), kPa",
  cap_dbm = "Controlled attenuation parameter (CAP), dB/m",
  homa_ir = "HOMA-IR",
  tyg_index = "TyG index",
  mets_ir = "METS-IR",
  tg_hdl_ratio = "Triglyceride-to-HDL ratio",
  aip = "Atherogenic index of plasma (AIP)",
  fib4 = "FIB-4",
  aar = "AST/ALT ratio (AAR)",
  hsi = "Hepatic steatosis index (HSI)",
  nfs = "NAFLD fibrosis score (NFS)",
  apri = "APRI",
  egfr = "eGFR, mL/min/1.73 m^2",
  egdr = "eGDR, mg/kg/min"
)

missing_labels <- setdiff(names(Data_out_raw), names(label_map))
if (length(missing_labels) > 0) {
  stop("以下列缺少 SCI 列名映射，请在 label_map 中补充：\n  ",
       paste(missing_labels, collapse = ", "))
}

# =========================================================
# 5) 基于 Data_out_raw 计算 S1：分布 + 原始缺失率
# =========================================================

# 自动判断变量类型
infer_type <- function(x) {
  x_no_na <- x[!is.na(x)]
  n_uniq  <- length(unique(x_no_na))
  if (is.numeric(x) || is.integer(x)) {
    if (n_uniq <= 2) {
      return("Binary")
    } else if (n_uniq <= 10) {
      return("Categorical")
    } else {
      return("Continuous")
    }
  } else {
    if (n_uniq <= 2) "Binary" else "Categorical"
  }
}

# 从 SCI 名称中抽取单位（逗号后）
extract_unit <- function(label) {
  if (str_detect(label, ",")) {
    unit <- str_trim(sub(".*?,", "", label))
    return(unit)
  } else {
    return("")
  }
}

# 编码说明（可按需要继续补充）
coding_dict <- list(
  dm = "1 = included diabetes cases; 0 = excluded",
  cvd_strict = "0 = no CHD/stroke; 1 = self-reported CHD, angina, MI, heart failure, or stroke",
  cvd_broad  = "0 = no major CVD; 1 = CVD_strict or CHF or CVD-related prescription or (optional) aspirin use",
  sex        = "Male / Female",
  smoking_status = "Never / Former / Current smoker",
  antihypertensive_med = "0 = not on antihypertensive medication; 1 = on antihypertensive medication",
  hypertension_dx_or_meas = "0 = no hypertension; 1 = self-reported diagnosis or SBP ≥140 mmHg or DBP ≥90 mmHg",
  pa_level   = "1 = inactive, 2 = moderate, 3 = vigorous (derived from PAQ)"
  # 你可以按需要继续补充 race、egfr 等
)

get_coding <- function(var_name) {
  if (var_name %in% names(coding_dict)) coding_dict[[var_name]] else ""
}

# NHANES 分布 + 范围
summarize_nhanes_var <- function(x, type) {
  x_no_na <- x[!is.na(x)]
  n <- length(x_no_na)
  if (n == 0) return("All values missing")
  
  if (type == "Continuous") {
    m   <- mean(x_no_na)
    sdv <- stats::sd(x_no_na)
    med <- stats::median(x_no_na)
    q1  <- stats::quantile(x_no_na, 0.25)
    q3  <- stats::quantile(x_no_na, 0.75)
    mn  <- min(x_no_na)
    mx  <- max(x_no_na)
    sprintf(
      "mean %.2f ± %.2f; median %.2f [%.2f–%.2f]; range %.2f–%.2f (n = %d)",
      m, sdv, med, q1, q3, mn, mx, n
    )
  } else {
    tab <- table(x_no_na)
    pct <- round(100 * tab / sum(tab), 1)
    parts <- paste0(names(tab), ": ", as.integer(tab),
                    " (", pct, "%)")
    paste0(paste(parts, collapse = "; "),
           sprintf(" (n = %d)", n))
  }
}

# ---------- 外部验证集（可选） ----------

external <- NULL
if (nzchar(path_external) && file.exists(path_external)) {
  external <- readr::read_csv(path_external, show_col_types = FALSE)
  message("外部验证集已读取：", path_external)
} else {
  message("未指定或未找到外部验证集，S1 中外部部分留空。")
}

vars_for_S1 <- names(Data_out_raw)

# 默认假设外部队列列名与 NHANES 一样；
# 若不同，可在此处修改，例如 external_name_map <- c(age = "Age", hba1c_pct = "HbA1c", ...)
external_name_map <- setNames(vars_for_S1, vars_for_S1)

summarize_external_var <- function(var_name_nhanes, type) {
  if (is.null(external)) return("")
  ext_name <- external_name_map[[var_name_nhanes]]
  if (is.null(ext_name) || !ext_name %in% names(external)) return("")
  x <- external[[ext_name]]
  x_no_na <- x[!is.na(x)]
  n <- length(x_no_na)
  if (n == 0) return("All values missing")
  
  if (type == "Continuous") {
    m   <- mean(x_no_na)
    sdv <- stats::sd(x_no_na)
    med <- stats::median(x_no_na)
    q1  <- stats::quantile(x_no_na, 0.25)
    q3  <- stats::quantile(x_no_na, 0.75)
    mn  <- min(x_no_na)
    mx  <- max(x_no_na)
    sprintf(
      "mean %.2f ± %.2f; median %.2f [%.2f–%.2f]; range %.2f–%.2f (n = %d)",
      m, sdv, med, q1, q3, mn, mx, n
    )
  } else {
    tab <- table(x_no_na)
    pct <- round(100 * tab / sum(tab), 1)
    parts <- paste0(names(tab), ": ", as.integer(tab),
                    " (", pct, "%)")
    paste0(paste(parts, collapse = "; "),
           sprintf(" (n = %d)", n))
  }
}

miss_rate_external <- function(var_name_nhanes) {
  if (is.null(external)) return(NA_real_)
  ext_name <- external_name_map[[var_name_nhanes]]
  if (is.null(ext_name) || !ext_name %in% names(external)) return(NA_real_)
  round(mean(is.na(external[[ext_name]])) * 100, 1)
}

# ---------- 逐变量生成 S1 行 ----------

# 若你已确定紧凑 K=6 模型的变量，可以在此填写（机器名），自动标记 (final model)
final_model_vars <- c(
  # 示例: "age", "hba1c_pct", "waist_cm", "homa_ir", "fib4", "lsm_kpa"
)

s1_table <- purrr::map_dfr(vars_for_S1, function(v) {
  x <- Data_out_raw[[v]]
  type <- infer_type(x)
  sci_label <- label_map[[v]]
  unit <- extract_unit(sci_label)
  coding <- get_coding(v)
  nhanes_dist <- summarize_nhanes_var(x, type)
  nhanes_miss <- miss_rate(x)
  
  ext_name  <- if (!is.null(external)) external_name_map[[v]] else ""
  ext_dist  <- summarize_external_var(v, type)
  ext_miss  <- miss_rate_external(v)
  
  is_final <- v %in% final_model_vars
  predictor_name <- if (is_final) {
    paste0(sci_label, " (final model)")
  } else {
    sci_label
  }
  
  tibble::tibble(
    Predictor                      = predictor_name,
    NHANES_variable_name          = v,
    External_variable_name        = ifelse(is.null(external), "", ext_name),
    Unit                          = unit,
    Type                          = type,
    Coding_or_derivation          = coding,
    NHANES_distribution_and_range = nhanes_dist,
    External_distribution_and_range = ext_dist,
    Missingness_NHANES_percent    = nhanes_miss,
    Missingness_External_percent  = ext_miss
  )
})

# ---------- 导出 S1 CSV（带 BOM，方便 Excel） ----------

write_csv_with_bom <- function(df, path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_excel_csv(df, path, na = "")
  } else {
    tmp <- tempfile(fileext = ".csv")
    utils::write.csv(df, tmp, row.names = FALSE,
                     fileEncoding = "UTF-8", na = "")
    con_in  <- file(tmp, open = "rb")
    raw_dat <- readBin(con_in, what = "raw",
                       n = file.info(tmp)$size)
    close(con_in)
    con_out <- file(path, open = "wb")
    writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con_out) # UTF-8 BOM
    writeBin(raw_dat, con_out)
    close(con_out)
  }
}

out_path <- "S1_variable_mapping.csv"
write_csv_with_bom(s1_table, out_path)

cat("\n=====================================\n")
cat("S1 变量映射表已生成：", out_path, "\n")
cat("变量数（行数）：", nrow(s1_table), "\n")
cat("Missingness_NHANES_percent 反映的是【未插补前】的缺失率。\n")
cat("请在 Excel 中打开后，手动补充 Coding 列和个别说明即可投稿使用。\n")
cat("=====================================\n")

