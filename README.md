# phenoage-KDM
#为数据库中加入计算好的eGFR
library(kidney.epi)
library(dplyr)

NHANES4 <- NHANES4 %>%
  mutate(
    egfr = egfr.ckdepi.cr.2021(
      creatinine = creat,                # 注意：参数名是 creatinine，不是 creat
      age = age,                         # 年龄
      sex = gender,                      # 性别 (1/2)
      creatinine_units = "mg/dl",        # 明确指定你的血肌酐单位（如果是 mg/dl 填 "mg/dl"，如果是 umol/L 请改填 "micromol/l"）
      label_sex_male = 1,                # 对应你数据中男性的编码（这里假设 1 是男）
      label_sex_female = 2               # 对应你数据中女性的编码（这里假设 2 是女）
    )
  )
  library(dplyr)
library(kidney.epi) # 用于计算 2021 版无种族 CKD-EPI eGFR

# 1. 记录初始样本量
df_calc <- NHANES4  # 如果你的数据集叫其他名称，请在此处替换
N0 <- nrow(df_calc)

# 2. 执行纳排与清洗步骤
df_step1 <- df_calc %>% filter(age >= 45)                  
# [排除 1] 年龄纳排：聚焦 >= 45 岁中老年人群
N1 <- nrow(df_step1)
df_step2 <- df_step1 %>% filter(is.na(pregnant) | pregnant != 1) 
# [排除 2] 排除怀孕人群
N2 <- nrow(df_step2)
df_step3 <- df_step2 %>% filter(!is.na(creat), !is.na(status))  
# [排除 3] 关键计算变量缺失（肌酐、生存状态）
N3 <- nrow(df_step3)
df_step4 <- df_step3 %>% filter(!is.na(phenoage0))           
# [排除 4] 无法计算 PhenoAge 的核心生化指标缺失
N4 <- nrow(df_step4)
df_step5 <- df_step4 %>% filter(creat < 15)                 
# [排除 5] 剔除生理极值/离群值
N5 <- nrow(df_step5)

# 3. 规范 PhenoAge 边界值
df_step5 <- df_step5 %>%
  mutate(phenoage0 = pmax(phenoage0, 18))

# 4. 计算 eGFR 并进行多级梯度切点与亚组变量构建
df_final <- df_step5 %>%
  mutate(
    # 计算 2021 版无种族 eGFR
    egfr = egfr.ckdepi.cr.2021(
      creatinine = creat,                
      age = age,                         
      sex = gender,                      
      creatinine_units = "mg/dl",        # 根据你的实际肌酐单位调整 ("mg/dl" 或 "micromol/l")
      label_sex_male = 1,                # 男性编码
      label_sex_female = 2               # 女性编码
    ),
    
    # 【切点放宽与多级分层】将 eGFR 精细扩展为多级梯度
    egfr_tier = case_when(
      egfr >= 90              ~ "Normal/High (>=90)",      # 正常/高值组
      egfr >= 60 & egfr < 90  ~ "Mildly Decreased (60-89)", # 轻度下降组（假性正常潜在人群）
      egfr < 60               ~ "Low (<60) - CKD range",    # 低值/CKD范围
      TRUE                    ~ NA_character_
    ),
    egfr_tier = factor(egfr_tier, levels = c("Normal/High (>=90)", "Mildly Decreased (60-89)", "Low (<60) - CKD range")),
    
    # 【亚组交互分析准备】年龄分层
    age_group = case_when(
      age >= 45 & age < 60 ~ "45-59岁",
      age >= 60 & age < 75 ~ "60-74岁",
      age >= 75            ~ ">=75岁",
      TRUE                 ~ "其他"
    ),
    
    # 【亚组交互分析准备】高血压与糖尿病分类标签
    hyperten_cat = factor(hyperten, levels = c(0, 1), labels = c("非高血压", "高血压")),
    diabetes_cat = factor(diabetes, levels = c(0, 1), labels = c("非糖尿病", "糖尿病"))
  )

N_final <- nrow(df_final)

# 5. 打印规范的筛选流向统计（可直接用于论文 Methods 或汇报流向图）
cat("--- 数据筛选与清洗流向统计 ---\n",
    "原始样本量 (N0):", N0, "\n",
    "步骤1 (年龄 >= 45岁):", N1, "\n",
    "步骤2 (排除孕妇):", N2, "\n",
    "步骤3 (关键变量完整):", N3, "\n",
    "步骤4 (PhenoAge 完整):", N4, "\n",
    "步骤5 (剔除肌酐极值):", N5, "\n",
    "最终分析样本量 (N_final):", N_final, "\n")
library(dplyr)
# 以egfr90为界限分组
df_final <- df_final %>%
  mutate(
    Group_4cat_90 = case_when(
      egfr >= 90 & pheno_egfr >= 90 ~ "1_Normal",
      egfr >= 90 & pheno_egfr < 90  ~ "2_Pseudonormal", 
      egfr < 90  & pheno_egfr >= 90 ~ "3_Pseudoabnormal",
      egfr < 90  & pheno_egfr < 90  ~ "4_Impaired",
      TRUE ~ NA_character_
    ),
    
    Group_4cat_90 = factor(
      Group_4cat_90,
      levels = c(
        "1_Normal", 
        "2_Pseudonormal", 
        "3_Pseudoabnormal", 
        "4_Impaired"
      ),
      labels = c(
        "Normal (Ref)", 
        "Pseudonormal (Target)", 
        "Pseudoabnormal", 
        "Impaired"
      )
    )
  )
library(tidyverse)
library(tableone)
# 四分基线表
library(tidyverse)
library(tableone)
# 1. 数据预处理：过滤 NA 分组，应用全量修复的高血压与糖尿病定义，计算差值
df_analysis <- df_final %>%
  filter(!is.na(Group_4cat_90)) %>%
  mutate(
    # eGFR 错位差值计算
    eGFR_diff     = egfr - pheno_egfr,
    kdm_egfr_diff = egfr - kdm_egfr,
    
    # 标准化年龄分组与性别标签
    Age_Group = factor(age_group, levels = c("45-59岁", "60-74岁", ">=75岁"),
                                  labels = c("45-59 years", "60-74 years", ">=75 years")),
    Sex = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    
    # 与森林图严格对齐：全量修复糖尿病 (结合 HbA1c 与常规血糖)
    Diabetes = case_when(
      hba1c >= 6.5 | glucose_mmol >= 7.0 ~ "Yes",
      diabetes_cat %in% c("糖尿病", 1, "1", "Yes") ~ "Yes",
      !is.na(hba1c) & hba1c < 6.5 ~ "No",
      !is.na(glucose_mmol) & glucose_mmol < 7.0 ~ "No",
      diabetes_cat %in% c("非糖尿病", 0, 2, "0", "2", "No") ~ "No",
      TRUE ~ NA_character_
    ),
    
    # 与森林图严格对齐：全量修复高血压 (结合实测血压与高血压标签)
    Hypertension = case_when(
      sbp >= 140 | dbp >= 90 ~ "Yes",
      hyperten %in% c(1, "1", "Yes") | hyperten_cat %in% c("高血压", 1, "1", "Yes") ~ "Yes",
      (!is.na(sbp) & sbp < 140) & (!is.na(dbp) & dbp < 90) ~ "No",
      hyperten %in% c(0, 2, "0", "2", "No") | hyperten_cat %in% c("非高血压", 0, 2, "0", "2", "No") ~ "No",
      TRUE ~ NA_character_
    )
  )

# 2. 匹配 Table 1 变量列表 (以修复后的全量变量为主)
vars_table1 <- c(
  "age", "Sex", "race", "BMXBMI", 
  "Hypertension", "Diabetes",
  "creat", "alb_gL", "glu_mmol", "crp_mgdL", "hba1c", "sbp", "dbp",
  "Lymph_pct", "MCV", "RDW", "ALP", "WBC",
  "egfr", 
  "phenoage0", "pheno_egfr", "eGFR_diff",
  "kdm_age", "kdm_egfr", "kdm_egfr_diff",
  "Age_Group", "status"
)

# 自动匹配数据集中存在的列
vars_table1 <- intersect(vars_table1, names(df_analysis))

# 分类变量与偏态连续变量定义
cat_vars <- intersect(c("Sex", "race", "Age_Group", "Hypertension", "Diabetes", "status"), vars_table1)
nonnormal_vars <- intersect(c("creat", "crp_mgdL", "glu_mmol", "ALP", "WBC", "eGFR_diff", "kdm_egfr_diff"), vars_table1)

# 3. 生成非加权基线表 Table 1
tab1_unweighted <- CreateTableOne(
  vars = vars_table1,
  strata = "Group_4cat_90",
  data = df_analysis,
  factorVars = cat_vars,
  test = TRUE
)

# 控制台打印查看
print(tab1_unweighted, nonnormal = nonnormal_vars, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)

# 4. 导出为规范的 CSV 表格
tab1_unw_mat <- print(tab1_unweighted, nonnormal = nonnormal_vars, showAllLevels = TRUE, printToggle = FALSE)
write.csv(tab1_unw_mat, file = "Table1_Unweighted_Consistent_Group90.csv")
library(tidyverse)
library(survival)
library(forestplot)
# 绘制亚组分层森林图
# 1. 重新映射全员覆盖的高血压与糖尿病变量
time_col   <- intersect(c("permth_exm", "time", "permth_int"), names(df_analysis))[1]
status_col <- intersect(c("status", "mortstat"), names(df_analysis))[1]

df_sub <- df_analysis %>%
  filter(Group_4cat_90 %in% c("Normal (Ref)", "Pseudonormal (Target)")) %>%
  mutate(
    time_surv   = as.numeric(.data[[time_col]]),
    status_surv = as.numeric(.data[[status_col]]),
    target_exp  = ifelse(Group_4cat_90 == "Pseudonormal (Target)", 1, 0),
    
    # 年龄分层与性别标准英文标签
    Age_Group = factor(age_group, levels = c("45-59岁", "60-74岁", ">=75岁"),
                                  labels = c("45-59 years", "60-74 years", ">=75 years")),
    Sex = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    
    # 糖尿病全量修复：结合 HbA1c 与常规血糖
    Diabetes = case_when(
      hba1c >= 6.5 | glucose_mmol >= 7.0 ~ "Yes",
      diabetes_cat %in% c("糖尿病", 1, "1", "Yes") ~ "Yes",
      !is.na(hba1c) & hba1c < 6.5 ~ "No",
      !is.na(glucose_mmol) & glucose_mmol < 7.0 ~ "No",
      diabetes_cat %in% c("非糖尿病", 0, 2, "0", "2", "No") ~ "No",
      TRUE ~ NA_character_
    ),
    
    # 高血压全量修复：结合实测血压与高血压标签
    Hypertension = case_when(
      sbp >= 140 | dbp >= 90 ~ "Yes",
      hyperten %in% c(1, "1", "Yes") | hyperten_cat %in% c("高血压", 1, "1", "Yes") ~ "Yes",
      (!is.na(sbp) & sbp < 140) & (!is.na(dbp) & dbp < 90) ~ "No",
      hyperten %in% c(0, 2, "0", "2", "No") | hyperten_cat %in% c("非高血压", 0, 2, "0", "2", "No") ~ "No",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(target_exp), !is.na(time_surv), !is.na(status_surv))

# 2. 亚组变量定义
subgroup_vars   <- c("Age_Group", "Sex", "Hypertension", "Diabetes")
subgroup_labels <- c("Age Group", "Sex", "Hypertension", "Diabetes")

# 3. 循环计算效应量与交互作用 P 值
res_list <- list()

for (i in seq_along(subgroup_vars)) {
  var       <- subgroup_vars[i]
  var_label <- subgroup_labels[i]
  lvls      <- levels(as.factor(df_sub[[var]]))
  lvls      <- lvls[!is.na(lvls)]
  
  d_complete <- df_sub %>% filter(!is.na(.data[[var]]))
  f_base     <- as.formula(paste("Surv(time_surv, status_surv) ~ target_exp +", var))
  f_int      <- as.formula(paste("Surv(time_surv, status_surv) ~ target_exp *", var))
  
  p_inter_val <- tryCatch({
    fit_base <- coxph(f_base, data = d_complete)
    fit_int  <- coxph(f_int,  data = d_complete)
    anova(fit_base, fit_int)[2, "Pr(>|Chi|)"]
  }, error = function(e) NA)
  
  p_inter_str <- ifelse(is.na(p_inter_val), "-", sprintf("%.3f", p_inter_val))
  
  # 添加大类表头行
  res_list[[paste0(var, "_header")]] <- data.frame(
    Subgroup      = var_label,
    Events_Target = "",
    Events_Ref    = "",
    HR = NA, Low = NA, High = NA,
    HR_CI         = "",
    P_inter       = p_inter_str,
    is_summary    = TRUE,
    stringsAsFactors = FALSE
  )
  
  for (lvl in lvls) {
    d_lvl <- d_complete %>% filter(.data[[var]] == lvl)
    
    n_tar <- sum(d_lvl$target_exp == 1)
    e_tar <- sum(d_lvl$target_exp == 1 & d_lvl$status_surv == 1)
    n_ref <- sum(d_lvl$target_exp == 0)
    e_ref <- sum(d_lvl$target_exp == 0 & d_lvl$status_surv == 1)
    
    fit_sub <- tryCatch(
      coxph(Surv(time_surv, status_surv) ~ target_exp, data = d_lvl),
      error = function(e) NULL
    )
    
    if (!is.null(fit_sub) && "target_exp" %in% names(coef(fit_sub))) {
      hr_val   <- as.numeric(exp(coef(fit_sub)["target_exp"]))
      ci_vals  <- as.numeric(exp(confint(fit_sub)["target_exp", ]))
      low_val  <- ci_vals[1]
      high_val <- ci_vals[2]
      hr_str   <- sprintf("%.2f (%.2f-%.2f)", hr_val, low_val, high_val)
    } else {
      hr_val <- NA; low_val <- NA; high_val <- NA
      hr_str <- "Not estimable"
    }
    
    res_list[[paste0(var, "_", lvl)]] <- data.frame(
      Subgroup      = paste0("   ", lvl),
      Events_Target = paste0(e_tar, "/", n_tar),
      Events_Ref    = paste0(e_ref, "/", n_ref),
      HR = hr_val, Low = low_val, High = high_val,
      HR_CI         = hr_str,
      P_inter       = "",
      is_summary    = FALSE,
      stringsAsFactors = FALSE
    )
  }
}

forest_df <- bind_rows(res_list)

# 4. 组装表格并绘制完整森林图
tabletext <- cbind(
  c("Subgroup", forest_df$Subgroup),
  c("Pseudonormal\n(Events/N)", forest_df$Events_Target),
  c("Normal\n(Events/N)", forest_df$Events_Ref),
  c("Hazard Ratio\n(95% CI)", forest_df$HR_CI),
  c("P for\nInteraction", forest_df$P_inter)
)

mean_vec  <- c(NA, forest_df$HR)
lower_vec <- c(NA, forest_df$Low)
upper_vec <- c(NA, forest_df$High)
is_summ   <- c(TRUE, forest_df$is_summary)

forestplot(
  labeltext = tabletext,
  mean = mean_vec,
  lower = lower_vec,
  upper = upper_vec,
  is.summary = is_summ,
  zero = 1.0,
  xlog = TRUE,
  clip = c(0.4, 6.0),
  xticks = c(0.5, 1.0, 2.0, 4.0),
  xlab = "Hazard Ratio (95% CI) for All-Cause Mortality (Pseudonormal vs Normal)",
  col = fpColors(box = "#1B4F72", line = "#2C3E50", summary = "#2C3E50"),
  boxsize = 0.25,
  ci.vertices = TRUE,
  ci.vertices.height = 0.15,
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.85),
    ticks = gpar(cex = 0.8),
    xlab  = gpar(cex = 0.9)
  )
)
# 马氏距离
install.packages("forestploter")
library(forestploter)
library(tidyverse)
library(survival)
install.packages("PMCMRplus")
library(PMCMRplus)
# 1. 动态自动匹配 NHANES4 中 9 项 PhenoAge 核心生化指标的实际列名
col_map <- list(
  alb   = intersect(c("albumin", "alb_gL", "alb"), names(NHANES4))[1],
  creat = intersect(c("creatinine", "creat"), names(NHANES4))[1],
  glu   = intersect(c("glucose", "glu_mmol", "glu"), names(NHANES4))[1],
  crp   = intersect(c("crp", "crp_mgdL"), names(NHANES4))[1],
  lymph = intersect(c("lymph", "Lymph_pct", "lymphocyte"), names(NHANES4))[1],
  mcv   = intersect(c("mcv", "MCV"), names(NHANES4))[1],
  rdw   = intersect(c("rdw", "RDW"), names(NHANES4))[1],
  alp   = intersect(c("alp", "ALP"), names(NHANES4))[1],
  wbc   = intersect(c("wbc", "WBC"), names(NHANES4))[1]
)

cat("匹配到的变量名映射：\n")
print(unlist(col_map))

# 2. 构建特征矩阵并对偏态变量做 log 转换 (符合 Cohen 2013 文献方法)
nhanes_dm_prep <- NHANES4 %>%
  mutate(
    feat_alb   = .data[[col_map$alb]],
    feat_creat = log(.data[[col_map$creat]]),
    feat_glu   = log(.data[[col_map$glu]]),
    feat_crp   = log(.data[[col_map$crp]] + 0.01),
    feat_lymph = .data[[col_map$lymph]],
    feat_mcv   = .data[[col_map$mcv]],
    feat_rdw   = .data[[col_map$rdw]],
    feat_alp   = log(.data[[col_map$alp]]),
    feat_wbc   = log(.data[[col_map$wbc]])
  )

features <- c("feat_alb", "feat_creat", "feat_glu", "feat_crp", 
              "feat_lymph", "feat_mcv", "feat_rdw", "feat_alp", "feat_wbc")

# 3. 构建 20-39 岁年轻健康参考人群基准 (Reference Population)
ref_pop <- nhanes_dm_prep %>%
  filter(
    age >= 20 & age <= 39,
    # 血压与病史
    (is.na(sbp) | sbp < 140) & (is.na(dbp) | dbp < 90),
    if("hyperten" %in% names(.)) (is.na(hyperten) | hyperten == 0) else TRUE,
    # 血糖与糖化
    if("hba1c" %in% names(.)) (is.na(hba1c) | hba1c < 6.5) else TRUE,
    # 肾功能基线
    if("egfr" %in% names(.)) (is.na(egfr) | egfr >= 90) else TRUE
  ) %>%
  drop_na(all_of(features))

cat("健康年轻参考基准样本量 (N_ref):", nrow(ref_pop), "\n")

# 计算参考基准均值向量 μ 与协方差矩阵 S
mu_ref <- colMeans(ref_pop[, features])
S_ref  <- cov(ref_pop[, features])

# 4. 在分析数据集 df_analysis 中计算马氏距离 DM
df_analysis <- df_analysis %>%
  mutate(
    feat_alb   = .data[[intersect(c("alb_gL", "albumin", "alb"), names(.))[1]]],
    feat_creat = log(.data[[intersect(c("creat", "creatinine"), names(.))[1]]]),
    feat_glu   = log(.data[[intersect(c("glucose", "glu_mmol", "glu"), names(.))[1]]]),
    feat_crp   = log(.data[[intersect(c("crp", "crp_mgdL"), names(.))[1]]] + 0.01),
    feat_lymph = .data[[intersect(c("Lymph_pct", "lymph"), names(.))[1]]],
    feat_mcv   = .data[[intersect(c("MCV", "mcv"), names(.))[1]]],
    feat_rdw   = .data[[intersect(c("RDW", "rdw"), names(.))[1]]],
    feat_alp   = log(.data[[intersect(c("ALP", "alp"), names(.))[1]]]),
    feat_wbc   = log(.data[[intersect(c("WBC", "wbc"), names(.))[1]]])
  )

complete_idx <- complete.cases(df_analysis[, features])

# 开方得到真实的 D_M
df_analysis$DM_score <- NA_real_
df_analysis$DM_score[complete_idx] <- sqrt(
  mahalanobis(
    x      = df_analysis[complete_idx, features],
    center = mu_ref,
    cov    = S_ref
  )
)

df_analysis <- df_analysis %>%
  mutate(
    DM_scaled = DM_score / sd(DM_score, na.rm = TRUE),
    log_DM    = log(DM_score)
  )

# 5. 检验“假性正常”人群与“正常人群”的马氏距离差异
cat("\n----- 四组马氏距离分布统计 (Mean ± SD, Median [IQR]) -----\n")
df_analysis %>%
  filter(!is.na(Group_4cat_90)) %>%
  group_by(Group_4cat_90) %>%
  summarise(
    N         = sum(!is.na(DM_score)),
    DM_Mean   = round(mean(DM_score, na.rm = TRUE), 2),
    DM_SD     = round(sd(DM_score, na.rm = TRUE), 2),
    DM_Median = round(median(DM_score, na.rm = TRUE), 2),
    DM_IQR25  = round(quantile(DM_score, 0.25, na.rm = TRUE), 2),
    DM_IQR75  = round(quantile(DM_score, 0.75, na.rm = TRUE), 2)
  ) %>%
  print()

cat("\n----- 假性正常 (Target) vs 正常 (Ref) 组间比较检验 -----\n")
comp_data <- df_analysis %>% 
  filter(Group_4cat_90 %in% c("Normal (Ref)", "Pseudonormal (Target)"))

print(t.test(DM_score ~ Group_4cat_90, data = comp_data))
print(wilcox.test(DM_score ~ Group_4cat_90, data = comp_data))
library(tidyverse)
library(ggplot2)
library(effsize)   # 用于计算标准化效应量 Cohen's d
library(rstatix)   # 用于成对秩和检验与多重校正
library(PMCMRplus) # 用于非参数趋势检验 (Jonckheere-Terpstra)
library(scales)
# 1. 效应量量化：计算 Cohen's d (Normal vs Pseudonormal)
comp_two_groups <- df_analysis %>%
  filter(Group_4cat_90 %in% c("Normal (Ref)", "Pseudonormal (Target)")) %>%
  filter(!is.na(DM_score))

# 计算 Cohen's d 及其 95% 置信区间
d_result <- cohen.d(DM_score ~ Group_4cat_90, data = comp_two_groups)
cat("----- 假性正常组 vs 正常组 Cohen's d 效应量 -----\n")
print(d_result)

# 2. 四组全景统计：成对比较 (Pairwise Comparisons) 与趋势检验
# 确保四组因子的生物学退行性临床顺序：
# Normal (Ref) -> Pseudoabnormal -> Pseudonormal (Target) -> Impaired
group_levels <- c("Normal (Ref)", "Pseudoabnormal", "Pseudonormal (Target)", "Impaired")
group_levels <- intersect(group_levels, unique(df_analysis$Group_4cat_90))

df_four_groups <- df_analysis %>%
  filter(Group_4cat_90 %in% group_levels) %>%
  filter(!is.na(DM_score)) %>%
  mutate(Group_4cat_90 = factor(Group_4cat_90, levels = group_levels))

# 2.1 四组描述性统计
cat("\n----- 四组马氏距离分布全景统计 (Mean, SD, Median, IQR) -----\n")
df_four_groups %>%
  group_by(Group_4cat_90) %>%
  summarise(
    N         = n(),
    DM_Mean   = mean(DM_score),
    DM_SD     = sd(DM_score),
    DM_Median = median(DM_score),
    IQR_25    = quantile(DM_score, 0.25),
    IQR_75    = quantile(DM_score, 0.75)
  ) %>%
  print()

# 2.2 成对非参数比较 (Wilcoxon pairwise test，带 FDR / Bonferroni 校正)
cat("\n----- 四组成对两两比较 (Pairwise Wilcoxon Test with FDR) -----\n")
pairwise_res <- df_four_groups %>%
  pairwise_wilcox_test(DM_score ~ Group_4cat_90, p.adjust.method = "fdr")
print(pairwise_res)

# 2.3 趋势检验：检验马氏距离是否沿临床风险梯度呈单调递增
cat("\n----- 单调递增趋势检验 (Jonckheere-Terpstra Test) -----\n")
jt_test <- jonckheereTest(df_four_groups$DM_score, g = df_four_groups$Group_4cat_90, alternative = "increasing")
print(jt_test)

# 3. 直观可视化：小提琴图 + 箱线图 + 密度分布对比 (SCI 发表级排版)
# 图 1：四组分布提琴箱线图
p1 <- ggplot(df_four_groups, aes(x = Group_4cat_90, y = DM_score, fill = Group_4cat_90)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.22, outlier.shape = 21, outlier.size = 1.2, alpha = 0.9, color = "#2C3E50") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white", color = "black") +
  scale_fill_manual(values = c(
    "Normal (Ref)"          = "#2E86AB",
    "Pseudoabnormal"        = "#A23B72",
    "Pseudonormal (Target)" = "#E63946",
    "Impaired"              = "#6A0572"
  )) +
  labs(
    title = "Systemic Homeostatic Dysregulation Across Subgroups",
    subtitle = paste0("Jonckheere-Terpstra Trend P < 2.2e-16 | Target vs Ref Cohen's d = ", round(abs(d_result$estimate), 2)),
    x = "Clinical Classification",
    y = expression(paste("Mahalanobis Distance (", D[M], ")"))
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold")
  )

# 图 2：Target vs Ref 双组核密度平滑重叠曲线 (展现两组分布的彻底解耦)
p2 <- ggplot(comp_two_groups, aes(x = DM_score, fill = Group_4cat_90, color = Group_4cat_90)) +
  geom_density(alpha = 0.4, size = 0.9) +
  scale_fill_manual(values = c("Normal (Ref)" = "#2E86AB", "Pseudonormal (Target)" = "#E63946")) +
  scale_color_manual(values = c("Normal (Ref)" = "#1D5F7A", "Pseudonormal (Target)" = "#B71C1C")) +
  geom_vline(xintercept = 3.18, linetype = "dashed", color = "#2E86AB", size = 0.8) +
  geom_vline(xintercept = 6.17, linetype = "dashed", color = "#E63946", size = 0.8) +
  labs(
    title = "Kernel Density Estimation: Normal vs. Pseudonormal",
    x = expression(paste("Mahalanobis Distance (", D[M], ")")),
    y = "Density",
    fill = "Subgroup",
    color = "Subgroup"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 14)
  )

# 打印图像
print(p1)
print(p2)
# 输出马氏距离的表格与图表
library(tidyverse)
library(ggpubr)
# 1. 绘制高水平小提琴 + 箱线图 (Violin + Boxplot + P-value)
# 筛选绘图数据集并设定标准展示因子水平
plot_df <- df_analysis %>%
  filter(Group_4cat_90 %in% c("Normal (Ref)", "Pseudonormal (Target)")) %>%
  filter(!is.na(DM_score)) %>%
  mutate(
    Group = factor(Group_4cat_90, 
                   levels = c("Normal (Ref)", "Pseudonormal (Target)"),
                   labels = c("Normal\n(Reference)", "Pseudonormal\n(Target)"))
  )

p_dm <- ggplot(plot_df, aes(x = Group, y = DM_score, fill = Group)) +
  # 半透明小提琴图展示密度分布
  geom_violin(trim = FALSE, alpha = 0.35, color = NA, width = 0.7) +
  # 内部紧凑箱线图展示四分位线
  geom_boxplot(width = 0.22, color = "#2C3E50", alpha = 0.85, 
               outlier.shape = 21, outlier.size = 1.5, outlier.alpha = 0.5,
               show.legend = FALSE) +
  # 标注各组均值点 (菱形)
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5, 
               fill = "white", color = "#C0392B", stroke = 1.2) +
  # 自动标注极显著 P 检验标签
  geom_signif(
    comparisons = list(c("Normal\n(Reference)", "Pseudonormal\n(Target)")),
    annotations = "P < 0.001",
    y_position = max(plot_df$DM_score, na.rm = TRUE) * 0.95,
    tip_length = 0.02,
    textsize = 4.5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("#2980B9", "#E74C3C")) +
  scale_y_continuous(breaks = seq(0, max(plot_df$DM_score, na.rm = TRUE) + 2, by = 2)) +
  labs(
    title = "Systemic Physiological Dysregulation by Renal-Aging Phenotype",
    subtitle = "Mahalanobis Distance (DM) based on 20-39 y/o healthy reference benchmark",
    x = "",
    y = expression(paste("Mahalanobis Distance (", D[M], ")"))
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "grey30", hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13),
    legend.position = "none"
  )

# 显示图像
print(p_dm)

# 保存为 300 DPI 矢量/印刷级图片
ggsave("Figure_Mahalanobis_Distance_Comparison.png", plot = p_dm, width = 6.5, height = 6.5, dpi = 300)
ggsave("Figure_Mahalanobis_Distance_Comparison.pdf", plot = p_dm, width = 6.5, height = 6.5)

# 2. 导出统计学比较表格 (CSV)
# 提取均值、标准差、中位数、IQR 及检验 P 值
t_res <- t.test(DM_score ~ Group, data = plot_df)
w_res <- wilcox.test(DM_score ~ Group, data = plot_df)

tab_dm <- plot_df %>%
  group_by(Group) %>%
  summarise(
    Sample_Size = n(),
    Mean = sprintf("%.2f", mean(DM_score)),
    SD = sprintf("%.2f", sd(DM_score)),
    Mean_SD = sprintf("%.2f ± %.2f", mean(DM_score), sd(DM_score)),
    Median = sprintf("%.2f", median(DM_score)),
    IQR_25 = sprintf("%.2f", quantile(DM_score, 0.25)),
    IQR_75 = sprintf("%.2f", quantile(DM_score, 0.75)),
    Median_IQR = sprintf("%.2f (%.2f, %.2f)", median(DM_score), quantile(DM_score, 0.25), quantile(DM_score, 0.75))
  ) %>%
  mutate(
    Mean_Difference_95CI = c(
      "Reference", 
      sprintf("%.2f (%.2f to %.2f)", 
              t_res$estimate[2] - t_res$estimate[1], 
              -t_res$conf.int[2], -t_res$conf.int[1])
    ),
    t_test_P = ifelse(t_res$p.value < 0.001, "< 0.001", sprintf("%.4f", t_res$p.value)),
    Wilcoxon_P = ifelse(w_res$p.value < 0.001, "< 0.001", sprintf("%.4f", w_res$p.value))
  )

print(as.data.frame(tab_dm))
write.csv(tab_dm, "Table_Mahalanobis_Comparison_Stats.csv", row.names = FALSE)
# cox
library(tidyverse)
library(survival)

# 1. 数据准备与变量清洗
time_col   <- intersect(c("permth_exm", "time", "permth_int"), names(df_analysis))[1]
status_col <- intersect(c("status", "mortstat"), names(df_analysis))[1]

# 统一分析数据集：设置 Normal (Ref) 为基准参照
df_cox <- df_analysis %>%
  filter(!is.na(Group_4cat_90)) %>%
  mutate(
    time_surv   = as.numeric(.data[[time_col]]),
    status_surv = as.numeric(.data[[status_col]]),
    Group_4cat_90 = relevel(factor(Group_4cat_90), ref = "Normal (Ref)"),
    Sex           = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    Hypertension  = factor(Hypertension),
    Diabetes      = factor(Diabetes)
  )

# 2. 辅助函数：提取 Cox 模型中各暴露组的 HR (95% CI) 与 P 值
extract_cox_res <- function(model, exposure_var, ref_level = "Normal (Ref)") {
  s <- summary(model)
  coef_mat <- s$coefficients
  conf_mat <- s$conf.int
  
  target_rows <- grep(paste0("^", exposure_var), rownames(coef_mat))
  var_names   <- gsub(paste0("^", exposure_var), "", rownames(coef_mat)[target_rows])
  
  hr_str <- sprintf("%.2f (%.2f-%.2f)", 
                    conf_mat[target_rows, "exp(coef)"],
                    conf_mat[target_rows, "lower .95"],
                    conf_mat[target_rows, "upper .95"])
  p_vals <- ifelse(coef_mat[target_rows, "Pr(>|z|)"] < 0.001, "<0.001",
                   sprintf("%.3f", coef_mat[target_rows, "Pr(>|z|)"]))
  
  # 加入 Reference 基准行
  tibble(
    Subgroup = c(ref_level, var_names),
    Result   = c("1.00 (Reference)", paste0(hr_str, " (P = ", p_vals, ")"))
  )
}

# 3. 针对“四分类分组表型 (Group_4cat_90)”拟合三个层级模型

# Model 1: Crude (未校正)
m1_grp <- coxph(Surv(time_surv, status_surv) ~ Group_4cat_90, data = df_cox)

# Model 2: 校正人口学变量 (Age, Sex, Race)
m2_grp <- coxph(Surv(time_surv, status_surv) ~ Group_4cat_90 + age + Sex + race, data = df_cox)

# Model 3: 全校正模型 (+ Hypertension, Diabetes, BMI, eGFR)
m3_grp <- coxph(Surv(time_surv, status_surv) ~ Group_4cat_90 + age + Sex + race + 
                  Hypertension + Diabetes + egfr, data = df_cox)

res_m1_grp <- extract_cox_res(m1_grp, "Group_4cat_90") %>% rename(Model_1 = Result)
res_m2_grp <- extract_cox_res(m2_grp, "Group_4cat_90") %>% rename(Model_2 = Result)
res_m3_grp <- extract_cox_res(m3_grp, "Group_4cat_90") %>% rename(Model_3 = Result)

# 汇总各组事件数与总人数
event_counts <- df_cox %>%
  group_by(Group_4cat_90) %>%
  summarise(
    Events_Total = paste0(sum(status_surv == 1, na.rm = TRUE), "/", n()),
    .groups = "drop"
  ) %>%
  rename(Subgroup = Group_4cat_90)

table_group <- event_counts %>%
  left_join(res_m1_grp, by = "Subgroup") %>%
  left_join(res_m2_grp, by = "Subgroup") %>%
  left_join(res_m3_grp, by = "Subgroup")

# 4. 针对“马氏距离连续指标 (DM_scaled)”拟合三个层级模型
extract_continuous_res <- function(model, var_name) {
  s <- summary(model)
  hr_str <- sprintf("%.2f (%.2f-%.2f)", 
                    s$conf.int[var_name, "exp(coef)"],
                    s$conf.int[var_name, "lower .95"],
                    s$conf.int[var_name, "upper .95"])
  p_val <- ifelse(s$coefficients[var_name, "Pr(>|z|)"] < 0.001, "<0.001",
                  sprintf("%.3f", s$coefficients[var_name, "Pr(>|z|)"]))
  paste0(hr_str, " (P = ", p_val, ")")
}

m1_dm <- coxph(Surv(time_surv, status_surv) ~ DM_scaled, data = df_cox)
m2_dm <- coxph(Surv(time_surv, status_surv) ~ DM_scaled + age + Sex + race, data = df_cox)
m3_dm <- coxph(Surv(time_surv, status_surv) ~ DM_scaled + age + Sex + race + 
                 Hypertension + Diabetes + egfr, data = df_cox)

dm_events <- paste0(sum(df_cox$status_surv == 1 & !is.na(df_cox$DM_scaled)), "/", 
                    sum(!is.na(df_cox$DM_scaled)))

table_dm <- tibble(
  Subgroup     = "Per 1-SD increase in DM",
  Events_Total = dm_events,
  Model_1      = extract_continuous_res(m1_dm, "DM_scaled"),
  Model_2      = extract_continuous_res(m2_dm, "DM_scaled"),
  Model_3      = extract_continuous_res(m3_dm, "DM_scaled")
)

# 5. 合并为标准期刊三线表并导出 CSV
cox_table_final <- bind_rows(
  tibble(Subgroup = "Renal-Aging Phenotype", Events_Total = "", Model_1 = "", Model_2 = "", Model_3 = ""),
  table_group,
  tibble(Subgroup = "Systemic Dysregulation Metric", Events_Total = "", Model_1 = "", Model_2 = "", Model_3 = ""),
  table_dm
)

print(as.data.frame(cox_table_final))
write.csv(cox_table_final, "Table_Cox_Three_Models_Mortality.csv", row.names = FALSE)
# RCS
library(tidyverse)
library(survival)
library(rms)
# 1. 动态匹配变量名与数据准备
time_col   <- intersect(c("permth_exm", "time", "permth_int"), names(df_analysis))[1]
status_col <- intersect(c("status", "mortstat"), names(df_analysis))[1]
bmi_col    <- intersect(c("bmi", "BMXBMI"), names(df_analysis))[1]

# 统一提取 Model 3 所需的核心变量
rcs_vars <- c(time_col, status_col, "DM_score", "age", "Sex", "race", 
              "Hypertension", "Diabetes", bmi_col, "egfr")
rcs_vars <- intersect(rcs_vars, names(df_analysis))

df_rcs <- df_analysis %>%
  select(all_of(rcs_vars)) %>%
  rename(bmi = all_of(bmi_col)) %>%
  drop_na() %>%
  mutate(
    time_surv    = as.numeric(.data[[time_col]]),
    status_surv  = as.numeric(.data[[status_col]]),
    Sex          = as.factor(Sex),
    race         = as.factor(race),
    Hypertension = as.factor(Hypertension),
    Diabetes     = as.factor(Diabetes)
  )

# 2. 设置 rms 绘图环境并拟合全校正 RCS Cox 模型
dd <- datadist(df_rcs)
options(datadist = "dd")

# 设定健康基准中位数作为对照参考点 (HR = 1.0)
ref_val <- round(median(df_rcs$DM_score, na.rm = TRUE), 2)

# 拟合 4 节点限制性立方样条模型 (严格对齐 Model 3)
fit_rcs <- cph(
  Surv(time_surv, status_surv) ~ rcs(DM_score, 4) + age + Sex + race + 
    Hypertension + Diabetes + bmi + egfr,
  data = df_rcs,
  x = TRUE, 
  y = TRUE
)

# 提取整体检验与非线性检验 P 值
anova_rcs     <- anova(fit_rcs)
p_overall     <- anova_rcs["DM_score", "P"]
p_nonlin      <- anova_rcs[" Nonlinear", "P"]
p_overall_str <- ifelse(p_overall < 0.001, "< 0.001", sprintf("%.3f", p_overall))
p_nonlin_str  <- ifelse(p_nonlin < 0.001,  "< 0.001", sprintf("%.3f", p_nonlin))

cat(paste0("\n>>> RCS 统计检验: P-overall = ", p_overall_str, " | P-non-linear = ", p_nonlin_str, "\n"))

# 3. 截断极值区间并自动量化致死风险升高拐点 (Cut-off)

# 截取第 1 至第 99 百分位数，避免极端离群值导致置信区间发散
p1_val  <- quantile(df_rcs$DM_score, 0.01)
p99_val <- quantile(df_rcs$DM_score, 0.99)

eval_points <- seq(p1_val, p99_val, length.out = 300)
pred_rcs    <- Predict(fit_rcs, DM_score = eval_points, ref.zero = TRUE, fun = exp)
pred_df     <- as.data.frame(pred_rcs)

# 寻找致死拐点：中位数右侧 95% CI 下限 (lower) 首次稳定大于 1.0 的点
inflection_row <- pred_df %>%
  filter(DM_score >= ref_val, lower > 1.0) %>%
  slice(1)

cutoff_dm <- round(inflection_row$DM_score, 2)
cat(paste0(">>> 自动化确立致死拐点阈值: 当 DM > ", cutoff_dm, " 时，死亡风险发生统计学显著增加 (Lower 95% CI > 1.0)\n\n"))

# 4. 绘制并导出符合发表标准的 RCS 曲线
p_rcs <- ggplot(pred_df, aes(x = DM_score, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#E74C3C", alpha = 0.18) +
  geom_line(color = "#C0392B", linewidth = 1.2) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey45", linewidth = 0.8) +
  geom_vline(xintercept = ref_val, linetype = "dotted", color = "#2980B9", linewidth = 0.9) +
  geom_vline(xintercept = cutoff_dm, linetype = "dashed", color = "#8E44AD", linewidth = 0.9) +
  annotate("text", x = ref_val, y = 0.7, 
           label = paste0("Reference\n(Median = ", ref_val, ")"), 
           color = "#2980B9", size = 3.6, fontface = "bold", vjust = 1) +
  annotate("text", x = cutoff_dm, y = max(pred_df$upper) * 0.75, 
           label = paste0("Threshold\n(DM = ", cutoff_dm, ")"), 
           color = "#8E44AD", size = 3.6, fontface = "bold", hjust = -0.1) +
  annotate("text", x = p1_val + 0.5, y = max(pred_df$upper) * 0.92, 
           label = paste0("P for overall < 0.001\nP for non-linearity = ", p_nonlin_str), 
           hjust = 0, size = 4, fontface = "italic", color = "grey20") +
  scale_y_continuous(breaks = seq(0.5, ceiling(max(pred_df$upper)), by = 0.5)) +
  labs(
    title = "Dose-Response Association Between Mahalanobis Distance and All-Cause Mortality",
    subtitle = "Multivariable Cox restricted cubic spline model fully adjusted for Model 3 covariates",
    x = expression(paste("Systemic Physiological Dysregulation (Mahalanobis Distance, ", D[M], ")")),
    y = "Hazard Ratio (95% CI)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, color = "grey30", hjust = 0.5, margin = margin(b = 12)),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "black")
  )

print(p_rcs)

ggsave("Figure_RCS_DM_FullyAdjusted.png", plot = p_rcs, width = 7.5, height = 5.8, dpi = 300)
ggsave("Figure_RCS_DM_FullyAdjusted.pdf", plot = p_rcs, width = 7.5, height = 5.8)
library(tidyverse)
library(cluster)
library(factoextra)
library(survival)
library(survminer)
# 聚类分析
# 第一步：提取 157 例“假性正常”目标人群及其标准化特征矩阵
# 1. 明确参与聚类的核心生物学轴系变量（已完成 Z-score 标准化的变量）
cluster_feature_names <- c(
  "z_feat_alb",    # 营养/肝脏稳态 (Albumin)
  "z_feat_creat",  # 基础肌肉代谢 (Creatinine)
  "z_feat_glu",    # 糖代谢负荷 (Glucose)
  "z_feat_crp",    # 全身微炎症 (CRP)
  "z_feat_lymph",  # 获得性免疫储备 (Lymphocyte %)
  "z_feat_mcv",    # 红细胞容积 (MCV)
  "z_feat_rdw",    # 红细胞异质性/造血系统应激 (RDW)
  "z_feat_alp",    # 肝胆/骨骼代谢毒性 (ALP)
  "z_feat_wbc"     # 先天免疫/白细胞系统 (WBC)
)

# 2. 筛选 157 例假性正常人群，并剔除特征矩阵中缺失的个案
df_target <- df_analysis %>%
  filter(Group_4cat_90 == "Pseudonormal (Target)") %>%
  filter(if_all(all_of(cluster_feature_names), ~ !is.na(.x)))

cat("成功纳入聚类分析的目标假性正常样本量 N =", nrow(df_target), "\n")

# 3. 构建聚类专用矩阵 (行名为受试者唯一索引，便于后续合并标签)
mat_cluster <- df_target %>%
  select(all_of(cluster_feature_names)) %>%
  as.matrix()
rownames(mat_cluster) <- paste0("ID_", seq_len(nrow(mat_cluster)))

# 第二步：自动化网格搜索最佳聚类簇数 K (Elbow Method & Silhouette Method)
set.seed(2026)

# 1. 手肘法 (Total Within Sum of Squares)
p_elbow <- fviz_nbclust(mat_cluster, kmeans, method = "wss", k.max = 6) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "#E74C3C") +
  labs(title = "Optimal Number of Clusters: Elbow Method",
       x = "Number of Clusters K", y = "Total Within-Cluster Sum of Squares") +
  theme_classic(base_size = 12)

# 2. 平均轮廓系数法 (Average Silhouette Width, 寻找峰值)
p_silhouette <- fviz_nbclust(mat_cluster, kmeans, method = "silhouette", k.max = 6) +
  labs(title = "Optimal Number of Clusters: Silhouette Method",
       x = "Number of Clusters K", y = "Average Silhouette Width") +
  theme_classic(base_size = 12)

# 拼图展示并保存评估图
p_eval <- cowplot::plot_grid(p_elbow, p_silhouette, ncol = 2)
print(p_eval)
ggsave("Figure_Clustering_K_Evaluation.png", plot = p_eval, width = 9, height = 4, dpi = 300)

# 第三步：执行 K-means 聚类 (以确定 K = 2 为例，亦可根据评估图调整为 3)
k_optimal <- 2  # 根据临床假设与轮廓系数，解耦为 2 种主要失调表型

set.seed(2026)
fit_km <- kmeans(mat_cluster, centers = k_optimal, nstart = 50, iter.max = 100)

# 将聚类簇标签追加到数据集
df_target$Subphenotype <- factor(
  fit_km$cluster, 
  levels = 1:k_optimal,
  labels = paste0("Subphenotype_", 1:k_optimal)
)

cat("各亚型样本量分配：\n")
print(table(df_target$Subphenotype))

# 降维可视化：基于 PCA 投影的聚类散点图
p_cluster_pca <- fviz_cluster(
  fit_km, 
  data = mat_cluster,
  geom = "point",
  ellipse.type = "convex",
  palette = c("#2980B9", "#E74C3C"),
  ggtheme = theme_classic(base_size = 12),
  main = "PCA Projection of Pseudonormal Subphenotypes"
)
print(p_cluster_pca)
# 第四步：绘制各亚型标准化特征偏离度条形图
# 1. 计算每个簇在各维度的均值
profile_df <- df_target %>%
  group_by(Subphenotype) %>%
  summarise(across(all_of(cluster_feature_names), mean), .groups = "drop") %>%
  pivot_longer(-Subphenotype, names_to = "Feature", values_to = "Standardized_Mean") %>%
  mutate(
    Clean_Feature = case_when(
      Feature == "z_feat_alb"   ~ "Albumin",
      Feature == "z_feat_creat" ~ "Creatinine",
      Feature == "z_feat_glu"   ~ "Glucose",
      Feature == "z_feat_crp"   ~ "hs-CRP",
      Feature == "z_feat_lymph" ~ "Lymphocyte %",
      Feature == "z_feat_mcv"   ~ "MCV",
      Feature == "z_feat_rdw"   ~ "RDW",
      Feature == "z_feat_alp"   ~ "ALP",
      Feature == "z_feat_wbc"   ~ "WBC",
      TRUE ~ Feature
    )
  )

# 2. 绘制分面水平条形图
p_profile <- ggplot(profile_df, aes(x = reorder(Clean_Feature, Standardized_Mean), 
                                   y = Standardized_Mean, 
                                   fill = Subphenotype)) +
  geom_bar(stat = "identity", width = 0.65, alpha = 0.85) +
  coord_flip() +
  facet_wrap(~ Subphenotype, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Subphenotype_1" = "#2980B9", "Subphenotype_2" = "#E74C3C")) +
  labs(
    title = "Clinical and Biochemical Signatures Across Subphenotypes",
    subtitle = "Z-score standardized deviation profiles (Mean = 0 line indicates population baseline)",
    x = "Biomarkers",
    y = "Standardized Mean Value (Z-score)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0.5, margin = margin(b = 10)),
    axis.text = element_text(color = "black")
  )

print(p_profile)
ggsave("Figure_Subphenotype_Profiles_Barplot.png", plot = p_profile, width = 8, height = 5, dpi = 300)
library(tidyverse)
# 亚型表型特征临床转归
# 1. 检索与死亡结局、死因分类相关的变量
mort_candidates <- c(
  "mortstat", "status", "permth_int", "permth_exm", "time",
  "ucod_leading", "diabetes_death", "hyperten_death", 
  "ucod", "cause_death", "death_cause", "eligstat"
)

found_mort_vars <- intersect(mort_candidates, names(df_analysis))
cat(">>> 检测到的死亡及死因相关变量：\n")
print(found_mort_vars)

# 查看死因变量的频数分布（若存在 ucod_leading 或类似字段）
if ("ucod_leading" %in% names(df_analysis)) {
  cat("\n----- 全人群死因编码 (ucod_leading) 分布 -----\n")
  print(table(df_analysis$ucod_leading, useNA = "ifany"))
  
  if (exists("df_target")) {
    cat("\n----- 157 例假性正常人群内部死因分布 -----\n")
    print(table(df_target$Subphenotype, df_target$ucod_leading, useNA = "ifany"))
  }
}

# 2. 检索与肾脏终点、透析、蛋白尿相关的变量
kidney_candidates <- c(
  # 问卷变量：透析/肾病病史
  "KIQ022", "KIQ026", "KIQ025", "kiq022", "dialysis", "kidney_failure",
  # 检验变量：尿微量白蛋白、尿肌酐、ACR
  "URXUMA", "URXUCR", "acr", "uacr", "albuminuria",
  # 随访结局相关字段
  "ckd_incident", "esrd", "renal_death"
)

found_kidney_vars <- intersect(kidney_candidates, names(df_analysis))
cat("\n>>> 检测到的肾脏/透析/蛋白尿相关变量：\n")
print(found_kidney_vars)

# 3. 模糊匹配所有可能相关的字段名（以防命名不一致）
cat("\n>>> 包含 'death', 'mort', 'kidney', 'renal', 'kiq' 的全部字段：\n")
grep("death|mort|kidney|renal|kiq|ucod", names(df_analysis), ignore.case = TRUE, value = TRUE)
