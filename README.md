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
