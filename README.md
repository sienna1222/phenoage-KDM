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
