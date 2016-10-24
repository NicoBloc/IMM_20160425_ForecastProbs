library(gdata)
library(dplyr)

d <- read.xls("../../Data/QC_Disk_Diffusion_2015_Ecol_ATCC_25922_Sirweb.xlsx")

d$lot <- c(rep('131224A', 23), rep('140818C', 18), rep('150205A', 14))

f <- select(d, TOB, lot) %>%
  group_by(lot) %>%
  mutate(a=(TOB - mean(TOB))) %>%
  ungroup

print(sd(f$a))
