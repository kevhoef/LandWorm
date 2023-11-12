library(readxl)
gp_oasys_rd=read_xlsx("./[Database]/[Raw datasets]/LandWorm data providers/Oasys GP/OASYS 2021 - sauvegarde 03.02.2022.xlsx",sheet="data pour plugin LMCU")
write.csv(gp_oasys_rd, "[Database]/[Raw datasets]/EW_datasets/gp_oasys_rd.csv", row.names = FALSE)
