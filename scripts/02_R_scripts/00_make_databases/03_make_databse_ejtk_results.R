
# Housekeeping ---------------------------------------------------------------
#
set.seed(420)
rm(list = ls())
#
## Load packages ----------
pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)
#
# set conflict preference (housekeeping to make sure functions work as expected)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")
#
## set parameters and thresholds --------
#
# gamma-pvalue threshold for inferring rhythmicity
gamma.pval = 0.05
#
## Path to save files ------
# # for supplementary files
# supp.path="~/University\ of\ Central\ Florida/Charissa\ De\ Bekker\ -\ Ant-Fungus-Clock-Interactions/04_manuscript/03_supplementary_files/"


# 01. Load database -----------------------------------------------------------

# 1. TC9_ejtk.db
#
# Desc: This database contains all ejtk-output for TC7
#
### Contents (last updated 19-Sep-21)
#
# a. pbar_ld_zscores_24h
# b. pbar_ld_zscores_12h
# c. pbar_ld_zscores_08h
#
# d. pbar_dd_zscores_24h
# e. pbar_dd_zscores_12h
# f. pbar_dd_zscores_08h
#
# Load the data
my.db <- dbConnect(RSQLite::SQLite(),
                   "./data/databases/TC9_ejtk.db")
# which tables are in the database
src_dbi(my.db)
#
# set path to eJTK output data
path_to_data = "./results/ejtk_output/"


# 01. pbar_ld ----------------------------------------------------------
#
## Period = 24h
#
# path to sample data
folder_name = "LD"
file_name = "zscores_noNAs_pbar_LD_cos24_ph0022by2_as0222by2_zscore_jtkout_GammaP.txt"
#
pbar.ld.24.zscore <- read.csv(paste0(path_to_data,folder_name,"/",file_name),
                          sep = "\t", header = T, stringsAsFactors = F)
pbar.ld.24.zscore %>%
  head()
# Save file to database
dbWriteTable(my.db, "pbar_ld_zscores_24h", pbar.ld.24.zscore)
#
#
## Period = 12h
#
# path to sample data
file_name = "zscores_noNAs_pbar_LD_cos12_ph0022by2_as0222by2_zscore_jtkout_GammaP.txt"
#
pbar.ld.12.zscore <- read.csv(paste0(path_to_data,folder_name,"/",file_name),
                          sep = "\t", header = T, stringsAsFactors = F)
pbar.ld.12.zscore %>%
  head()
# Save file to database
dbWriteTable(my.db, "pbar_ld_zscores_12h", pbar.ld.12.zscore)
#
#
## Period = 8h
#
# path to sample data
file_name = "zscores_noNAs_pbar_LD_cos08_ph0022by2_as0222by2_zscore_jtkout_GammaP.txt"
#
pbar.ld.08.zscore <- read.csv(paste0(path_to_data,folder_name,"/",file_name),
                          sep = "\t", header = T, stringsAsFactors = F)
pbar.ld.08.zscore %>%
  head()
# Save file to database
dbWriteTable(my.db, "pbar_ld_zscores_08h", pbar.ld.08.zscore)
#
#
# rm(pbar.ld.08.zscore, pbar.ld.12.zscore, pbar.ld.24.zscore)
## End.


# 02. pbar_dd ----------------------------------------------------------
#
## Period = 24h
#
# path to sample data
folder_name = "DD"
file_name = "zscores_noNAs_pbar_DD_cos24_ph0022by2_as0222by2_zscore_jtkout_GammaP.txt"
#
pbar.dd.24.zscore <- read.csv(paste0(path_to_data,folder_name,"/",file_name),
                                   sep = "\t", header = T, stringsAsFactors = F)
pbar.dd.24.zscore %>%
  head()
# Save file to database
dbWriteTable(my.db, "pbar_dd_zscores_24h", pbar.dd.24.zscore, overwrite=T)
#
#
## Period = 12h
#
# path to sample data
file_name = "zscores_noNAs_pbar_DD_cos12_ph0022by2_as0222by2_zscore_jtkout_GammaP.txt"
#
pbar.dd.12.zscore <- read.csv(paste0(path_to_data,folder_name,"/",file_name),
                                   sep = "\t", header = T, stringsAsFactors = F)
pbar.dd.12.zscore %>%
  head()
# Save file to database
dbWriteTable(my.db, "pbar_dd_zscores_12h", pbar.dd.12.zscore, overwrite=T)
#
#
## Period = 8h
#
# path to sample data
file_name = "zscores_noNAs_pbar_DD_cos08_ph0022by2_as0222by2_zscore_jtkout_GammaP.txt"
#
pbar.dd.08.zscore <- read.csv(paste0(path_to_data,folder_name,"/",file_name),
                                   sep = "\t", header = T, stringsAsFactors = F)
pbar.dd.08.zscore %>%
  head()
# Save file to database
dbWriteTable(my.db, "pbar_dd_zscores_08h", pbar.dd.08.zscore, overwrite=T)
#
# rm(pbar.dd.08.zscore, pbar.dd.12.zscore, pbar.dd.24.zscore)
#
## End.
