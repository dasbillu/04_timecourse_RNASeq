### Contents (last updated: 26 April 2023)

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
## set parameters and thresholds
#


# 00. Databases -----------------------------------------------------------

# 1. TC9_data.db
#
# Desc: This database contains all normalized expression data collected for TC7
#
#
# a. pbar_ld_expressed_genes
# b. pbar_ld_fpkm
# c. pbar_ld_log2fpkm
# d. pbar_ld_zscores
#
# e. pbar_dd_expressed_genes
# f. pbar_dd_fpkm
# g. pbar_dd_log2fpkm
# h. pbar_dd_zscores
#
# Load the database
my.db <- dbConnect(RSQLite::SQLite(),
                   "./data/databases/TC9_data.db")
# which tables are in the database
src_dbi(my.db)
#
# set path to raw FPKM data
path_to_data = "./results/gene_exp/"


##-##-##-##-##-##-##-##-##-##-##-##-
#
# Dataset: LD  ##-##-##-##-##-##-##-
#
##-##-##-##-##-##-##-##-##-##-##-##-

# 01. pbar_ld ----------------------------------------------------------

##-##-##-##-##-##-##-##-##-##-##-##-
# SETUP 
##-##-##-##-##-##-##-##-##-##-##-##-
# path to sample data
folder_name = "LD"
file_name = "normalized_fpkm_pogo_LD.csv"

# read data
pbar.ld <- read.csv(paste0(path_to_data,folder_name,"/",file_name),
                   header = T, stringsAsFactors = F, na.strings = c(NA, "", " "))
# Save file to database
dbWriteTable(my.db, "pbar_ld_fpkm", pbar.ld)


# Expressed genes ---------------------------------------------------------

##-##-##-##-##-##-##-##-##-##-##-##-
# list of all pbar_ld "Expressed" genes 
# (>1 FPKM during the 24h period)
##-##-##-##-##-##-##-##-##-##-##-##-
expressed <-
  pbar.ld %>%
  na.omit() %>%
  filter_at(vars(starts_with("Z")), any_vars(. > 1)) %>%
  pull(gene_name) %>%
  unique()
# append this information to the pbar.ld data
expressed.pbar.ld <-
  pbar.ld %>%
  select(gene_name) %>%
  mutate(expressed = ifelse(gene_name %in% expressed,"yes","no"))
# Save file to databse
dbWriteTable(my.db, "pbar_ld_expressed_genes", expressed.pbar.ld)

# log2-expression --------------------------------------------------------

##-##-##-##-##-##-##-##-##-##-##-##-
# log2-transform the data
##-##-##-##-##-##-##-##-##-##-##-##-
gene.names <- pbar.ld %>% pull(gene_name)
log2.pbar.ld <- log2(pbar.ld[-1] + 1)
log2.pbar.ld$gene_name <- gene.names
log2.pbar.ld <- log2.pbar.ld %>% select(gene_name, everything())
# check the log2-transformed data
log2.pbar.ld %>% head()
# Save file to database
dbWriteTable(my.db, "pbar_ld_log2fpkm", log2.pbar.ld)

# zscore-log2-expression --------------------------------------------------

##-##-##-##-##-##-##-##-##-##-##-##-
# z-score the data
##-##-##-##-##-##-##-##-##-##-##-##-
gene.names <- log2.pbar.ld %>% pull(gene_name)
sample.names <- names(log2.pbar.ld[-1])
zscores.pbar.ld <- log2.pbar.ld %>%
  # create a gene x exp matrix
  select(-1) %>%
  as.matrix() %>%
  # use the scale function for each row to calculate z-scores
    # scale calculates (x-mean(X))/sd(X)
    # 1 indicates row-wise
  apply(., 1, scale) %>%
  # the output needs to be transposed
  t() %>%
  # make it a dataframe
  as.data.frame()
# add the column names
names(zscores.pbar.ld) <- sample.names
zscores.pbar.ld$gene_name <- gene.names
zscores.pbar.ld <- zscores.pbar.ld %>% select(gene_name, everything())
# check the z-score transformed dataset
zscores.pbar.ld %>% head()
# Save file to database
dbWriteTable(my.db, "pbar_ld_zscores", zscores.pbar.ld)

##-##-##-##-##-##-##-##-##-##-##-##-
# Save a csv with the zscores
##-##-##-##-##-##-##-##-##-##-##-##-
pbar.ld.zscore <- tbl(my.db, "pbar_ld_zscores") %>% collect()

pbar.ld.zscore.noNAs <-
  pbar.ld.zscore %>%
    # remove the rows with NAs (ejtk-does-not-like-NAs)
    na.omit()

# format the column headings for ejtk
colnames(pbar.ld.zscore.noNAs)[1] <- "#"
colnames(pbar.ld.zscore.noNAs)[-1] <- paste0("ZT",seq(1,23,2))

# check if all looks good
pbar.ld.zscore.noNAs %>%
# save the file
write.csv(.,
          file = paste0(path_to_data,"../zscore/",folder_name,"/","zscores_noNAs_pbar_LD.csv"),
          row.names = F)




##-##-##-##-##-##-##-##-##-##-##-##-
#
# Dataset: DD  ##-##-##-##-##-##-##-
#
##-##-##-##-##-##-##-##-##-##-##-##-

# 01. pbar_dd ----------------------------------------------------------

##-##-##-##-##-##-##-##-##-##-##-##-
# SETUP 
##-##-##-##-##-##-##-##-##-##-##-##-
# path to sample data
folder_name = "DD"
file_name = "normalized_fpkm_pogo_DD.csv"

# read data
pbar.dd <- read.csv(paste0(path_to_data,folder_name,"/",file_name),
                    header = T, stringsAsFactors = F, na.strings = c(NA, "", " "))
# Save file to database
dbWriteTable(my.db, "pbar_dd_fpkm", pbar.dd)


# Expressed genes ---------------------------------------------------------

##-##-##-##-##-##-##-##-##-##-##-##-
# list of all pbar_dd "Expressed" genes 
# (>1 FPKM during the 24h period)
##-##-##-##-##-##-##-##-##-##-##-##-
expressed <-
  pbar.dd %>%
  na.omit() %>%
  filter_at(vars(starts_with("Z")), any_vars(. > 1)) %>%
  pull(gene_name) %>%
  unique()
# append this information to the pbar.dd data
expressed.pbar.dd <-
  pbar.dd %>%
  select(gene_name) %>%
  mutate(expressed = ifelse(gene_name %in% expressed,"yes","no"))
# Save file to databse
dbWriteTable(my.db, "pbar_dd_expressed_genes", expressed.pbar.dd)

# log2-expression --------------------------------------------------------

##-##-##-##-##-##-##-##-##-##-##-##-
# log2-transform the data
##-##-##-##-##-##-##-##-##-##-##-##-
gene.names <- pbar.dd %>% pull(gene_name)
log2.pbar.dd <- log2(pbar.dd[-1] + 1)
log2.pbar.dd$gene_name <- gene.names
log2.pbar.dd <- log2.pbar.dd %>% select(gene_name, everything())
# check the log2-transformed data
log2.pbar.dd %>% head()
# Save file to database
dbWriteTable(my.db, "pbar_dd_log2fpkm", log2.pbar.dd)

# zscore-log2-expression --------------------------------------------------

##-##-##-##-##-##-##-##-##-##-##-##-
# z-score the data
##-##-##-##-##-##-##-##-##-##-##-##-
gene.names <- log2.pbar.dd %>% pull(gene_name)
sample.names <- names(log2.pbar.dd[-1])
zscores.pbar.dd <- log2.pbar.dd %>%
  # create a gene x exp matrix
  select(-1) %>%
  as.matrix() %>%
  # use the scale function for each row to calculate z-scores
  # scale calculates (x-mean(X))/sd(X)
  # 1 indicates row-wise
  apply(., 1, scale) %>%
  # the output needs to be transposed
  t() %>%
  # make it a dataframe
  as.data.frame()
# add the column names
names(zscores.pbar.dd) <- sample.names
zscores.pbar.dd$gene_name <- gene.names
zscores.pbar.dd <- zscores.pbar.dd %>% select(gene_name, everything())
# check the z-score transformed dataset
zscores.pbar.dd %>% head()
# Save file to database
dbWriteTable(my.db, "pbar_dd_zscores", zscores.pbar.dd)

##-##-##-##-##-##-##-##-##-##-##-##-
# Save a csv with the zscores
##-##-##-##-##-##-##-##-##-##-##-##-
pbar.dd.zscore <- tbl(my.db, "pbar_dd_zscores") %>% collect()

pbar.dd.zscore.noNAs <-
  pbar.dd.zscore %>%
  # remove the rows with NAs (ejtk-does-not-like-NAs)
  na.omit()

# format the column headings for ejtk
colnames(pbar.dd.zscore.noNAs)[1] <- "#"
colnames(pbar.dd.zscore.noNAs)[-1] <- paste0("ZT",seq(1,23,2))

# check if all looks good
pbar.dd.zscore.noNAs %>%
  # save the file
  write.csv(.,
            file = paste0(path_to_data,"../zscore/",folder_name,"/","zscores_noNAs_pbar_DD.csv"),
            row.names = F)



