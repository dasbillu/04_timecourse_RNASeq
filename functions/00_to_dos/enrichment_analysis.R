# Enrichment test (GO and PFAM)


#' # 1. GO enrichment - Function --------------------------------------------------
#' go_enrichment <- function(geneset,
#'                           function.dir = ".",
#'                           data.dir = ".",
#'                           org = "cflo", 
#'                           bg = "all", 
#'                           atleast = 5, 
#'                           enriched.terms="over") {
#' 
#'   # save the input list of genes for enrichment test
#'   genes <- geneset
#' 
#'   ## Load the required libraries
#'   library(tidyverse)
#' 
#'   ## load the selected annotation file
#'     if (org=="ophio_cflo"){
#'       
#'       print("Loading annotation file for Ophiocordyceps camponoti-floridani")
#'       all_genes <- read.csv(paste0(function.dir,"/functions/func_data/ophio_cflo_annots_robin_ncbi.csv"), 
#'                             header=T, stringsAsFactors = F, na.strings = c(NA,""," "))
#'       # define the separator
#'       separator = "; "
#'       
#'       print("Done.")
#'       
#'     } else if (org=="ophio_kim"){
#'       
#'       print("Loading annotation file for Ophiocordyceps kimflemingae")
#'       all_genes <- read.csv(paste0(function.dir,"/functions/func_data/ophio_kim_annots_robin_ncbi.csv"), 
#'                             header=T, stringsAsFactors = F, na.strings = c(NA,""," "))
#'       # define the separator
#'       separator = ";"
#'       
#'       print("Done.")
#'       
#'     } else if (org=="cflo"){
#'       
#'       print("Loading annotation file for Camponotus floridanus")
#'       all_genes <- read.csv(paste0(function.dir,"/functions/func_data/cflo_annots.csv"), 
#'                             header=T, stringsAsFactors = F, na.strings = c(NA,""," "))
#'       # define the separator
#'       separator = "; "
#'       
#'       print("Done.")
#'       
#'     } else {
#'       
#'       print("Invalid option for argument org.")
#'       print("Select from: ophio_cflo, ophio_kim, cflo")
#'       stop()
#'     }
#'   
#'   ## make the flattened gene x annotation file
#'     # Select only the columns that we need
#'     all_genes <- all_genes[,c("gene_name","GOs")]
#'     #' Let's replace the NAs in GOs and pfams with "no_annot"
#'     all_genes[is.na(all_genes)] <- "no_annot"
#'     #' Let's flatten the files
#'     all_genes_gos <- all_genes %>%
#'       dplyr::mutate(go_split = str_split(GOs, separator)) %>%
#'       unnest() %>%
#'       #dplyr::select(-GO) %>%
#'       separate(go_split, c("GO","GO_desc"), sep = "([\\|])", extra = "drop") %>%
#'       dplyr::select(gene_name, GO, GO_desc) %>%
#'       unique()
#'   
#'   ## define the background geneset to run enrichment against
#'     if(bg == "all") {
#'       # save the background data frame with gene_name, GO term, and GO description in an object
#'       background <- all_genes_gos %>%
#'         arrange(gene_name)
#'       
#'     } else if (bg == "expressed") {
#'       
#'       ## If-else statement to load the expressed geneset for the organism
#'       
#'       if (org=="ophio_cflo"){
#'         
#'         foo <- tbl(dbConnect(RSQLite::SQLite(),paste0(data.dir,"/data/databases/TC6_fungal_data.db")), 
#'                    "ophio_cflo_expressed_genes") %>% 
#'                       filter(expressed=="yes") %>% 
#'                       collect() %>% pull(gene_name) %>% as.character()
#'         background <- all_genes_gos %>%
#'           # filter and keep user specified background geneset 
#'           filter(gene_name %in% foo) %>%
#'           arrange(gene_name)
#'       }
#'       
#'     } else (
#'       background <- all_genes_gos %>%
#'         # filter and keep user specified background geneset 
#'         filter(gene_name %in% as.character(bg)) %>%
#'         arrange(gene_name)
#'     )
#'   
#'   ## Need a GO term to GO description file
#'     go_to_desc <- dplyr::distinct(as.data.frame(all_genes_gos[-1]))
#' 
#'   ## Enrichment to be tested for all GO terms that are:
#'     ## 1. Present in the test geneset
#'     ## 2. all unique GO terms each present in at least x number of genes
#'       annot_terms <- background %>%
#'         # Keep only the genes in my test geneset
#'         filter(gene_name %in% genes) %>%
#'         group_by(GO) %>%
#'         summarize(num_genes = n()) %>%
#'         arrange(num_genes) %>%
#'         # Keep only the GO terms that are annotated in at least 5 genes
#'         filter(num_genes >= atleast) %>%
#'         pull(GO)
#' 
#'   # Let's get all the genes in our geneset for each GO term
#'   # Test geneset dataframe
#'   df.test <- background %>%
#'     filter(gene_name %in% genes)
#'   go_to_genes <- aggregate( .~ GO, df.test, function(x) toString(unique(x)))
#' 
#'   ## Make an empty list that can save your results for each GO term
#'   df.list <- list()
#' 
#'   ## Print the summary stats for the enrichment test
#'   print(paste0("Number of genes in background geneset: ", background %>% distinct(gene_name) %>% nrow()))
#'   print(paste0("Number of genes in the test set: ", length(genes)))
#'   print("--------------------------------")
#'   print(paste0("Number of GO terms in background geneset: ", background %>% distinct(GO) %>% nrow()))
#'   print(paste0("Number of GO terms (at least ", atleast, "genes) in background geneset: ", background %>% group_by(GO) %>% summarise(num_genes = n()) %>% filter(num_genes >= atleast) %>% nrow()))
#'   print(paste0("Number of GO terms (at least ", atleast, "genes) in test set: ",length(annot_terms)))
#' 
#'   if(length(annot_terms) == 0) {
#'     return(as.character("There are no GO terms to test enrichment for."))
#'   } else if (length(annot_terms) >= 1)
#'   {
#'     print("Testing for enrichment...")
#'     # Test the enrichment for each of the GO terms
#'     for (i in 1:length(annot_terms))
#'     {
#'       # get the GO term to be tested for enrichment
#'       annot <- annot_terms[i]
#' 
#'       # number of DEGs (or genes of interest)
#'       n_DEG <- length(genes)
#' 
#'       # Number of genes annotated with the GO term in the background gene set
#'       n_GO <- background %>%
#'         filter(GO == annot) %>%
#'         nrow()
#' 
#'       # Number of genes NOT annotated with the GO term in the background gene set
#'       n_not_GO <- length(unique(background[[1]])) - n_GO
#' 
#'       # Number of genes annotated with the GO term in the test set
#'       n_GO_DEG <-  df.test %>%
#'         filter(GO == annot) %>%
#'         nrow()
#' 
#'       pval <- dhyper(n_GO_DEG,
#'                      n_GO,
#'                      n_not_GO,
#'                      n_DEG, log=F)
#' 
#'       df.list[[i]] <- data.frame(GO = annot_terms[i],
#'                                  sam_freq = round(n_GO_DEG/n_DEG, 5),
#'                                  back_freq = round(n_GO/sum(n_GO, n_not_GO), 5),
#'                                  n_GO_DEG = n_GO_DEG,  # x
#'                                  n_DEG = n_DEG,      # k
#'                                  n_GO = n_GO,       # m
#'                                  #n_not_GO = n_not_GO,  # n
#'                                  n_all_genes = sum(n_GO,n_not_GO),  # n_all_genes
#'                                  pVal = pval)
#' 
#' 
#'     }
#' 
#'   ## Make the output table:
#'     df.enriched <- bind_rows(df.list, .id = "column_label") %>%
#'       dplyr::select(-column_label) %>%
#'       arrange(pVal) %>%
#'       mutate(adj_pVal = p.adjust(pVal, "BH")) %>%
#'       # keeps only the GO terms that are found in the test set
#'       filter(n_GO_DEG != 0) %>%
#'       #filter(pVal < 0.1 | adj_pVal < 0.1) %>%   ## I am not filtering anything yet. Return the whole file.
#'       #left_join(go_to_desc, by="GO") %>%
#'       left_join(go_to_genes, by="GO") %>%
#'       mutate(over_under = ifelse(sam_freq > back_freq, "over", "under")) %>%
#'       dplyr::select(GO, GO_desc, over_under, adj_pVal, everything()) %>%
#'       arrange(over_under, adj_pVal)
#' 
#'     # if (terms == "over") {
#'     #   df.enriched <- df.enriched %>%
#'     #     filter(over_under == "over")
#'     # }
#' 
#' 
#'     return(df.enriched);
#'   }
#' 
#'   
#' }

### Usage:
# #Let's test the enichment for some mock geneset
# # Let's load the core dataset to sample some gene names
# load("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_core_datasets.RData")
# mock.genes <- sample(cflo.annots.exp[[1]], 1000)
# mock.enriched.gos <- cflo_go_enrichment(geneset = mock.genes)
# # View the first 5 rows of the results
# mock.enriched.gos[1:5,1:6]

# #Let's try the set of rhythmic genes (RAIN)
# # load datasets and get all the sig. (FDR 5%) rhythmic genes
# load("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_RAIN_datasets.RData")
# for.rhy.genes <- for.rain.atleast4fpkm %>% 
#   mutate(gene_name = rownames(for.rain.atleast4fpkm)) %>% 
#   filter(pVal < 0.05) %>% 
#   pull(gene_name)
# # call the enrichment function
# for.rhy.enriched.gos <- cflo_go_enrichment(geneset = for.rhy.genes)
# for.rhy.enriched.gos %>% 
#   # only keep the over-enriched terms
#   filter(over_under == "over") %>% 
#   # set FDR at 1%
#   filter(adj_pVal < 0.01) %>% 
#   select(1:6)


# # 2. plotting GO enrichments ----------------------------------------------
# 
# go_enrichment_plot <- function(data, 
#                                function.dir = ".",
#                                category, 
#                                fdr=5, 
#                                clean="yes") {
# 
#   source(paste0(function.dir,"/functions/theme_publication.R"))
# 
#   #Save the data to an object
#   df <- data
# 
#   # Let's read the file with all Cflo GO terms and their categories
#   cflo_gos <- read.csv(paste0(function.dir,"/functions/func_data/distinct_gos_namespace.csv"), header = T, stringsAsFactors = F)
#   cflo_gos <- cflo_gos %>%
#     dplyr::select(1:3) %>%
#     dplyr::select(GO = "GOTerm.identifier", GO_category = "GOTerm.namespace")
# 
#   cflo_gos[cflo_gos$GO_category=="biological_process",]$GO_category <- "BP"
#   cflo_gos[cflo_gos$GO_category=="cellular_component",]$GO_category <- "CP"
#   cflo_gos[cflo_gos$GO_category=="molecular_function",]$GO_category <- "MF"
# 
#   col.scheme <- c("#143740",
#                   "#286E80",
#                   "#3BA4BF",
#                   "#4FDBFF",
#                   "#47C5E6")
# 
# 
#   #Format the dataframe
#   df <- df %>%
#     # only keep the over-enriched terms
#     filter(over_under == "over") %>%
#     # remove the NA term
#     filter(GO_desc != "NA") %>%
#     # set FDR at 1%
#     filter(adj_pVal < (fdr/100)) %>%
#     # add the rhythmicity scores
#     mutate(score = -log(adj_pVal)) %>%
#     # add a column containing the GO categories
#     left_join(cflo_gos, by="GO") %>%
#     # reorder the cols
#     dplyr::select(GO, GO_category, everything())
# 
# 
#   goplot <- ggplot(df) +
#     # set overall appearance of the plot
#     theme_Publication() +
#     # Define the dependent and independent variables
#     aes(x = reorder(GO_desc, score), y = score, fill=GO_category) +
#     # Add a line showing the alpha = 0.01 level
#     geom_hline(yintercept = -log(0.01), size = 1.2, color = "#C7D7D9", alpha=0.8) +
#     # From the defined variables, create a vertical bar chart
#     geom_col(position = "dodge", alpha=0.9, size = 1) +
#     # Set main and axis titles
#     # ggtitle(paste0("GO enrichment | FDR = ",fdr,"%")) +
#     xlab("GOs") +
#     # add caption
#     labs(
#       # title = paste0("GO enrichment | FDR = ",fdr,"%"),
#       # subtitle = sub,
#       caption = "*vertical line denotes FDR of 1%") +
#     ylab(expression(-log[10]*' '*q[enriched])) +
#     # # add annotations
#     # ggrepel::geom_label_repel(aes(label = paste(round((n_GO_DEG/n_GO*100),2), "%", paste(" of ",n_GO,sep=""), sep="")),
#     #                           fill = "transparent",
#     #                           color = 'black',
#     #                           size = 3,
#     #                           direction = "x",
#     #                           ylim=c(10,max(df$score)),
#     #                           point.padding = 0.25,
#     #                           label.padding = 0.25,
#     #                           segment.color = 'transparent',
#     #                           # get rid of the outline for the label
#   #                           label.size = NA) +
#   theme(legend.position = "bottom") +
#     ylim(c(0,max(df$score)+2)) +
#     # Add facetting for each GO category
#     facet_grid(GO_category ~ ., scales = "free_y", space = "free_y") +
#     # Shorten very long labels (GO descriptions)
#     scale_x_discrete(label = function(x) stringr::str_trunc(x, 35)) +
#     # Flip the x and y axes
#     coord_flip() +
#     scale_fill_manual(values= setNames(col.scheme, levels(df$GO_category))) +
#     # scale_fill_manual(values = setNames(c("lightblue", "darkgreen"), levels(tstat$Hemisphere)))
#     theme(strip.background = element_blank(), strip.text = element_blank(), # get rid of facet grid labels
#           plot.title = element_text(hjust = 0.5),
#           axis.line.y = element_line(colour = "transparent",
#                                      size=1),
#           legend.title = element_blank(),
#           # legend.position = "None",
#           plot.caption = element_text(hjust=1),
#           axis.title.y = element_blank()) +
#     guides(
#       fill = guide_legend(
#         title = "Legend Title",
#         override.aes = aes(label = "")))
# 
#   if (clean == "yes") {
#     goplot <- goplot
#   }
# 
#   else(
#     goplot <- goplot +
#       # add annotations
#       ggrepel::geom_label_repel(aes(label = paste(round((n_GO_DEG/n_GO*100),2), "%", paste(" of ",n_GO,sep=""), sep="")),
#                                 fill = "transparent",
#                                 color = 'black',
#                                 size = 3,
#                                 direction = "x",
#                                 ylim=c(10,max(df$score)),
#                                 point.padding = 0.25,
#                                 label.padding = 0.25,
#                                 segment.color = 'transparent',
#                                 # get rid of the outline for the label
#                                 label.size = NA)
#   )
# 
# 
#   return(goplot)
# }


# # 3. Plotting heat maps ------------------------------------------------------
# 
# cflo_heatmap <- function(geneset, 
#                          cluster.r=F, cluster.c=F,
#                          cutree_cols = 1,
#                          title="Heatmap (z-score)", 
#                          show_rownames=F, show_colnames=F,
#                          annotation =T) {
#   
#   # load libraries
#   library(tidyverse)
#   library(pheatmap)
#   library(viridis)
#   
#   # Let's load the datasets:
#   load(file = "./functions/func_data/TC5_core_datasets.RData")
#   load(file = "./functions/func_data/gene_to_annot.RData")
#   
#   oldnames.for <- names(cflo.zscores.for[,-1])
#   oldnames.nur <- names(cflo.zscores.nur[,-1])
#   newnames.for <- c("2F","4F","6F","8F","10F","12F","14F","16F","18F","20F","22F","24F")
#   newnames.nur <- c("2N","4N","6N","8N","10N","12N","14N","16N","18N","20N","22N","24N")
#   
#   genes = as.character(geneset) %>% unique()
#   
#   ## Load the data - zscores
#   cflo.rhy.exp.for <- cflo.zscores.for %>% 
#     #select(gene_name, X2F:X24F) %>% 
#     #mutate(rownames(cflo.annots.exp) = gene_name) %>% 
#     filter(gene_name %in% genes) %>% 
#     rename_at(vars(oldnames.for), ~ newnames.for)
#   cflo.rhy.exp.nur <- cflo.zscores.nur %>% 
#     #select(gene_name, X2F:X24F) %>% 
#     #mutate(rownames(cflo.annots.exp) = gene_name) %>% 
#     filter(gene_name %in% genes) %>% 
#     rename_at(vars(oldnames.nur), ~ newnames.nur)
#   
#   cflo.exp <- cflo.rhy.exp.for %>% 
#     left_join(cflo.rhy.exp.nur, by="gene_name") %>% 
#     na.omit()
#   
#   # Decide if you want gene name or blast annotation for the gene to be displayed
#   if(annotation == T) {
#     
#     cflo.exp.2 <- cflo.exp %>% left_join(gene_to_annot, by="gene_name")
#     cflo.exp.2 <- 
#       cflo.exp.2 %>% 
#       mutate(annot2 <- paste0(gene_name,"| ",annot))
#     
#     rownames(cflo.exp) = cflo.exp.2$annot2
#     cflo.exp <- data.matrix(cflo.exp[-1])
#     
#   }
#   
#   else if (annotation == F) {
#     rownames(cflo.exp) = cflo.exp$gene_name
#     cflo.exp <- data.matrix(cflo.exp[-1])
#   }
#   
#   
#   # I’ll add some column annotations and create the heatmap.
#   # Annotations for:
#   # 1. Is the sample collected during the light or dark phase? 
#   my_sample_col <- data.frame(caste = rep(c("foragers", "nurses"), c(12,12)),
#                               phase = rep(rep(c("light", "dark", "light"), c(5,6,1)),2))
#   
#   row.names(my_sample_col) <- colnames(cflo.exp)
#   # my_sample_col.nur <- data.frame(phase = rep(c("light", "dark", "light"), c(5,6,1)))
#   # row.names(my_sample_col.nur) <- colnames(cflo.rhy.exp.nur)
#   
#   ## Manual color palette
#   my_colour = list(
#     caste = c(foragers = "#F23030", nurses = "#1A80D9"),
#     phase = c(light = "#F2E205", dark = "#010440")
#     #cluster = c("#D97941", "#F2D5C4", "#BF7C63", "#A6401B"),
#     #overlap = c(no = "white", yes = "#D92818")
#   )
#   # my_colour.nur = list(
#   #   phase = c(light = "#FFF640", dark = "grey60")
#   #   #cluster = c("#D97941", "#F2D5C4", "#BF7C63", "#A6401B"),
#   #   #overlap = c(no = "white", yes = "#D92818")
#   #   )
#   
#   # Color scale
#   my.breaks = seq(min(cflo.exp), max(cflo.exp), by=0.06)
#   
#   heat <- pheatmap(cflo.exp, show_rownames = show_rownames, show_colnames = show_colnames,
#                    #annotation_row = my_gene_col.for[,c("cluster","overlap")], 
#                    annotation_col = my_sample_col,
#                    # cutree_rows = 4,
#                    cutree_cols = cutree_cols,
#                    annotation_colors = my_colour,
#                    border_color=FALSE,
#                    cluster_cols = cluster.c,
#                    # user can specify if clustering occurs or not
#                    cluster_rows = cluster.r,
#                    breaks = my.breaks,
#                    ## color scheme borrowed from: 
#                    color = inferno(length(my.breaks) - 1),
#                    ## remove annotation legend all together
#                    annotation_legend = T,
#                    # treeheight_row = 0, 
#                    # treeheight_col = 0,
#                    # remove the color scale or not
#                    legend = T,
#                    main = as.character(title))
#   
#   # nur.rhy.heat <- pheatmap(cflo.rhy.exp.nur, show_rownames = F, show_colnames = F,
#   #                          # annotation_row = my_gene_col.nur[,c("cluster","overlap")], 
#   #                          annotation_col = my_sample_col.nur,
#   #                          # cutree_rows = 8,
#   #                          # cutree_cols = 2,
#   #                          annotation_colors = my_colour.nur,
#   #                          border_color=FALSE,
#   #                          cluster_cols = F,
#   #                          breaks = my.breaks,
#   #                          ## color scheme borrowed from: 
#   #                          color = inferno(length(my.breaks) - 1),
#   #                          ## remove annotation legend all together/or keep it
#   #                          annotation_legend = T,
#   #                          # treeheight_row = 0, 
#   #                          # treeheight_col = 0,
#   #                          # remove the color scale or not
#   #                          legend = T,
#   #                          main = "Nurses (z-scores)")
#   
#   ## Removing specific elements from the figure
#   # Code borrowed from: https://www.biostars.org/p/351551/
#   
#   # library(grid)
#   # grid.ls(grid.force())
#   # # removing the column annotation 
#   # grid.gedit("col_annotation", gp = gpar(col = "white", fill = "white", text = ""))
#   # # removing the "phase" legend title from the annotation legend
#   # grid.gedit("GRID.text.1225", # THIS NUMBER WILL CHANGE EVERY TIME YOU PLOT
#   #            gp = gpar(col="white"))
#   
#   
#   # arrange the two heatmaps in a single figure
#   # library(gridExtra)
#   # for.nur.rhy.heat <- grid.arrange(grobs = list(for.rhy.heat[[4]], nur.rhy.heat[[4]]), ncol=2)
#   # 
#   # # return the heatmaps as a list
#   # return(list(for.nur.rhy.heat, for.rhy.heat[[4]], nur.rhy.heat[[4]]))
#   
#   return(heat)
# }
# 
# 
# 

# 4. Checking sig. overlap ----------------------------------------------

check_overlap <- function(list1, # first list of genesets
                          list2, # second list of genesets that will be compared to list1 in a pair-wise manner
                          tot_genes, # list of all genes or its length to be used as a background
                          text.col = "grey80",   
                          function.dir = ".",
                          fdr=5) {
  
  ## CHECK FOR OVERLAP
  pacman::p_load(GeneOverlap)
  
  # how many genes are there in Cflo genome
  nGenesCflo = as.numeric(ifelse(is.character(tot_genes),length(tot_genes),tot_genes))
  
  ## make a GOM object
  gom <- newGOM(list1, list2, genome.size = nGenesCflo)
  
  gom.plot <- drawHeatmap(gom,
              adj.p=T,
              cutoff=(fdr/100),
              what="odds.ratio",
              # what="Jaccard",
              log.scale = T,
              note.col = text.col)
  
  print(gom.plot)
  
  gom.stats <- list()
  
  # length of genesets in each list
  gom.stats[[1]] <- c(lapply(list1, length) %>% unlist(), lapply(list2, length) %>% unlist())
  
  # stats from running GeneOverlap
  gom.stats[[2]] <- getMatrix(gom, name = "intersection")
  gom.stats[[3]] <- getMatrix(gom, name = "odds.ratio")
  gom.stats[[4]] <- getMatrix(gom, name = "pval")
  names(gom.stats) <- c("length.geneset","intersection","odds.ratio","pval")
  
  return(gom.stats)
  print(gom.stats)
  
}
  

# 5. TC7_annotator --------------------

TC7_annotator <- function(genes, gamma.pval=0.05, FDR=5) {
  
  set.seed(420)
  # load libraries
  pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
  pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)
  # set conflict preference
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
  conflict_prefer("layout", "plotly")
  
  # load database
  db <- dbConnect(RSQLite::SQLite(), "/Users/biplabendudas/Documents/GitHub/R-scripts_zombie_ant_lab/RSQLite/sql_dbs/TC5_data.db")
  # writeLines("Loading 'TC7_ejtk' database that contains all ejtk-output for TC7")
  ejtk.db <- dbConnect(RSQLite::SQLite(),
                       paste0(path_to_repo, "/data/databases/TC7_ejtk.db"))
  # writeLines("Loading 'TC7_data' database that contains all expression data for TC7")
  data.db <- dbConnect(RSQLite::SQLite(),
                       paste0(path_to_repo, "/data/databases/TC7_data.db"))
  
  # Retrieve the gene annotations:
  cflo.annots <-
    tbl(db, "annot_fpkm") %>% 
    select(gene_name, blast_annotation = old_annotation, signalP, TMHMM, GOs, pfams) %>% 
    collect()
  
  ## Specify sample.names
  sample.name <- c("cflo_control","cflo_ophio-infected","cflo_beau-infected")
  
  ## Retrive eJTK data for rhythmicity in controls v. ocflo-inf v. beau-inf
  rhy.24 <- list()
  rhy.12 <- list()
  rhy.08 <- list()
  
  # rhy.genes <- list()
  
  for (i in 1:length(sample.name)) {
    
    ## Load all the rhythmic genesets 
    ## Note, ordered according to their p-value; highly rhythmic at the top.
    
    # Circadian genes (period = 24h)
    ## get the gene-names for sig. 24h-rhythmic genes
    rhy.24[[i]] <-
      tbl(ejtk.db, paste0(sample.name[i],"_zscores_24h")) %>% 
      filter(GammaP < gamma.pval) %>% 
      select(ID, GammaP) %>% collect() %>% arrange(GammaP) %>%
      select(ID) %>% pull()
    
    # Ultradian genes (period = 12h)
    ## get the gene-names for sig. 12h-rhythmic genes
    rhy.12[[i]] <-
      tbl(ejtk.db, paste0(sample.name[i],"_zscores_12h")) %>%
      filter(GammaP < gamma.pval) %>%
      select(ID, GammaP) %>% collect() %>% arrange(GammaP) %>%
      select(ID) %>% pull()
    
    # Ultradian genes (period = 08h)
    ## get the gene-names for sig. 08h-rhythmic genes
    rhy.08[[i]] <-
      tbl(ejtk.db, paste0(sample.name[i],"_zscores_08h")) %>%
      filter(GammaP < gamma.pval) %>%
      select(ID, GammaP) %>% collect() %>% arrange(GammaP) %>%
      select(ID) %>% pull()
    
    # writeLines(paste0("How many sig. rhythmic genes are there in ", sample.name[i],"  heads? (GammaP < ",gamma.pval, ")"))
    # rhy.genes[[i]] <- list(rhy.24[[i]],rhy.12[[i]],rhy.08[[i]])
    # names(rhy.genes[[i]]) <- paste0(sample.name[i], c("_24h", "_12h", "_08h"))
    # print(sapply(rhy.genes[[i]],length))
    
  }
  
  
  # # Retrieve the eJTK data
  # ejtk.24 <- tbl(db, "ejtk_all") %>% collect()
  # ejtk.12 <- tbl(db, "ejtk_12h_all") %>% collect()
  # ejtk.8 <- tbl(db, "ejtk_8h_all") %>% collect()
  # ## make gene lists:
  # # Circadian genes - 24-hour
  # for24 <- ejtk.24 %>% filter(GammaP < 0.05) %>% filter(caste == "for") %>% arrange(GammaP) %>% pull(gene_name)
  # nur24 <- ejtk.24 %>% filter(GammaP < 0.05) %>% filter(caste == "nur") %>% arrange(GammaP) %>% pull(gene_name)
  # # Ultradian genes - 12-hour
  # for12 <- ejtk.12 %>% filter(GammaP < 0.05) %>% filter(caste == "for") %>% arrange(GammaP) %>% pull(gene_name)
  # nur12 <- ejtk.12 %>% filter(GammaP < 0.05) %>% filter(caste == "nur") %>% arrange(GammaP) %>% pull(gene_name)
  # # Ultradian genes - 8-hour
  # for8 <- ejtk.8 %>% filter(GammaP < 0.05) %>% filter(caste == "for") %>% arrange(GammaP) %>% pull(gene_name)
  # nur8 <- ejtk.8 %>% filter(GammaP < 0.05) %>% filter(caste == "nur") %>% arrange(GammaP) %>% pull(gene_name)
  # 
  # # Retrive the DEG (for-nur) data
  # tc5.deg <- tbl(db,"TC5_degs_all") %>% 
  #   filter(significant == "yes") %>% 
  #   collect() %>% 
  #   mutate(FC_DOL = round(2^abs(logFC), 2)) %>% 
  #   rename(DOL_direction = upregulation) %>% 
  #   select(gene_name,FC_DOL,DOL_direction)
  # 
  # # Retrieve DEG (control-biting) data {Will et al. 2020}
  # ophio.deg <-
  #   tbl(db,"ophio_biting_control") %>% 
  #   # only the alive vs control 
  #   filter(sample_1 == "Alive" & sample_2 == "Control") %>% 
  #   # only the significant DEGs
  #   filter(significant == "yes") %>% 
  #   # double check that at least one of the values was ≥ 1 fpkm
  #   filter(value_1 >= 1 | value_2 >= 1) %>%
  #   mutate(logFC = as.numeric(log2.fold_change.)) %>% 
  #   filter(abs(logFC) >= 1) %>% 
  #   # add a column to indicate up or down regulation
  #   mutate(ophio_direction =  ifelse(logFC > 0, "down", "up")) %>% 
  #   select(gene, logFC, ophio_direction) %>% 
  #   collect() %>% 
  #   mutate(FC_ophio = round(2^abs(logFC),2)) %>% 
  #   select(-logFC) %>% 
  #   rename(gene_name = gene)
  # 
  # # Let's get gene-cluster information:
  # ## load zscore datasets
  # cflo.zscores.for <- tbl(db, "zscore_for") %>% collect()
  # cflo.zscores.nur <- tbl(db, "zscore_nur") %>% collect()
  # # Filter the zscores to keep only circadian genes
  # cflo.rhy.exp.for <- cflo.zscores.for %>% 
  #   filter(gene_name %in% for24) %>% as.data.frame()
  # cflo.rhy.exp.nur <- cflo.zscores.nur %>% 
  #   filter(gene_name %in% nur24) %>% as.data.frame()
  # 
  # # Set genes as rownames and convert it into a matrix
  # rownames(cflo.rhy.exp.for) = cflo.rhy.exp.for$gene_name
  # cflo.rhy.exp.for <- as.matrix(cflo.rhy.exp.for[-1])
  # rownames(cflo.rhy.exp.nur) = cflo.rhy.exp.nur$gene_name
  # cflo.rhy.exp.nur <- as.matrix(cflo.rhy.exp.nur[-1])
  # 
  # # Hierarchical clustering of the for24 and nur24 genesets
  # my_hclust_gene.for <- hclust(dist(cflo.rhy.exp.for), method = "complete")
  # my_hclust_gene.nur <- hclust(dist(cflo.rhy.exp.nur), method = "complete")
  # 
  # # Make annotations for the gene-clusters 
  # ## Foragers
  # my_gene_col.for <- cutree(tree = as.dendrogram(my_hclust_gene.for), k = 4)
  # my_gene_col.for <- data.frame(for_cluster = my_gene_col.for)
  # my_gene_col.for$gene_name <- row.names(my_gene_col.for)
  # my_gene_col.for[my_gene_col.for$for_cluster == "1",]$for_cluster <- "night-peaking"
  # my_gene_col.for[my_gene_col.for$for_cluster == "2",]$for_cluster <- "day-peaking"
  # ## Nurses
  # my_gene_col.nur <- cutree(tree = as.dendrogram(my_hclust_gene.nur), k = 4)
  # my_gene_col.nur <- data.frame(nur_cluster = my_gene_col.nur)
  # my_gene_col.nur$gene_name <- row.names(my_gene_col.nur)
  # my_gene_col.nur[my_gene_col.nur$nur_cluster == "1",]$nur_cluster <- "day-peaking"
  # my_gene_col.nur[my_gene_col.nur$nur_cluster == "2",]$nur_cluster <- "night-peaking"
  # my_gene_col.nur[my_gene_col.nur$nur_cluster == "3",]$nur_cluster <- "night-peaking"
  # my_gene_col.nur[my_gene_col.nur$nur_cluster == "4",]$nur_cluster <- "night-peaking"
  # 
  # # # Retrieve the orthology data
  # # load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/reciprocal_blast/TC5_cflo_blastp_summ.RData")
  # 
  
  ## Let's make the dataframe to return
  my_df <- 
    
    # get the annotations for each gene
    cflo.annots %>% 
    # select the genes and their annotation
    filter(gene_name %in% genes) %>% 
    # keep the order of the gene list
    arrange(match(gene_name, genes)) %>% 
    
    # # add differential rhythmicity data
    # mutate(
    #   rhy_for = ifelse(
    #     gene_name %in% for24, "24", ifelse(
    #       gene_name %in% for12, "12", ifelse(
    #         gene_name %in% for8, "8", ""
    #       )
    #     )
    #   ),
    #   rhy_nur = ifelse(
    #     gene_name %in% nur24, "24", ifelse(
    #       gene_name %in% nur12, "12", ifelse(
    #         gene_name %in% nur8, "8", ""
    #       )
    #     )
    #   )
    # ) %>% 
    # 
    # # add cluster information
    # left_join(my_gene_col.for, by="gene_name") %>% 
    # left_join(my_gene_col.nur, by="gene_name") %>% 
    # 
    # # add differential expression data (for v. nur)
    # left_join(tc5.deg, by = "gene_name") %>% 
    # 
    # # add differential expression data (control v. ophio-biting)
    # left_join(ophio.deg, by = "gene_name") %>% 
    
    ## Add rhythmicity data
    mutate(rhy24_control = ifelse(gene_name %in% rhy.24[[1]], "yes", "no"),
           rhy24_ocflo_inf = ifelse(gene_name %in% rhy.24[[2]], "yes", "no"),
           rhy24_beau_inf = ifelse(gene_name %in% rhy.24[[3]], "yes", "no")) %>% 
  
    # order them accordingly
    select(gene_name, blast_annotation, 
           rhy24_control, rhy24_ocflo_inf, rhy24_beau_inf,
           # rhy_for, rhy_nur, 
           # for_cluster, nur_cluster,
           # FC_DOL, DOL_direction,
           # FC_ophio, ophio_direction,
           everything()) %>% 
    
    # remove all the NAs and replace them with empty space
    mutate_each(funs(replace(., which(is.na(.)), ""))) #%>% 
    
    # # add dmel-orthologs
    # left_join(
    #   (cflo.dmel.summ %>% select(gene_name = cflo_gene, dmel_gene)), by = "gene_name"
    # )
  
  
  
  return(my_df)
  
}
