
# Table of contents -------------------------------------------------------
### Function 1; 
# name: zplot(gene_names, set); 
# desc: plot zscores 
### Function 2;
# name: multi.plot(plotlist, rows, cols, multipages); 
# desc: plot multiple ggplots into one figure

# LOAD PACKAGES
# load the theme_Publication() 
path_to_repo = "~/Documents/05_postdoc/Deborah_Gordon/04_research/2022/04_timecourse_RNASeq"
source(paste0(path_to_repo,"/functions/theme_publication.R"))
source(paste0(path_to_repo,"/functions/summarySE.R"))
pacman::p_load(RSQLite, dbplyr)
pacman::p_load(gridExtra, tidyverse, ggplot2, conflicted)

# LOAD DATA
## set path to your working directory
path_to_repo = "~/Documents/05_postdoc/Deborah_Gordon/04_research/2022/04_timecourse_RNASeq"
## specify sample name
sample.name <- c("pbar_ld","pbar_dd")
## load databases
## 1. TC9_ejtk.db
ejtk.db <- dbConnect(RSQLite::SQLite(),
                     paste0(path_to_repo, "/data/databases/TC9_ejtk.db"))
## 2. TC9_data.db
data.db <- dbConnect(RSQLite::SQLite(),
                     paste0(path_to_repo, "/data/databases/TC9_data.db"))

## load zscore dataset
zscore.dat <- list()
for (i in 1:length(sample.name)){
  zscore.dat[[i]] <- data.db %>% tbl(., paste0(sample.name[i],"_zscores")) %>% collect()
}
## load FPKM dataset
fpkm.dat <- list()
for (i in 1:length(sample.name)){
  fpkm.dat[[i]] <- data.db %>% tbl(., paste0(sample.name[i],"_fpkm")) %>% collect()
}

# Function 1A: zplot() -----------------------------------
## name: zplot(gene_names, set); 
## desc: plot zscores 
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
## 
## gene_names: a character vector; list of names of genes
## set: a character; possible values include 
##        set = "ld" (plot only ld gene expression),
##        set = "dd" (plot only dd gene expression), or
##        set = "both" (plot both, ld and dd, gene expression)
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-

zplot <- function(gene_names, set = "both", lwd=1.5, alpha=0.9) {
  
  # save the list of gene names in a character vector
  g <- gene_names
  
  # Read the zscores for ld and dd datasets
  ld_zscores <- zscore.dat[[1]]
  dd_zscores <- zscore.dat[[2]]
  
  col.scheme <- c("#D91A2A", # red - LD
                  "#2182BF") # blue - DD
  
  # Transforming the data to be able to use ggplot
  # Let's try the logic on a dummy subset
  dummy.ld <- ld_zscores %>%
    filter(gene_name %in% g) %>% 
    gather(ZT, zscore, -1) %>% 
    # arrange(gene_name) %>% 
    arrange(match(gene_name, g)) %>% 
    mutate(set = "ld") %>% 
    mutate(ZT = readr::parse_number(ZT))
  
  dummy.dd <- dd_zscores %>%
    filter(gene_name %in% g) %>% 
    gather(ZT, zscore, -1) %>% 
    # arrange(gene_name) %>% 
    arrange(match(gene_name, g)) %>% 
    mutate(set = "dd") %>% 
    mutate(ZT = readr::parse_number(ZT))
  
  # make an if else statement to modify the dataset as per the "set" parameter
  if (set == "ld" | set == "ld" | set == "ld" | set == "f"){
    dummy <- dummy.ld
    col.scheme <- col.scheme[1]
  } else if (set == "dd" | set == "dd" | set == "dd" | set == "n") {
    dummy <- dummy.dd
    col.scheme <- col.scheme[2]
  } else if (set == "all" | set == "both" | set == "a" | set == "b") {
    dummy <- rbind(dummy.ld, dummy.dd)
    col.scheme <- col.scheme
  } else {
    print("Invalid value for set. Use one of the following options: for, dd, all.")
    # stop the function and give the user an error.
    stop();
  }
  
  # make the gene_name and set columns as factors
  dummy[[1]] <- as.factor(dummy[[1]])
  dummy[[4]] <- as.factor(dummy[[4]])
  # dummy[[5]] <- as.factor(dummy[[5]])
  
  # Initialize a list to save the plots
  l <- list()
  
  # Let's plot
  pd <- position_dodge(0.2)
  l <- lapply(unique(dummy[[1]]), function(i) {
    ggplot(dummy[dummy$gene_name==i,], 
           aes(x=as.numeric(as.character(ZT)), y=zscore)) + 
      #geom_errorbar(aes(ymin=Corrected_Exp-se, ymax=Corrected_Exp+se), width=.1, position=pd) +
      geom_hline(yintercept=0, color = "grey60", size = 1, lty=1, alpha=0.7) +
      ## if you need highlighting parts of the graph (dark phase in my case)
      geom_rect(aes(xmin = 11.5, xmax = 23.5, ymin = -Inf, ymax = Inf),    
                fill = "lightgrey", alpha = 0.02, color=NA) +
      #facet_grid(~gene_name, scales = "free") +
      facet_wrap(~ gene_name) +
      #ggtitle("TC5 - Z-scores") +
      xlab("") +
      ylab("expression (z-scores)") +
      theme_Publication() +
      scale_x_continuous(breaks = seq(1,23,2)) +
      # scale_y_continuous(limits = c(min(dummy[[3]]),max(dummy[[3]])),n.breaks = 4) +  
      scale_y_continuous(limits = c(-3.5, 3.5), breaks = c(-3,0,3)) +  
      geom_line(data = dummy[dummy$gene_name==i,], position=pd, 
                aes(col=as.factor(set)), size=lwd, alpha=alpha) +
      scale_color_manual(values=col.scheme) + 
      theme(text = element_text(size = 20, colour = 'black'),
            legend.position = "none",
            axis.title.x=element_blank(),
            axis.text.x=element_blank()) +
      # set transparency
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA))
  })
  
  # Let's name the plots with their resp. gene names
  names(l) <- g;
  
  # return the list with all the plots and the data frame containing   
  return(l);
  
}



# Compiling the multiple plots into one figure
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
## For plotting the multiple gene plots into one or multiple pages, 
## use the multi.plot() function 
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-



# Function 2: multi.plot() -----------------------------------

## name: multi.plot(plotlist, rows, cols, multipages); 
## desc: plot multiple ggplots into one figure 
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
## 
## plotlist = list of ggplots
## rows = number of plots in each column of the page (deafult is 3)
## cols = number of plots in each row of the page (default is 4)
## multipages = binary; 0 indicates print all plots in one page, 1 indicates multi-page plot
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-

multi.plot <- function(plotlist, rows=3, cols=4, multipages="0") {
  # load required packages
  library(gridExtra)
  # coerce the value for multipages into character 
  multi <- as.character(multipages)
  # plot them in 1 page
  if(multipages == "0") {
    # number of total plots?
    tplots <- length(plotlist);
    # number of possible plots given rows and cols?
    pplots <- rows*cols;
    if(tplots <= pplots){
      plots <- do.call(grid.arrange, c(plotlist, ncol=cols, nrow=rows))
    } else {
      print("Can't fit all plots in one page. Change rows or/and cols.")
      stop()
    }
    
  } else if(multipages == "1"){
    plots <- marrangeGrob(plotlist, nrow=rows, ncol=cols)
  } else {
    print("Invalid value for multipages./n 0=one page; 1=multiple pages")
    stop()
  }
  # return the plots object which is of class "arrangelist" "list"
  return(plots)
}

# multi.plot() returns one object of the class "arrangelist" "list"

### Usage for multi.plot() function
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
##
## Call the multi.plot function
## a <- multi.plot(plotlist = rando.plots$plots, rows=2, cols=2, multipages = 1)
## a
## Save the file
## ggsave("name_of_file.png", a, bg = "transparent") # didn't run yet, need to check
## If multipages = 1, you can save the multiple page file into a pdf 
## ggsave("multipage.pdf", a)
####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-


# Function 3: exp.plot() -----------------------------------

####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-####-
exp.plot <- function(gene_names, set = "both", log=T, lwd=2, alpha=0.9) {
  
  # save the list of gene names in a character vector
  g <- gene_names
  
  # Read the FPKM for set
  ld <- fpkm.dat[[1]]
  dd <- fpkm.dat[[2]]
  
  # Transforming the data to be able to use ggplot
  dummy.ld <- ld %>%
    filter(gene_name %in% g) %>% 
    gather(ZT, exp, -1) %>% 
    # mutate(log.exp = log10(exp+1)) %>% 
    # dplyr::select(-exp) %>% 
    arrange(gene_name) %>% 
    mutate(set = "ld") %>% 
    mutate(ZT = readr::parse_number(ZT))
  
  dummy.dd <- dd %>%
    filter(gene_name %in% g) %>% 
    gather(ZT, exp, -1) %>% 
    # mutate(log.exp = log10(exp+1)) %>% 
    # dplyr::select(-exp) %>% 
    arrange(gene_name) %>% 
    mutate(set = "dd") %>% 
    mutate(ZT = readr::parse_number(ZT))
  
  col.scheme <- c("#D91A2A", # red - LD
                  "#2182BF") # blue - DD
  
  # make an if else statement to modify the dataset as per the "set" parameter
  if (set == "ld" | set == "ld" | set == "ld" | set == "f"){
    dummy <- dummy.ld
    col.scheme <- col.scheme[1]
  } else if (set == "dd" | set == "dd" | set == "dd" | set == "n") {
    dummy <- dummy.dd
    col.scheme <- col.scheme[2]
  } else if (set == "all" | set == "both" | set == "a" | set == "b") {
    dummy <- rbind(dummy.ld, dummy.dd)
    col.scheme <- col.scheme
  } else {
    print("Invalid value for set. Use one of the following options: for, dd, all.")
    # stop the function and give the user an error.
    stop();
  }
  
  # make the gene_name and set columns as factors
  dummy[[1]] <- as.factor(dummy[[1]])
  dummy[[4]] <- factor(dummy[[4]],
                         levels = unique(dummy[[4]]))
  
  # Initialize a list to save the plots
  l <- list()
  # Let's plot
  library(ggplot2)
  pd <- position_dodge(0.2)
  
  if(log==T) {
  
  l <- lapply(sort(unique(dummy[[1]])), function(i) {
    ggplot(dummy[dummy$gene_name==i,], 
           aes(x=as.numeric(as.character(ZT)), y=log2(exp+1))) + 
      ## if you need highlighting parts of the graph (dark phase in my case)
      geom_rect(aes(xmin = 11.5, xmax = 23.5, ymin = -Inf, ymax = Inf),    
                fill = "lightgrey", alpha = 0.02, color=NA) +
      #facet_grid(~gene_name, scales = "free") +
      facet_wrap(~ gene_name, scales = "free_y") +
      #ggtitle("TC5 - Z-scores") +
      xlab("") +
      ylab("") +
      theme_Publication() +
      scale_x_continuous(breaks = seq(1,23,2)) +
      scale_y_continuous(n.breaks = 3) +
      geom_line(data = dummy[dummy$gene_name==i,], position=pd, 
                aes(col=as.factor(set)), size=lwd, alpha=alpha) +
      scale_color_manual(values=col.scheme) + 
      theme(text = element_text(size = 15, colour = 'black'),
            legend.position = "none",
            axis.title.x=element_blank(),
            axis.title.y = element_blank(),
            axis.text.x=element_blank()) +
      # set transparency
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA))
  })
  }
  
  else if (log==F) {
    l <- lapply(sort(unique(dummy[[1]])), function(i) {
      ggplot(dummy[dummy$gene_name==i,], 
             aes(x=as.numeric(as.character(ZT)), y=exp)) + 
        ## if you need highlighting parts of the graph (dark phase in my case)
        geom_rect(aes(xmin = 11.5, xmax = 23.5, ymin = -Inf, ymax = Inf),    
                  fill = "lightgrey", alpha = 0.02, color=NA) +
        #facet_grid(~gene_name, scales = "free") +
        facet_wrap(~ gene_name, scales = "free_y") +
        #ggtitle("TC5 - Z-scores") +
        xlab("") +
        ylab("") +
        theme_Publication() +
        scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
        scale_y_continuous(n.breaks = 3) +
        geom_line(data = dummy[dummy$gene_name==i,], position=pd, 
                  aes(col=as.factor(set)), size=lwd, alpha=alpha) +
        # geom_point(position=pd, 
        #            aes(fill=as.factor(set), shape = as.factor(set)), size=5, show.legend = F,
        #            color="black", pch=21) +
        #scale_fill_manual(values = c("#F2CB05","#0FBF67")) +
        scale_color_manual(values=col.scheme) + 
        theme(text = element_text(size = 15, colour = 'black'),
              legend.position = "none",
              axis.title.y = element_blank(),
              axis.title.x=element_blank(),
              axis.text.x=element_blank()) +
        # set transparency
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA))
    })
  }
  
  # return the list with all the plots and the data frame containing   
  return(l);
  
}


# Function 4: stacked.zplot() -----------------------------------

# Function name: stacked.zplot
# USage: to plot multiple diel gene expression (as z-scores) stacked on top of each other;
# Application: could be used to report similar diel expression patterns for multiple genes

stacked.zplot <- 
  function(gene_names, set = "all", lwd=1.5, alpha=0.75, bg.alpha = 0.1, text_size=20) {
    
    # save the list of gene names in a character vector
    g <- gene_names
    
    ld_zscores <- zscore.dat[[1]]
    dd_zscores <- zscore.dat[[2]]
    
    # Transforming the data to be able to use ggplot
    # Let's try the logic on a dummy subset
    dummy.ld <-
      ld_zscores %>%
      filter(gene_name %in% g) %>%
      gather(ZT, zscore, -1) %>%
      arrange(match(gene_name, g)) %>% 
      mutate(set = sample.name[1]) %>% 
      mutate(ZT = readr::parse_number(ZT))
    
    dummy.dd <- dd_zscores %>%
      filter(gene_name %in% g) %>% 
      gather(ZT, zscore, -1) %>% 
      arrange(match(gene_name, g)) %>% 
      mutate(set = sample.name[2]) %>% 
      mutate(ZT = readr::parse_number(ZT))
    
    
    col.scheme <- c("#D91A2A", # red - LD
                    "#2182BF") # blue - DD
    
    # make an if else statement to modify the dataset as per the "set" parameter
    if (set == "ld" | set == sample.name[1]){
      dummy <- dummy.ld
      
    } else if (set == "dd" | set == sample.name[2]) {
      dummy <- dummy.dd
    
    } else if (set == "all" | set == "both") {
      dummy <- rbind(dummy.ld, dummy.dd)
      
    } else {
      print("Invalid value for set. Use one of the following options: ld, dd, both")
      # stop the function and give the user an error.
      stop();
    }
    
    # make the gene_name and set columns as factors
    dummy[[1]] <- as.factor(dummy[[1]]) # gene name
    dummy[[4]] <- factor(dummy[[4]],
                            levels=c(sample.name[1],
                                     sample.name[2])) # set
    
    dummy.summary <-
      dummy %>% 
      group_by(set, ZT) %>% 
      summarise(mean_zscore = mean(zscore, na.rm = T))
    
    # Initialize a list to save the plots
    l <- list()
    
    # Let's plot
    pd <- position_dodge(0.2)
    
    l <- lapply(sort(unique(dummy[[4]])), function(i) {
      
      # Define the dataset
      ggplot() + 
        
        # indicate light-dark phase
        geom_rect(aes(xmin = 11.5, xmax = 23.5, ymin = -Inf, ymax = Inf),
                  fill = "lightgrey", alpha = 0.3, color=NA) +
        
        # Plot the individual gene expressions (zscores)
        geom_line(data=dummy[dummy$set==i,], 
                  aes(x=as.numeric(as.character(ZT)), y=zscore, group = gene_name),
                  size=0.5, alpha=bg.alpha, 
                  col = ifelse(i==sample.name[1], c(col.scheme[[1]]),
                               ifelse(i==sample.name[2], c(col.scheme[[2]]),
                                             c(col.scheme[[1]], col.scheme[[2]])))) +
        
        # Plot the mean gene expression for the given gene list
        geom_line(data = dummy.summary[dummy.summary$set == i,],
                  aes(x=as.numeric(as.character(ZT)), y=mean_zscore), 
                  col = "#3C3C40",
                  size = lwd, alpha = alpha) +
        
        
        # Set the theme
        xlab("") +
        ylab("zscore") +
        
        # Set titles for individual plots
        ggtitle(paste0(i)) +
        
        theme_Publication() +
        scale_x_continuous(breaks = c(1,13,23)) +
        # scale_y_continuous(limits = c(min(dummy[[3]]),max(dummy[[3]])), n.breaks = 3) + 
        scale_y_continuous(breaks = c(-3,0,3), limits = c(-3.5,3.5)) + 
        theme(text = element_text(size = text_size, colour = 'black'),
              axis.title.x=element_blank(),
              # axis.text.x=element_blank(),
              legend.position = "none") +
          # set transparency
        theme( 
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA)) 
    })
    
    # return the list with all the plots and the data frame containing   
    return(l);
    
}




#  Function 5: amplitude.plot() ---------------------------------
# Function name: 
# USage: 
# Application: 

amplitude.plot <- 
  function(gene_names, set = "all", lwd=0.75, tc5=F, mean=T, stats=F, ci=T, ...) {
    
    ## set conflict preference
    conflict_prefer("select", "dplyr")
    conflict_prefer("filter", "dplyr")
    conflict_prefer("sd", "stats")
    

    ## PREP THE DATASETS ##
    
    ## save the list of gene names in a character vector
    g <- gene_names
    
    
    # # Read the annotation file for description of the genes
    # all_genes <- read.csv(paste0(path_to_repo,"/functions/func_data/cflo_annots.csv"), header = T, stringsAsFactors = F)
    # all_genes_annots <- all_genes %>% 
    #   dplyr::select(gene_name, annot = old_annotation)
    
    # Transforming the data to be able to use ggplot
    # Let's try the logic on a dummy subset
    dummy.ld <-
      zscore.dat[[1]] %>% 
      filter(gene_name %in% g) %>% 
      gather(ZT, zscore, -1) %>% 
      arrange(match(gene_name, g)) %>% 
      mutate(set = sample.name[1]) %>% 
      mutate(ZT = readr::parse_number(ZT)) %>% 
      na.omit()
    
    dummy.dd <- 
      zscore.dat[[2]] %>%
      filter(gene_name %in% g) %>% 
      gather(ZT, zscore, -1) %>% 
      arrange(match(gene_name, g)) %>% 
      mutate(set = sample.name[2]) %>% 
      mutate(ZT = readr::parse_number(ZT)) %>% 
      na.omit()
    
    
    ## Specify color scheme
    col.scheme <- c("#AD212F","#5A829F")
    
    ## Specify different samples
    sets <- sample.name
    
    # make an if else statement to modify the dataset as per the "set" parameter
    if (set == "ld" | set == sample.name[1]){
      dummy <- dummy.ld
      col.scheme.sub <- col.scheme[1]
      sets.sub <- sets[1]
      
    } else if (set == "dd" | set == sample.name[2]){
      dummy <- dummy.dd
      col.scheme.sub <- col.scheme[2]
      sets.sub <- sets[2]
      
    } else if (set == "all" | set == "both") {
      dummy <- rbind(dummy.ld, dummy.dd)
      col.scheme.sub <- col.scheme[1:2]
      sets.sub <- sets[1:2]
      
    } else {
      print("Invalid value for set. Use one of the following options: ld, dd, both.")
      # stop the function and give the user an error.
      stop();
    }
    
    if (tc5==T){
      dummy <- rbind(dummy.ld, dummy.dd, dummy.cflo.brains)
      col.scheme.sub <- col.scheme[1:3]
      sets.sub <- c(sample.name, "cflo_ld")
    }
    
    # make the gene_name and set columns as factors
    dummy[[1]] <- as.factor(dummy[[1]]) # gene name
    dummy[[4]] <- factor(dummy[[4]], levels = sets.sub) # set
    
    # dummy[[5]] <- as.factor(dummy[[5]]) # annotation column
    
   ### CALCULATE AMPLITUDE
    plot.dat <-
      dummy %>% 
      select(gene_name,set,zscore) %>%
      group_by(gene_name,set) %>% 
      na.omit() %>% 
      mutate(amp=max(zscore)-min(zscore)) %>% 
      mutate(cv=sd(zscore)/mean(zscore)) %>% 
      select(set, amp) %>% 
      distinct()
    
    if (stats==T){
      ## Perform a Kruskal 
      kruskal.res <- ggpubr::compare_means(amp~set, data = na.omit(plot.dat), paired = T, method = "kruskal.test")
      if (kruskal.res$p.adj<0.05) {
        writeLines("\nKruskal wallis test found significant differences in means (q<0.05)\n")
        
        #calculate the summary stats
        writeLines("\nSummary stats:\n")
        plot.dat %>% 
          # make the summary table
          summarySE(.,
                    # specify your measurevar (prop_time)
                    measurevar= "amp", 
                    groupvars=c("set")) %>% 
          select(set, mean_amp=amp,sd:ci) %>% 
          print()
        
        # perform pairwise tests
        writeLines("\nPerforming pairwise Wilcox test...")
        ggpubr::compare_means(amp~set, data = plot.dat, paired = F, 
                                              method = "wilcox.test") %>% 
          select(-p) %>% 
          # select(-p.ldmat) %>% 
          select(-p.signif) %>% 
          print()
      } else {
        writeLines("\nKruskal wallis test did not find any evidence for difference in means (q > 0.05)\n")
      }
    } 
    

    plot <- 
      plot.dat %>% 
      ggplot() +
      geom_boxplot(aes(x=set,y=amp, col=set), size=lwd, outlier.size = 0.3) +
      # facet_grid(~set) +
      theme_Publication(base_size = 15) +
      scale_color_manual(values = col.scheme.sub) +
      scale_fill_manual(values = col.scheme.sub) +
      theme(axis.text.x = element_text(color = NA))
    
    ## MEAN ± SD
    
    ## plots with mean ± 95%CI
    pd <- position_dodge(0.1)
    plot.2 <-
      plot.dat %>%
        
        # make the summary table
        summarySE(.,
                  # specify your measurevar (prop_time)
                  measurevar= "amp", 
                  groupvars=c("set")) %>% 
        
        # make the value column
        mutate(value=amp) %>% 
        
        ggplot(aes(y=amp, x=set, col=set)) +
        
        # Add error bar here
        geom_errorbar(aes(ymin=value-sd, ymax=value+sd),
                      width=.4, position=pd, col="black", alpha = 0.7) +
        
        # Add the points on top of the error bars
        geom_point(position=pd, size=7,
                   aes(fill=set),
                   # col="black", fill="black",
                   show.legend = F, pch=21, alpha=1) +
        
        theme_Publication(base_size = 15) +
        scale_color_manual(values = col.scheme.sub) +
        scale_fill_manual(values = col.scheme.sub) +
        
        # scale_y_continuous(breaks = c(0,.2,.4,.6,.8,1),
        #                    labels = c("0", "", "0.4","","0.8",""),
        #                    limits = c(0,1)) +
        
        # ylab("proportion returned at end of Day 1")
        # ylab("cumulative prop. returned at end of Day 2")
        ylab("coefficient of variation\n (mean ± SD)") +
        # coord_flip() 
        theme(axis.text.x = element_text(color = NA))
    
    
    ## MEAN ± 95% CI
    plot.3 <-
      plot.dat %>%
      
      # make the summary table
      summarySE(.,
                # specify your measurevar (prop_time)
                measurevar= "amp", 
                groupvars=c("set")) %>% 
      
      # make the value column
      mutate(value=amp) %>% 
      
      ggplot(aes(y=amp, x=set, col=set)) +
      
      # Add error bar here
      geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                    width=.4, position=pd, col="black", alpha = 0.7) +
      
      # Add the points on top of the error bars
      geom_point(position=pd, size=7,
                 aes(fill=set),
                 # col="black", fill="black",
                 show.legend = F, pch=21, alpha=1) +
      
      theme_Publication(base_size = 15) +
      scale_color_manual(values = col.scheme.sub) +
      scale_fill_manual(values = col.scheme.sub) +
      
      # scale_y_continuous(breaks = c(0,.2,.4,.6,.8,1),
      #                    labels = c("0", "", "0.4","","0.8",""),
      #                    limits = c(0,1)) +
      
      # ylab("proportion returned at end of Day 1")
      # ylab("cumulative prop. returned at end of Day 2")
      ylab("coefficient of variation\n (mean ± 95% CI)") +
      # coord_flip() 
      theme(axis.text.x = element_text(color = NA))
        
        
        
    if (mean==T & ci==F){
      return(plot.2)  
    } else if (mean==F & ci==F) {
      return(plot)
    } else {
      return(plot.3)
    }
    
    
    
    
  }
