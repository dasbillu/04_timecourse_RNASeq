rm(list = ls())

# Load libraries
pacman::p_load(tidyverse, WaveletComp, ggplot2, conflicted)
pacman::p_load(patchwork)
conflict_prefer("rename", "plyr")
conflict_prefer("filter", "dplyr")

# Specify pathnames
path_to_figures <- paste0(getwd(),"/results/figures/")

# Load functions
## ggplot theme
source(file = "./functions/theme_publication.R")
## summary function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else  length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(as.numeric(datac$N))  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
## saving ggplots
bd.saveplot <- function(plot,name="myplot",width,height) {
  savingto <- paste0(path_to_figures,name,".png")
  writeLines(paste0("Plot is being saved to \n",savingto))
  png(savingto, 
      width = width, height = height, units = "cm", res = 300)
  print(plot)
  trash <- dev.off()
}



# Part 1 ------------------------------------------------------------------

## Format the data ------------------------------------------------------
# let's look at what the data structure looks like
dat <- read.csv("./data/01_activity_data/pogo_C01_activity_tube_21Nov22.csv", 
                header = T, na.strings = c(" ","NA", "", "-"), 
                stringsAsFactors = F)

## Selecting only the columns that I am interested in & summary of data
dat <- 
  dat %>% 
  select(day, phase, zt, tube_1:tube_6) %>% 
  # filter(day %in% c(1:7)) %>% 
  filter(!is.na(day)) %>% 
  mutate(phase=ifelse(phase=="lights off", "dark", "light")) %>% 
  mutate(phase=factor(phase, levels = c("light","dark"))) %>% 
  mutate(day=as.factor(day)) %>% 
  mutate(zt=as.numeric(zt))

# First thing first, how much data do we got?
dat %>% 
  # summarize the data to see how many days of data do we have 
  group_by(day, phase) %>% 
  arrange(day, phase) %>% 
  summarise(datapoints = sum(!is.na(tube_1)))

dat %>% head()
summary(dat) ## We have 84 data points with NAs; all night-time data is missing

# Let's plot the data --------------

# Plot raw ----------------------------------------------------------------
dat.light <- 
  dat %>% 
  filter(phase == "light") %>% 
  pivot_longer(cols = tube_1:tube_6) %>% 
  group_by(day, phase, zt) %>% 
  mutate(total=sum(value, na.rm = T),
         mean=mean(value, na.rm = T),
         median=median(value, na.rm = T),
         max=max(value, na.rm = T)) %>% 
  select(-name, -value) %>% 
  distinct()

foo <-
  dat.light %>% 
  right_join(dat %>% select(day, phase, zt)) %>% 
  group_by(day,phase,zt) %>% 
  mutate(total = ifelse(is.na(total),0,total)) %>% 
  mutate(mean = ifelse(is.na(mean),0,mean)) %>% 
  mutate(median = ifelse(is.na(median),0,median)) %>% 
  mutate(max = ifelse(is.na(max),0,max)) %>% 
  mutate(time2 = (as.numeric(day)-1)*24+as.numeric(zt))

# bd.saveplot(
#   name = "activity_daybyday_total",
#   width = 35, height = 10,
  foo %>% 
    ggplot(aes(x=time2,y=total)) +
    # ggplot(aes(x=1:nrow(foo),y=total)) +
    geom_line(size=2) +
    theme_Publication() +
    scale_x_continuous(breaks = seq(0,max(foo$time2), 12)) +
    labs(
      x="time (hour)",
      y="total activity",
      caption = "12:12h light-dark cycles \nLight was on at 0h and off at 12h, \nand on again at 24h."
    )
# )

# zscore data ---------------------------------------------------------------
bar <- 
  dat.light %>% 
  ungroup() %>% 
  group_by(day) %>% 
  mutate(ztotal=(total-mean(total))/sd(total)) %>% 
  mutate(zmedian=(median-mean(median))/sd(median)) 
  

# Subset ------------------------------------------------------------------
  ## filter to keep data for days 1 through 7
bar <- bar |> filter(day %in% c(1:7))
  

# Plot average foraging ---------------------------------------------------

## we will calculate the mean (±SE) activity of the ants at each time point
pd <- position_dodge(0.1)

## only for days 1 through 7
p1 <- 
  bar %>% 
  # make the summary table
  summarySE(.,
            # specify your measurevar
            measurevar= "ztotal", 
            groupvars=c("zt","phase")) %>%
  
  # make the value column
  mutate(value=ztotal) %>% 
  
  # Plot
  ggplot(aes(x= as.numeric(as.character(zt)), y=value)) +
  theme_Publication() +
  scale_colour_Publication() +
  xlab("") +
  ylab("") +
  ggtitle("") +
  theme(
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    legend.position='none',
    legend.text = element_blank()
  ) +
  theme(plot.title = element_blank()) +
  theme(panel.grid.major.x = element_line(colour = "#808080", size=0.1),
        panel.grid.major.y = element_line(colour = "#808080", size=0.2)) +
  
  ## plot the line connecting the dots
  geom_line(position=pd,
            col="#F2CB05", size=2, alpha=1) +
  ## Add error bar here
  geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                width=.2, position=pd, col="black", alpha = 0.7) +
  # Add the points on top of the error bars
  geom_point(position=pd, size=3.5,
             col="black", fill="black",
             show.legend = F, color="black", pch=21, alpha=0.9) +
  
  theme(text = element_text(size = 25, colour = 'black'),
        legend.position = "none") +
  # forcing the y-axis to start at 0
  expand_limits(x = 0, y = 0) +
  # setting breaks in the x-axis
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  
  # title/caption
  labs(y = "Activity (mean ± SE)",
       x = "Zeitgeber Time (h)"
       # subtitle="Total number of ants",
       # caption = "\n Total number of ants observed in the tube,\n at each time point,\n over a 3 min window"
       )

## Save the plot
bd.saveplot(
  name = "total_ants_summarized",
  width = 20, height = 20,
  p1
)



# Wavelet analyses --------------------------------------------

set.seed(420)

# use data only for daytime
w <- analyze.wavelet(
  bar, 
  ## Analyze which data: Total/Max/Median/Mean
  "ztotal",
  loess.span = 0,
  method = "white.noise",
  dt = 1, dj = 1/100,
  lowerPeriod = 2, upperPeriod = 14,
  make.pval = TRUE, n.sim = 100
)

### Wavelet power (dominant periods with accurate timeline)
# bd.saveplot(
#   name = "wavelet_decomposition_power",
#   width = 15, height = 15,
wt.image(
  w, 
  color.key = "i", 
  n.levels = 6,
  legend.params = list(
    n.ticks = 7
  ),
  siglvl = 0.05,
  main = "Total activity",
  spec.time.axis = list(
    at = c(1:7)*12,
    labels = paste0(
      "0-11h; D",
      c(1:7)
    ),
    las = 1
  ),
  spec.period.axis = list(
    at = c(3,6,12), 
    las=1
  ),
)
# )
### Average wavelet power
# bd.saveplot(
#   name = "wavelet_decomposition_power_average",
#   width = 8, height = 15,
wt.avg(
  w, 
  # maximum.level = maximum.level, 
  spec.period.axis = list(at = c(3,6,12), las=1),
  periodtck = 1, periodtcl = NULL,
  lwd = 4, 
  sigcex = c(2,1.5), exponent = 1,
  siglvl = c(0.05),
  sigcol = c("red"),
  show.legend = F,
  spec.avg.axis	= list(
    at = c(0.1, 0.2, 0.3)
  ),
  # legend.coords = "bottomleft",
  main = "Total activity"
) 
# )



## Reconstruct using all relevant periods
reconstruct(w, plot.waves = F, lwd = c(1,2),
            legend.coords = "none"
            # spec.time.axis = list(at = seq(1,45, by=3),
            #                       labels = seq(3,47, by=3))
            )
## Reconstruct using periods 5-7
reconstruct(w, sel.period = c(5:7), plot.waves = F, lwd = c(1,2),
            legend.coords = "none",
            spec.time.axis = list(at = seq(1,45, by=3),
                                  labels = seq(3,47, by=3)))


# END ---------------------------------------------------------


