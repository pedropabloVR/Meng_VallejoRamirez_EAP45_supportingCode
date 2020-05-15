# Supplementary code accompanying the manuscript: "Distinct domain requirements for EAP45 in HIV budding, late endosomal recruitment, and cytokinesis"

# Pedro Vallejo Ramirez
# Laser Analytics Group, University of Cambridge
# Created: 2019-05-10
# Updated: 2020-04-07

# This code imports the results from the nearest-neighbour analysis performed in 
# matlab (.csv files), and plots the percentage of associated Gag particles, i.e.
# those with a another protein (ALIX, EAP45) within a 1-pixel distance, for:

# - Hela cell lines expressing ALIX, FL EAP45, and the functionally compromised mutants 
#   (deltaG deltaH, delta PTAP, delta PTAP deltaG deltaH).

# - HAP-1EAP45 KO cells with FL EAP45 and delta G delta H EAP45. 

# This script can be used to recreate the plots in following figures:
  
# - Figure 3C
# - Figure 4E
# - Figure 4J
# - Supplementary Figure S2
# This code also includes the statistical analysis accompanying each of these figures.


##----------Import libraries and data files ----------##
{
library(tidyverse)    # Main library for R functions
library(ggplot2)      # For making nice plots
library(reshape2)     # For reshaping data frames 
library(ggbeeswarm)   # for beeswarm dot plots (symmetrically jittered)
library(ggquiver)     # for quiver plots
library(ggpubr)       # for arranging more than one plot in one page
library(RColorBrewer) # Good color palettes for plotting
#require(cowplot)      # for arranging more than one plot in one page
library(dabestr)      # for plotting the differences between distributions as an estimation graphic
#library(multcomp)     # for running Dunnett's test
#library(DescTools)    # for running Dunnett's test
}

# Set random number seed to get reproducibe plots
set.seed(2019)      

# For a full page a4 figure, w = 160 mm, h = ?
# For a half-page a4 figure, w = 80 mm, h = 70 mm (no legend)

# Use custom plotting function
prism_multi <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = "white"
  color.grid.major = "white"
  color.axis.text = "black"
  color.axis.title = "black"
  color.title = palette[9]
  font_size <- 10
  title_size <- 10
  
  # Begin construction of chart
  theme_bw(base_size=11) +
    
    # Format the legend, uncomment position if need be
    theme(legend.text  = element_text(size=font_size, color=color.axis.title))+
    theme(legend.title = element_text(size=font_size, color=color.axis.title))+
    theme(legend.position="bottom")+

    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=font_size, vjust=1.25, hjust = 0.4, face = "bold")) +
    theme(axis.text.x=element_text(size=font_size,color=color.axis.text, margin = unit(c(0.1, 0.5, 0.3, 0.5), "cm"))) +
    theme(axis.text.y=element_text(size=font_size,color=color.axis.text, margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm"))) +
    theme(axis.title.x=element_text(size=font_size,color=color.axis.title, 
                                    vjust=1.25, margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))) +
    theme(axis.title.y=element_text(size=font_size,color=color.axis.title, 
                                    vjust=-4, hjust = 0.5, margin = unit(c(0.2, 0.2, 0.2, 
                                                                                          0.2), "cm"))) 
  #theme(axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0.10)) 
} # plotting function

# Import .csv files with results
{
# HeLa cell data   
# full length (wildtype) EAP45 
data_wt_exp6 <-read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp6_slide3well7_20190526_weka_wtEAP45.csv')
data_wt_exp7 <- read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp7_slide1well6_20190526_weka_wtEAP45.csv')
data_wt_exp8 <- read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp8_slide3well6_20190528_weka_wtEAP45.csv')

# delta G delta H EAP45 mutant
data_mutant_exp7_slide1well7 <- read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp7_slide1well7_20190526_weka_mtEAP45.csv')
data_mutant_exp7_slide2well7 <- read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp7_slide2well7_20190526_weka_mtEAP45.csv')
data_mutant_exp8_slide3well7 <- read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp8_slide3well7_20190715_weka_mtEAP45.csv')

#double mutant (delta PTAP and delta G delta H EAP45)
data_doublemutant_slide2well7   <-read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp9_slide2well7_20190624_weka_doublemtEAP45.csv')
data_doublemutant_slide2well7_2 <-read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp9_slide2well7_20190803_weka_doublemtEAP45.csv')

# delta PTAP mutant 
data_PTAPmutant_slide2well6 <-read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp9_slide2well6_20190701_PTAP_mtEAP45.csv')

# ALIX data 
data_ALIX         <-read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp9_slide4well3_20190625_weka_ALIX.csv')

# Nanobooster data 
data_nanobooster  <-read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp6_slide1well3_20190517_weka_nanobody.csv')

# HAP1-EAP45KO cell data 
data_KOcells_wtEAP45 <- read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp10_slide1well6_20190816_weka_KOcells_wtEAP45.csv')
data_KOcells_mtEAP45 <- read.csv('/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/co-occurrence results/results_exp10_slide1well7_20190816_weka_KOcells_mtEAP45.csv')

}

# Tidy up data and organize into master data frames 
{
pix = 117 # pixel size 
x_vector <- c("x0","x1","x2","x3","x4","x5") # vector with labels for x axis

# Aggregating HeLa cell data from the different experiments for each condition  (MODIFY THIS IF WANT TO REPLOT INDIVIDUAL DATASETS OR CONDITIONS)
rows_wt       <- nrow(data_wt_exp6)+nrow(data_wt_exp7)+nrow(data_wt_exp8)
rows_mt       <- nrow(data_mutant_exp7_slide1well7)+nrow(data_mutant_exp7_slide2well7)+nrow(data_mutant_exp8_slide3well7)
rows_doublemt <- nrow(data_doublemutant_slide2well7)+nrow(data_doublemutant_slide2well7_2)
rows_PTAP     <- nrow(data_PTAPmutant_slide2well6)
rows_ALIX     <- nrow(data_ALIX)

# Tidy up data and aggregate the associated EAP45 and Gag particles into two master data frames with all conditions tested 

# HeLa cells
df_EAP45_aggregated_allconditions     <- data.frame(type   = c(rep('wt',rows_wt),
                                                               rep('EAP45 mt',rows_mt),
                                                               rep('double mt',rows_doublemt),
                                                               rep('PTAP mt',rows_PTAP),
                                                               rep('ALIX',rows_ALIX)),
                                                    x      = c(rep(x_vector,(rows_wt)/6),rep(x_vector,rows_mt/6),rep(x_vector,rows_doublemt/6),rep(x_vector,rows_PTAP/6),rep(x_vector,rows_ALIX/6)),
                                                    value  = c(data_wt_exp6[,7],data_wt_exp7[,7],data_wt_exp8[,7],
                                                               data_mutant_exp7_slide1well7[,7],data_mutant_exp7_slide2well7[,7],data_mutant_exp8_slide3well7[,7],
                                                               data_doublemutant_slide2well7[,7],data_doublemutant_slide2well7_2[,7],
                                                               data_PTAPmutant_slide2well6[,7],
                                                               data_ALIX[,7]))
  
df_Gag_aggregated_allconditions       <- data.frame(type   = c(rep('wt',rows_wt),
                                                               rep('EAP45 mt',rows_mt),
                                                               rep('double mt',rows_doublemt),
                                                               rep('PTAP mt',rows_PTAP),
                                                               rep('ALIX',rows_ALIX)),
                                                    x      = c(rep(x_vector,(rows_wt)/6),rep(x_vector,rows_mt/6),rep(x_vector,rows_doublemt/6),rep(x_vector,rows_PTAP/6),rep(x_vector,rows_ALIX/6)),
                                                    value  = c(data_wt_exp6[,8],data_wt_exp7[,8],data_wt_exp8[,8],
                                                               data_mutant_exp7_slide1well7[,8],data_mutant_exp7_slide2well7[,8],data_mutant_exp8_slide3well7[,8],
                                                               data_doublemutant_slide2well7[,8],data_doublemutant_slide2well7_2[,8],
                                                               data_PTAPmutant_slide2well6[,8],
                                                               data_ALIX[,8]))
# HAP1-EAP45KO cells
df_KOcells_GagAssociated              <- data.frame(type   = c(rep('wt',nrow(data_KOcells_wtEAP45)),
                                                               rep('mt',nrow(data_KOcells_mtEAP45))),
                                                    x      = c(rep(x_vector,nrow(data_KOcells_wtEAP45)/6),
                                                               rep(x_vector,nrow(data_KOcells_mtEAP45)/6)),
                                                    value  = c(data_KOcells_wtEAP45[,8],data_KOcells_mtEAP45[,8]))


# Choose only the first search radius (x1 = 117 nm)
df_EAP45_percentageAssociated       <- subset(df_EAP45_aggregated_allconditions, x=='x1')
df_Gag_percentageAssociated         <- subset(df_Gag_aggregated_allconditions, x=='x1')
df_KOcells_Gag_percentageAssociated <- subset(df_KOcells_GagAssociated, x == 'x1')
}

##----------Plots for figures ----------##

# Figure 3C: Scatter and boxplot for comparing the associated Gag particles with ALIX and FL EAP45 in Hela cells
{
  
  # Choose the ALIX and FL EAP45 datasets only from the master data set
  df_GagAssociated_ALIX_FLEAP45 <- subset(df_Gag_percentageAssociated, type!='EAP45 mt' & type!='PTAP mt' & type!='double mt')
  
  # This works to plot the mean instead of the median as a measure of the data centre
  ggplot(df_GagAssociated_ALIX_FLEAP45,aes(x = type, y = value, colour = type))+
    geom_quasirandom(alpha = 0.8,dodge.width=.8,size = 0.3)+
    geom_boxplot(aes(fill = type), alpha = 0.3,outlier.shape = NA,fatten = NULL)+
    #scale_x_discrete(labels = c("x0"="0","x1"="117","x2" = "234","x3"="351","x4"="468","x5"="585"))+
    labs(y='% associated Gag spots')+ #where # of associated spots refers to spots which have at least 1 neighbour within R_search
    prism_multi()+
    theme(legend.position = "none")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    scale_fill_manual(values=c("#00BFC4", "#F8766D"))+
    scale_color_manual(values=c("#00BFC4", "#F8766D"))
  #stat_summary(fun.y=mean, geom="point", shape=17, size=4)
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/results association gag-eap45 20190510/wt_ALIX_association_vs_radius_gag_20191028.pdf",width = 55, height = 50, units = "mm")
  
}

# Figure 4F: Scatter and boxplot for comparing the associated Gag particles in HeLa cells with FL EAP45 and its functionally compromised mutants.
{
# First, choose FL EAP45 and its mutants from master data frame 
  df_GagAssociated_EAP45_allMutants      <- subset(df_Gag_percentageAssociated, type!='ALIX')

# Order the different conditions to make it easier to display
  df_GagAssociated_EAP45_allMutants$type <- factor(df_GagAssociated_EAP45_allMutants$type, levels = c("wt","EAP45 mt","PTAP mt","double mt"))
  
# Plot the data as a scatter plot with an overlaid boxplot
ggplot(df_GagAssociated_EAP45_allMutants,aes(x = type, y = value, colour = type))+
  geom_quasirandom(alpha = 0.8,dodge.width=.8,size = 0.3)+
  geom_boxplot(aes(fill = type), alpha = 0.3,outlier.shape = NA,fatten=NULL)+
  #scale_x_discrete(labels = c("x0"="0","x1"="117","x2" = "234","x3"="351","x4"="468","x5"="585"))+ # add this line if more than 1 search radius is used
  labs(y='% Associated Gag spots')+ #where # of associated spots refers to spots which have at least 1 neighbour within R_search
  theme(legend.position = "none")+
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  prism_multi()

ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/results association gag-eap45 20190510/ch2-ch1_association_vs_radius_gag_20191028_mean.pdf",width = 105, height = 90, units = "mm")
}

# Figure 4E: Scatter and boxplot for comparing the associated Gag particles in HAP1-EAP45KO cells with FL EAP45 and a delG delH EAP45 mutant. 
{
  # Order the data such that it shows full length EAP45 first followed by the mutant
  df_KOcells_Gag_percentageAssociated$type <- factor(df_KOcells_Gag_percentageAssociated$type, levels = c("wt","mt"))
  
  # Plot the data as a scatter plot with an overlaid boxplot, with the center line showing the mean
  ggplot(df_KOcells_Gag_percentageAssociated,aes(x = type, y = value, colour = type))+
    geom_quasirandom(alpha = 0.8,dodge.width=.8,size = 0.3)+
    geom_boxplot(aes(fill = type), alpha = 0.3,outlier.shape = NA, fatten = NULL)+
    #scale_x_discrete(labels = c("x0"="0","x1"="117","x2" = "234","x3"="351","x4"="468","x5"="585"))+
    labs(y='% Associated Gag spots')+ #where # of associated spots refers to spots which have at least 1 neighbour within R_search
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    prism_multi()+
    scale_fill_manual(values=c("#F8766D", "#7CAE00"))+
    scale_color_manual(values=c("#F8766D", "#7CAE00"))+ # use same colours as in Figure 4F
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = "none")
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/results association gag-eap45 20190510/ch2-ch1_association_vs_radius_aggregate_20191028_KOcells_green_mean.pdf",width = 105, height = 90, units = "mm")

}

# Supplementary Figure S2
{
  # data frames for nano booster data
  df_nanobooster <- data.frame(type   = c(rep('GFP-NB',nrow(data_nanobooster))),
                               x      = c(rep(x_vector,nrow(data_nanobooster)/6)),
                               value  = c(data_nanobooster[,6]))
  
  df_nanobooster_perc <- data.frame(type   = c(rep('GFP-NB',nrow(data_nanobooster))),
                                    x      = c(rep(x_vector,nrow(data_nanobooster)/6)),
                                    value  = c(data_nanobooster[,8]))
  
  df_nanobooster_perc1           <- subset(df_nanobooster_perc, x=='x1' | x=='x2' | x=='x3')
  
  #plotting nanobooster data
  
  ggplot(df_nanobooster_perc1,aes(x = x, y = value, colour = type))+
    geom_quasirandom(alpha = 0.8)+
    geom_boxplot(aes(fill = type), alpha = 0.3,outlier.shape = NA)+
    scale_x_discrete(labels = c("x0"="0","x1"="117","x2" = "234","x3"="351","x4"="468","x5"="585"))+
    labs(x='Search radius (nm)',y='% Associated virus particles')+ #where # of associated spots refers to spots which have at least 1 neighbour within R_search
    ylim(0,100)+
    prism_multi()+
    theme(legend.position = "none")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/results association gag-eap45 20190510/ch2-ch1_association_vs_radius_GFPNanobooster_percentage.pdf",width = 90, height = 73, units = "mm")
}

# using estimation graphics to look at how different mutants compare to wildtype
{

unpaired_mean_diff_all <- dabest(df_Gag_percentageAssociated,type,value,
                             idx = c("wt","ALIX" , "EAP45 mt","PTAP mt", "double mt"),
                             paired = FALSE)

plot(unpaired_mean_diff_all, rawplot.ylabel = "% associated Gag spots found")

}

##----------Statistical Analysis ----------##

# Figure 3C: t-test to compare between associated Gag particles with ALIX and FL EAP45 in HeLa cells 
{wildtype  = subset(df_GagAssociated_ALIX_FLEAP45,type=='wt')
ALIX      = subset(df_GagAssociated_ALIX_FLEAP45,type=='ALIX')
x = wildtype %>% select(value)
y = ALIX     %>% select(value)
t.test(x,y)
}

# Figure 4F: ANOVA to compare between associated Gag particles with ALIX and FL EAP45 in HeLa cells 
{
# order factors to compare FL EAP45 vs mutants
df_GagAssociated_EAP45_allMutants$type <- factor(df_GagAssociated_EAP45_allMutants$type,levels = c("wt","EAP45 mt","PTAP mt","double mt"),ordered = TRUE)                                                   

# ANOVA with Dunnett grouping 
res.aov <- aov(value ~ type, df_GagAssociated_EAP45_allMutants)
# BEFORE RUNNING THE NEXT LINE, IMPORT THE MULTCOMP AND DESCTOOLS LIBRARIES
summary(glht(aov, linfct=mcp(type="Dunnett")))

#summary(glht(aov, linfct=mcp(type= c("ALIX - wt >= 0", 
                                     #"EAP45 mt - wt >= 0",
                                     #"double mt - wt >= 0",
                                     #"PTAP mt - wt >= 0"))))

DunnettTest(value ~ type, data = data_stat)
summary(res.aov) # obtain summary from one-way ANOVA
TukeyHSD(res.aov)# interpret F value using the Tukey metric

#1. Check homogeneity of variances 
plot(res.aov,1)
# Use Levene's test to check the homogeneity of the variances in the data
leveneTest(value ~ condition, data = df_rc)
#2. Check normality condition 
plot(res.aov,2)
}

# Figure 4J: t-test to compare associated Gag particles with FL EAP45 and mutants in HAP-1EAP45KO cells. 
{
wildtype  = subset(df_KOcells_Gag_percentageAssociated,type=='wt')
mutant    = subset(df_KOcells_Gag_percentageAssociated,type=='mt')
x = wildtype  %>% select(value)
y = mutant    %>% select(value)
t.test(x,y)
}


