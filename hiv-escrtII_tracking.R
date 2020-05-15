# Supplementary code accompanying the manuscript: "Distinct domain requirements for EAP45 in HIV budding, late endosomal recruitment, and cytokinesis"

# Pedro Vallejo Ramirez
# Laser Analytics Group, University of Cambridge
# Created: 2019-10-25
# Updated: 2020-04-07

# This script takes the results from the co-moving frame analysis in Matlab (the distances between
# the centroids of the EAP45 and Gag particles over time) and produces:

# - Line profile plots of the distance between the particles over time
# - A histogram of the percentage of the total observation time the particles spend at given distances
#   from each other. 
# - A kernel density plot showing the maximum consecutive time the EAP45 and Gag particle spend 
#   closer than 100 nm of each other. 

# This script can be used to recreate the plots in following figures:

# - Figure 5B (Line plots)
# - Figure 5C (Kernel density plot)
# - Figure 5D (Percentages of total observation time as a function of binned distances between EAP45 and Gag)
# - Supplementary Figure S4 (Line plots)

##----------Import libraries and data files ----------##

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggbeeswarm)

# Set random number seed to get reproducibe plots
set.seed(2019) 
time <- 5.22/60; # conversion factor for time

# Use custom plotting functions
prism_multi  <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = "white"
  color.axis.text = "black"
  color.axis.title = "black"
  color.title = palette[9]
  font_size <- 10
  title_size <- 12
  
  # Begin construction of chart
  theme_bw(base_size=11) +
    
    # Format the legend, uncomment position if need be
    theme(legend.text = element_text(size=font_size, color=color.axis.title)) +
    theme(legend.title =element_text(size=font_size, color=color.axis.title))+
    theme(legend.position="none")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=font_size, vjust=1.25, hjust = 0.4, face = "bold")) +
    theme(axis.text.x=element_text(size=font_size,color=color.axis.text, margin = unit(c(0.1, 0.5, 0.3, 0.5), "cm"))) +
    theme(axis.text.y=element_text(size=font_size,color=color.axis.text, margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm"))) +
    theme(axis.title.x=element_text(size=font_size,color=color.axis.title, 
                                    vjust=1.25, margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))) +
    theme(axis.title.y=element_text(size=font_size,color=color.axis.title, 
                                    vjust=-4, hjust = 0.5, margin = unit(c(0.2, 0.2, 0.2, 
                                                                           0.2), "cm"))) 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 
  #0.25)) 
}

prism_multi2 <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = "white"
  color.axis.text = "black"
  color.axis.title = "black"
  color.title = palette[9]
  font_size <- 10
  title_size <- 12
  
  # Begin construction of chart
  theme_bw(base_size=11) +
    
    # Format the legend, uncomment position if need be
    theme(legend.text = element_text(size=font_size, color=color.axis.title)) +
    theme(legend.title =element_text(size=font_size, color=color.axis.title))+
    theme(legend.position="none")+
    #theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=font_size, vjust=1.25, hjust = 0.4, face = "bold")) +
    theme(axis.text.x=element_text(size=font_size,color=color.axis.text, margin = unit(c(0.1, 0.5, 0.3, 0.5), "cm"))) +
    theme(axis.text.y=element_text(size=font_size,color=color.axis.text, margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm"))) +
    theme(axis.title.x=element_text(size=font_size,color=color.axis.title, 
                                    vjust=1.25, margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))) +
    theme(axis.title.y=element_text(size=font_size,color=color.axis.title, 
                                    vjust=-4, hjust = 0.5, margin = unit(c(0.2, 0.2, 0.2, 
                                                                           0.2), "cm"))) 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.25)) 
} 

# Helper functions to plot the percentage of total observation time spent by the EAP45 at given distances from the Gag particle 
is.even      <- function(x) {
  x %% 2 == 0
} # function to check for even digits

percentage_histogramValues <- function(values){
  # test to get percentage of frames from 0-2, 2-4, 4-6
  
  max_value = floor(max(values$distance)+1)
  even <- is.even(max_value)
  if (even == FALSE) {max_value <- max_value + 1}
  
  my.bin.width <-2
  foo   <-hist(values$distance,breaks=seq(0,max_value,by=my.bin.width))
  total <- sum(foo$counts)
  perc <- (foo$counts/total)*100
  df <- data.frame(frame = foo$breaks[1:length(foo$breaks)-1],
                   value = perc)
  plot(foo$breaks[1:length(foo$breaks)-1],perc, 
       xlab="Distance in pixels", ylab="Percentage of total frames")
  return(df)
  
} # function to get y axis values from a histogram with binwidth 2, given a 2-column table 


# Combined histogram of all the results for the 14 ROIs analyzed so far
# Live8
live8_1_ROI1 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live8_slide1well2_1_ROI1\\live8_slide1well2_1_488_registered_ROI1distance_vs_time_results.csv")
live8_1_ROI2 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live8_slide1well2_1_ROI2\\live8_slide1well2_1_488_registered_ROI2distance_vs_time_results.csv")
live8_4_ROI1 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live8_slide1well2_4_ROI1\\live8_slide1well2_4_488_registered_ROI1distance_vs_time_results.csv")
live8_4_ROI3 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live8_slide1well2_4_ROI3\\live8_slide1well2_4_488_registered_ROI3distance_vs_time_results.csv")

# Live6
live6_12_ROI1 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live6_slide1well3_12_ROI1\\live6_slide1well3_12_488_registered_ROI1distance_vs_time_results.csv");
live6_12_ROI3 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live6_slide1well3_12_ROI3\\live6_slide1well3_12_488_registered_ROI3distance_vs_time_results.csv");
live6_12_ROI5 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live6_slide1well3_12_ROI5\\live6_slide1well3_12_488_registered_ROI5distance_vs_time_results.csv");
live6_12_ROI6 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live6_slide1well3_12_ROI6\\live6_slide1well3_12_488_registered_ROI6distance_vs_time_results.csv");

# Live7
live7_4_ROI1 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live7_slide1well2_4_ROI1\\live7_slide1well2_4_488_registered_ROI1distance_vs_time_results.csv")
live7_4_ROI2 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live7_slide1well2_4_ROI2\\live7_slide1well2_4_488_registered_ROI2distance_vs_time_results.csv")
live7_5_ROI1 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live7_slide1well2_5_ROI1\\live7_slide1well2_5_488_registered_ROI1distance_vs_time_results.csv")
live7_7_ROI1 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live7_slide1well2_7_ROI1\\live7_slide1well2_7_488_registered_ROI1distance_vs_time_results.csv")
live7_8_ROI3 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live7_slide1well2_8_ROI3\\live7_slide1well2_8_488_registered_ROI3distance_vs_time_results.csv")
live7_8_ROI4 <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\live7_slide1well2_8_ROI4\\live7_slide1well2_8_488_registered_ROI4distance_vs_time_results.csv")

##----------Generate plots ----------##

# Figure 5D Percentage of the total observation time the EAP45 spends within given distances from the Gag particle
{
# call the percentage histogram function in a for loop and fill all the percentage values for a given bin width

# List of all open data sets with the prefix "live"
fileList = ls(pattern = "live")

for (i in 1:length(fileList)){
  file = get(fileList[i])
  if (i == 1){aggregated_percentage <- percentage_histogramValues(file)}
  else{aggregated_percentage = rbind(aggregated_percentage,percentage_histogramValues(file))}
    
}

max_distance = 14;

aggregated_percentage2 <- subset(aggregated_percentage, frame< max_distance)

# Plotting the percentage of the total observation time each particle EAP45 particle spends a given distance from the centroid of a virus particle
ggplot(aggregated_percentage2,aes(x = frame,y = value))+
  geom_point(alpha = 0.8,colour = "#00BFC4",size = 0.3)+
  geom_boxplot(aes(group = frame), alpha = 0.2,outlier.shape = NA,colour = "#00BFC4",fill = "#00BFC4",fatten=NULL)+
  scale_x_continuous(name = "Distance from centroid of virus particle (pixels)",
                     breaks=seq(0,max_distance,2),labels = c("0-2","2-4","4-6","6-8","8-10","10-12","12-14","14-16"))+
  ylab("Percentage of total observation time")+
  theme(legend.position = "none")+
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 1.5, size = 1, linetype = "solid",colour = "#00BFC4")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  prism_multi()
ggsave("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\PercentageTime_distance_histogram_14ROIs_withMean.pdf",width = 104, height = 60, units = "mm")

# concatenate results into single data frame
dfMerged <- do.call("rbind", list(live8_1_ROI1, live8_1_ROI2,live8_4_ROI1,live8_4_ROI3,
                                  live6_12_ROI1,live6_12_ROI3,live6_12_ROI5,live6_12_ROI6,
                                  live7_4_ROI1,live7_4_ROI2,live7_5_ROI1,live7_7_ROI1,live7_8_ROI3,live7_8_ROI4))

# plot histogram of merged data
ggplot(dfMerged,aes(distance))+
  geom_histogram(bins = 8,alpha = 0.6,colour = "black",binwidth = 2,fill = "#00BFC4")+
  #geom_density(aes(y=2 * ..count..),alpha = 0.2, fill = "#00BFC4",color = "black")+
  xlab("Distance (pixels)")+
  ylab("Number of candidates")+
  prism_multi()
ggsave("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\Concatenated_distance_histogram_14ROIs.pdf",width = 80, height = 60, units = "mm")
}

# Figure 5C: Consecutive time the EAP45 particle spends within 117 nm of the Gag particle
{
# scatter plot for maximum number of consecutive frames
consecutive <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\max_consecutive_frames.csv")
consecutive <- read.csv("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/tracking data/max_consecutive_frames.csv")
ggplot(consecutive,aes(x = filler, y = time/60)) + 
  geom_quasirandom(alpha = 0.8, size = 1,dodge.width = 0.8)+
  geom_violin(alpha = 0.4,fill = "#00BFC4")+
  ylab("Consecutive time period (min)")+
  prism_multi()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\MaxConsecutiveFrames_14ROIs.pdf",width = 80, height = 60, units = "mm")
ggsave("/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/tracking data/MaxConsecutiveFrames_14ROIs_violin.pdf",width = 75, height = 55, units = "mm")
}

# Figure 5B Line Profiles 
{
# Helper functions to generate line profiles of the distance between EAP45 and Gag particles versus time 
extractName <- function(v1) 
{
  v1_string <-deparse(substitute(v1))
  return(v1_string)
}

makeLineProfile <- function(inputData)
{
  # use a ggline plot0
  df_lineProfile <- data.frame(frame = inputData[,1]*time,
                               distance = inputData[,2],
                               threshold = rep(2,nrow(inputData)))
  filename <- extractName(sys.call())
  root_dir <- "D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\"
  path_dir <- paste(root_dir,"distance_vs_frame_",filename,".pdf",sep="")
  
  ggplot(df_lineProfile,aes(x = frame))+
    geom_line(aes(y = distance),size = 0.5)+
    geom_line(aes(y = threshold),size = 0.5,linetype="dashed")+
    labs(x='Time (min) ',y='Distance (pixels)')+
    prism_multi()
  #ggsave(path_dir,width = 50, height = 30, units = "mm") # tried to save the profiles automatically but it was tricky - left as an exercise to the reader
  
}

# Call line profile function for each of the data sets individually and save the output
makeLineProfile(live6_12_ROI6)
ggsave("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\tracking data\\distance_vs_frame_live6_12_ROI6.pdf",width = 80, height = 60, units = "mm")
}

