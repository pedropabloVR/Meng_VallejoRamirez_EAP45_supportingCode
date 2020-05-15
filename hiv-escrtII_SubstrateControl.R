# Supplementary code accompanying the manuscript: "Distinct domain requirements for EAP45 in HIV budding, late endosomal recruitment, and cytokinesis"

# Pedro Vallejo Ramirez
# Laser Analytics Group, University of Cambridge
# Created: 2020-01-07
# Updated: 2020-04-07

# Script to create two line plots from the fluorescence intensity cross-sections of two images. 
# The .csv files from a line plot in Fiji are imported, and then plotted. 

# This script can be used to recreate the plots in the Supplementary Figure S1


##----------Import libraries and data files ----------##
library(tidyverse)
library(RColorBrewer)
library(ggbeeswarm)


prism_multi <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = "white"
  color.grid.major = "white"
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
} # plotting function

# read in data
data <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\images for substrate vs no substrate figure\\results_intensity_substrate_NoSubstrate.csv")
data_lineProfile <- read.csv("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\images for substrate vs no substrate figure\\substrate_lineProfile.csv")


##----------Plot data ----------##
ggplot(data, aes(x = type, y = value, color = type,fill = type)) +
  geom_boxplot(alpha = 0.4)+
  geom_beeswarm()+
  labs(y='Fluorescence intensity (a.u.)')+
  prism_multi()+
  theme(axis.title.x=element_blank())
ggsave("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\images for substrate vs no substrate figure\\results_intensity_plot.pdf",width = 80, height = 60, units = "mm")


ggplot(data_lineProfile,aes(x = x, y = value))+
  geom_line()+
  labs(x='Distance (pixels)',y='Fluorescence intensity (a.u.)')+
  prism_multi()
ggsave("D:\\OneDriveNew\\OneDrive - University of Cambridge\\lag\\microscopy work\\hiv-eap45\\images for substrate vs no substrate figure\\lineProfile_Substrate.pdf",width = 80, height = 60, units = "mm")

