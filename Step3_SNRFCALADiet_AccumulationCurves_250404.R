###################
##=============================
# Sympatric Sierra Nevada Red Fox and Coyote Diet
# Step 3 - # produce rarefraction/accumulation curves for  scat species diversity 
#=============================
#### details on data collection and source code scrubbed for double-blind review
###################

rm(list=ls());gc() #clear the memory

#library(plyr)
library(dplyr)
library(ggplot2)

###  Load from csv
setwd("C:/Users/yourwdhere")

rm(list=ls());gc() #I like to start with these commands to clear the memory

#####
#load packages
require(iNEXT)
require(vegan)
require(ggplot2)
require(plyr)
require(tidyr) #for spread and gather

#####
data.all <- read.csv("ProcessedRedFox_PhyloSummaryData_250328.csv")
foxdata <- data.all
foxdata <- foxdata[!duplicated(foxdata[,c("lineage","fieldid.x")]),]
fox.wide <- spread(foxdata[,c("fieldid.x","lineage","speciesreads")], lineage, speciesreads, fill=0) #dims: 196 rows, 42 species, should be 128 rows
fox.wide[is.na(fox.wide)] <- 0 #replace NAs with 0
fox.wide[,2:23][fox.wide[,2:23] > 0] <- 1

###############
#create and plot curves
data.inext <- list(Fox = t(fox.wide[,!names(fox.wide) %in% c("fieldid.x"),]))
hcl.colors(5, palette = "ag_GrnYl")
"#255668" "#007E7D" "#17A77E" "#82CC6C" "#EDEF5C"

out.inext <- iNEXT(x=data.inext, datatype="incidence_raw", endpoint=100) #data needs to be rows = species, columns = samples
p <- ggiNEXT(out.inext, color.var = "Assemblage",type=1) + 
  xlab("Number of scat samples") + 
  ylab("Taxonomic richness") +
  xlim(0,100)+
  #ylim(0,70)+
  scale_color_manual(values=c("#404080"))+
  scale_fill_manual(values = c("#404080"))+
  ggtitle("Red fox")+
  #theme(axis.title = element_text(size = 16))+
  #theme(legend.position="right")+
  #theme(legend.title = element_text(colour="black", size=16))+
  #theme(legend.text = element_text(colour="black", size = 16))+
  #theme(axis.text=element_text(size=10))
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_line(color = "lightgrey"),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))
p
ggsave("prelimresults/SNRFDiet_Accum_250402.jpeg", width=8,height=6, units = 'in', dpi= 300)

### coyote data 
yotedata <- read.csv("ProcessedCoyoteRedFoxStudyArea_PhyloSummaryData_250328.csv")
yotedata <- yotedata[!duplicated(yotedata[,c("lineage","fieldid.x")]),]
yote.wide <- spread(yotedata[,c("fieldid.x","lineage","speciesreads")], lineage, speciesreads, fill=0) #dims: 196 rows, 42 species, should be 128 rows
yote.wide[is.na(yote.wide)] <- 0 #replace NAs with 0
yote.wide[,2:21][yote.wide[,2:21] > 0] <- 1

###############
#create and plot curves
data.inext2 <- list(Coyote = t(yote.wide[,!names(yote.wide) %in% c("fieldid.x"),]))
hcl.colors(5, palette = "ag_GrnYl")
"#255668" "#007E7D" "#17A77E" "#82CC6C" "#EDEF5C"
"#69b3a2", "#404080"
out.inext2 <- iNEXT(x=data.inext2, datatype="incidence_raw", endpoint=100) #data needs to be rows = species, columns = samples
yote <- ggiNEXT(out.inext2, color.var = "Assemblage",type=1) + 
  xlab("Number of scat samples") + 
  ylab("Taxonomic richness") +
  xlim(0,100)+
  ylim(0,65)+
  scale_color_manual(values=c("#69b3a2"))+
  scale_fill_manual(values = c("#69b3a2"))+
  ggtitle("Coyote")+
  #theme(axis.title = element_text(size = 16))+
  #theme(legend.position="right")+
  #theme(legend.title = element_text(colour="black", size=16))+
  #theme(legend.text = element_text(colour="black", size = 16))+
  #theme(axis.text=element_text(size=10))
  theme(axis.text = element_text(size=18),
        axis.title= element_text(size = 22),
        panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_line(color = "lightgrey"),
        panel.background = element_blank(),
        strip.text.x=element_text(size=18, margin = margin(2,2,2,2)),
        legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        axis.line = element_line(colour = "black"))
yote
ggsave("prelimresults/SNRFDiet_Accum_250402.jpeg", width=8,height=6, units = 'in', dpi= 300)


ggsave(ggarrange(p,
                 yote,
                 nrow=1, ncol = 2), 
       filename=paste("prelimresults/SNRFCALA_AccuCurves_", Sys.Date(), ".jpg", sep=""), 
       height=4, width=10, units="in", dpi=600)
