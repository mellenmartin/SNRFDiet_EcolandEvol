###################
##=============================
# Sympatric Sierra Nevada Red Fox and Coyote Diet
# Step 2 - phylogenetic plots, bar charts to visualize diet patterns, and calculating diet overlap
#=============================
#### Diet data collected as part of carnivore occurrence/diet study funded by Katie Moriarty + Taal Levi and collected by Rogue Detection Teams
###################

rm(list=ls());gc() #clear the memory

#library(plyr)
library(dplyr)
library(ggplot2)

###  Load from csv
setwd("C:/Users/yourdrivehere/DietData")

diet <- read.csv("ProcessedRedFox_DietGeoDataSimplified_250404.csv") #SNRF diet
yotediet <- read.csv("ProcessedCoyoteRedFoxStudyArea_DietGeoDataSimplified_250404.csv") # CALA diet
names(diet) = tolower(names(diet))
names(yotediet) = tolower(names(yotediet))

##### phylo plots #####
require(dplyr)
require(reshape2)
require(metacoder)
require(ggplot2)
require(ggpubr)
require(lubridate)
require(plyr)

#create functions for metacoder
plots <- function(season, meta.obj, phylum, tree_color_range, seednum, nobs = 0)
{
  set.seed(seednum)
  temp <- meta.obj %>%
    filter_taxa(n_obs > nobs, reassign_obs = F)
  ht <- heat_tree(temp,
                  #node_size=n_obs,
                  node_size = temp$data$carnivore_occ[["total"]],
                  node_size_range = c(0.001, 0.01),
                  node_color = temp$data$carnivore_occ[[season]],
                  node_color_range = tree_color_range,
                  edge_size_range = c(0.005, 0.015),
                  tree_color_range = tree_color_range,
                  layout="davidson-harel",
                  title = phylum,
                  #margin_size = c(-0.1, -0.1, 0.05,0),
                  node_size_axis_label = "All samples",
                  edge_size_axis_label = "No. samples",
                  node_legend_title	= NA,
                  node_color_axis_label	= "Species samples",
                  node_color_digits	= 1,
                  node_size_digits = 1,
                  edge_color_digits	= 1,
                  edge_size_digits = 1)
  return(ht)
}

plots2 <- function(season, meta.obj, phylum, tree_color_range, seednum, nobs = 0)
{
  set.seed(seednum)
  temp <- meta.obj %>%
    filter_taxa(n_obs > nobs, reassign_obs = F)
  ht <- heat_tree(temp,
                  #node_size=n_obs,
                  node_size = temp$data$carnivore_occ[["total"]],
                  node_size_range = c(0.001, 0.01),
                  node_color = temp$data$carnivore_occ[[season]],
                  node_color_range = tree_color_range,
                  edge_size_range = c(0.005, 0.015),
                  tree_color_range = tree_color_range,
                  layout="davidson-harel",
                  title = phylum,
                  node_color_axis_label="",
                  node_size_axis_label="")
  return(ht)
}

key_plot <- function(meta.obj, phylum, tree_color_range, seednum, nobs = 0)
{
  set.seed(seednum)
  temp <- meta.obj %>%
    filter_taxa(n_obs > nobs, reassign_obs = F)
  ht <- heat_tree(temp,
                  node_label = taxon_names, 
                  node_label_size_range = c(0.022, 0.05),
                  #node_size=n_obs,
                  node_size = temp$data$carnivore_occ[["total"]],
                  node_size_range = c(0.005, 0.015),
                  node_color = temp$data$carnivore_occ[["total"]],
                  node_color_range = tree_color_range,
                  edge_size_range = c(0.005, 0.01),
                  tree_color_range = tree_color_range,
                  repel_force=1,
                  overlap_avoidance = 1.2,
                  margin_size = c(0.05, 0.05, 0.0, 0),
                  title = phylum,
                  layout="davidson-harel",
                  node_legend_title	= "",
                  node_size_axis_label = " ",
                  node_color_axis_label=" ",
                  edge_size_axis_label = "No. samples",
                  node_color_digits	=1,
                  node_size_digits = 1,
                  edge_color_digits	= 1,
                  edge_size_digits = 1)
  return(ht)
}

key_plot2 <- function(meta.obj, phylum, tree_color_range, seednum, nobs = 0)
{
  set.seed(seednum)
  temp <- meta.obj %>%
    filter_taxa(n_obs > nobs, reassign_obs = F)
  ht <- heat_tree(temp,
                  node_label = taxon_names, 
                  node_label_size_range = c(0.02, 0.04),
                  #node_size=n_obs,
                  node_size = temp$data$carnivore_occ[["total"]],
                  node_size_range = c(0.001, 0.015),
                  node_color = temp$data$carnivore_occ[["total"]],
                  node_color_range = tree_color_range,
                  edge_size_range = c(0.005, 0.01),
                  tree_color_range = tree_color_range,
                  overlap_avoidance = 1,
                  margin_size = c(0, 0.15, 0, 0),
                  title = phylum,
                  layout="davidson-harel",
                  node_size_axis_label = " ",
                  node_color_axis_label="Occurrence")
  return(ht)
}

####### format SNRF data ######
library(dplyr)
detach(package:plyr)

data.all2 <- diet %>% dplyr::group_by(fieldid.x, 
                                      defecator, 
                                      speciessimpleid.x,
                                      final_scientific,
                                      species, 
                                      genus, 
                                      family, 
                                      order, 
                                      class, 
                                      phylum) %>% 
            dplyr::summarize(reps = sum(replicates),
                      speciesreads = mean(spec_reads),
                      totalreads = mean(total_reads),
                      rra = mean(spec_reads)/mean(total_reads),
                      pidmatch = mean(best_match_pid.mean)) %>% ungroup

# create phylogenetic lineage for the phyloplot
data.all2$lineage <- paste(data.all2$phylum, data.all2$class, data.all2$order, data.all2$family, data.all2$genus, data.all2$species, sep=",") #make lineages line up across metabarcoding data and manual sort data
data.all2$lineage <-  gsub(data.all2$lineage, pattern=",{2,}", replacement="")
write.csv(data.all2, "ProcessedRedFox_PhyloSummaryData_250404.csv", row.names = FALSE)

####### format CALA data ######
data.all3 <- yotediet %>% dplyr::group_by(fieldid.x, 
                                      defecator, 
                                      speciessimpleid.x,
                                      final_scientific,
                                      species, 
                                      genus, 
                                      family, 
                                      order, 
                                      class, 
                                      phylum) %>% 
  dplyr::summarize(reps = sum(replicates),
                   speciesreads = mean(spec_reads),
                   totalreads = mean(total_reads),
                   rra = mean(spec_reads)/mean(total_reads),
                   pidmatch = mean(best_match_pid.mean)) %>% ungroup

# create phylogenetic lineage for the phyloplot
data.all3$lineage <- paste(data.all3$phylum, data.all3$class, data.all3$order, data.all3$family, data.all3$genus, data.all3$species, sep=",") #make lineages line up across metabarcoding data and manual sort data
data.all3$lineage <-  gsub(data.all3$lineage, pattern=",{2,}", replacement="")
write.csv(data.all3, "ProcessedCoyoteRedFoxStudyArea_PhyloSummaryData_250404.csv", row.names = FALSE)

# combine red fox and coyote
data.all <- rbind(data.all2, data.all3)

#vertebrate data
vertdataall <- data.all %>% filter(phylum == "Chordata") # drop the 1-2 incidental red fox samples from low elevation sites on the OR coasts
#vertebrate data
vertdata <- data.all2 %>% filter(defecator == "Vulpes vulpes" & phylum == "Chordata")

#vertebrate data
vertdatayote <- data.all3 %>% filter(defecator == "Canis latrans" & phylum == "Chordata")

##### phyloplots ####
library(plyr)
library("ggsci")
library("ggplot2")
library("gridExtra")
library(grDevices)

set.seed(123)
hcl.colors(6, palette = "teal", alpha = 0.8, rev = TRUE)

vert_color <- c("#D2EEEACC", "#A3D6D6CC", "#7BB8C1CC", "#5599ABCC", "#387893CC", "#2A5676CC")

sk <- vertdataall %>% dplyr::select(fieldid.x, lineage)
sk <- melt(sk, id=c("lineage","fieldid.x"))
sk <- dcast(data=sk, lineage~fieldid.x, fun.aggregate=length)

obj1 <- parse_tax_data(sk, class_cols="lineage", class_sep=",")
obj1$data$carnivore_abund <- calc_taxon_abund(obj1, 'tax_data', cols=unique(vertdataall$fieldid.x))
obj1$data$carnivore_occ <- calc_n_samples(obj1, 'carnivore_abund', cols=unique(vertdataall$fieldid.x), groups=unique(vertdataall[,c("fieldid.x","defecator")])$defecator)
obj1$data$carnivore_occ$total <- obj1$data$carnivore_occ$`Canis latrans` + obj1$data$carnivore_occ$`Vulpes vulpes`

f <- plots("Vulpes vulpes", meta.obj=obj1, phylum="Red fox", tree_color_range=vert_color, seednum=1) # plt subset to fox FOO 
c <- plots("Canis latrans", meta.obj=obj1, phylum="Coyote", tree_color_range=vert_color, seednum=1) # plot subset to coyote FOO
vkey <- key_plot(meta.obj=obj1, phylum="All samples", tree_color_range=vert_color, seednum=1) 

f
c
vkey

indplots <- ggarrange(f,c)
ggarrange(vkey, indplots, nrow = 2, heights = c(1.7, 1))          
ggsave(ggarrange(vkey, indplots, nrow = 2, heights = c(1.5, 1)), 
       filename="prelimresults/SNRFCALADiet_Figure3DietPhylo_250404.jpeg",  
       height=16, width=10, units="in", dpi=600)
######bar charts######
#tabulate number of occurences and mean rra
###########
library(plyr)
library("ggsci")
library("ggplot2")
library("gridExtra")

#summarize number of verts per scat
vertcount2 <- ddply(vertdata, .(fieldid.x), nrow)
#
samples <- unique(vertdata[,c("fieldid.x")])

skfox <- vertdata %>% dplyr::select(fieldid.x, lineage)
skfox <- melt(skfox, id=c("fieldid.x", "lineage"))
skfox <- dcast(data=skfox, lineage~fieldid.x, fun.aggregate=length)

obj <- parse_tax_data(skfox, class_cols="lineage", class_sep=",")
obj$data$carnivore_abund <- calc_taxon_abund(obj, 'tax_data', cols=unique(vertdata$fieldid.x))
obj$data$carnivore_occ <- calc_n_samples(obj, 'carnivore_abund', cols=unique(vertdata$fieldid.x), groups=unique(vertdata[,c("fieldid.x","studyregion_broad")])$studyregion_broad)
obj$data$carnivore_occ$total <- obj$data$carnivore_occ$Fox 


#summarize number of verts per scat
vertcountyote <- ddply(vertdatayote, .(fieldid.x), nrow)
#
samplesyote <- unique(vertdatayote[,c("fieldid.x")])

skyote <- vertdatayote %>% dplyr::select(fieldid.x, lineage)
skyote <- melt(skyote, id=c("fieldid.x", "lineage"))
skyote <- dcast(data=skyote, lineage~fieldid.x, fun.aggregate=length)

objyote <- parse_tax_data(skyote, class_cols="lineage", class_sep=",")
objyote$data$carnivore_abund <- calc_taxon_abund(objyote, 'tax_data', cols=unique(vertdatayote$fieldid.x))
objyote$data$carnivore_occ <- calc_n_samples(objyote, 'carnivore_abund', cols=unique(vertdatayote$fieldid.x), groups=unique(vertdatayote[,c("fieldid.x","studyregion_broad")])$studyregion_broad)
objyote$data$carnivore_occ$total <- objyote$data$carnivore_occ$Fox 

fox.vert <- data.frame(sample=samples$fieldid.x)
fox.vert$vert <- 1
nrow(fox.vert) # 28 snrf scats w/ non fox vertebrates
yote.vert <- data.frame(sample=samplesyote$fieldid.x)
yote.vert$vert <- 1
nrow(yote.vert) # 34 snrf scats w/ non coyote vertebrates

library(wesanderson)
names(wes_palettes)
vert.colors <- wes_palette("Darjeeling2", 4, type = c("discrete"))

#calculate overall frequency of occurrence and relative read abundance
#hcl.colors(7, palette = "Viridis")
"#4B0055" "#353E7C" "#007094" "#009B95" "#00BE7D" "#96D84B" "#FDE333"
require(plyr)
#vert.colors <- c("#4B0055","#007094","#009B95","#00BE7D")
legend.params <- theme(legend.position = c(0.95, 0.025),
                       legend.justification = c("right", "bottom"),
                       legend.box.just = "right",
                       legend.margin = margin(6, 6, 6, 6),
                       legend.title=element_text(size=24), 
                       legend.text=element_text(size=22))

###  Load from csv
setwd("C:/Users/sean.matthews/Documents/MEM/OSU/Projects/Moriarty_SNRFDietPaper/SNRFDiet_CodeData")
#### niche overlap stuff #####
diet <- read.csv("ProcessedRedFox_DietGeoDataSimplified_250404.csv")
yotediet <- read.csv("ProcessedCoyoteRedFoxStudyArea_DietGeoDataSimplified_250404.csv")
names(diet) = tolower(names(diet))
names(yotediet) = tolower(names(yotediet))
sealD <- rbind(diet,yotediet)

evalq(tapply(fieldid.x,defecator,function(x) length(unique(x))),sealD)
length(unique(sealD$speciessimpleid.x))
# 23
#[1] 71
# Per sample
n.tps=evalq(tapply(speciessimpleid.x,fieldid.x,function(x) length(unique(x))),sealD)
hist(n.tps)
mean(n.tps)

seq.depth=evalq(tapply(rep_reads,fieldid.x,sum),sealD)
summary(seq.depth)

sum(seq.depth<50)
names(seq.depth[seq.depth<50])
hist(seq.depth)
# which samples have less than 50 reads
sealD$sample %in% names(seq.depth[seq.depth<50])

sealD50 <- sealD
sealD50$defcommonname <- ifelse(sealD50$defecator=="Vulpes vulpes", "Redfox","Coyote")

Redfox = evalq(defcommonname =="Redfox", sealD50)
Coyote = evalq(defcommonname=="Coyote", sealD50)

IndexNAme= c(unique(sealD50$defcommonname))

length(IndexNAme) # This is how many subsets of data we are using
# Create subsets
t=0
for( i in IndexNAme){
  t=t+1
  nam =paste("Subset",formatC(t,width = 2, format = "d", flag = "0"),sep="")
  temp=evalq(tapply(rep_reads,list(speciessimpleid.x,fieldid.x),sum),sealD50[eval(as.name(i)),])
  temp[is.na(temp)]=0
  temp50=  temp[,apply(temp,2,sum)>50]
  assign(nam,temp50 )
}
## Look at one subset in detail
Subset01 #first subset
dim(Subset01) # number of prey and number of samples

#### Redfox ####
#Look at the counts of prey identified in this subset
Redfox=evalq(tapply(rep_reads,list(speciessimpleid.x,fieldid.x),sum),sealD50[Redfox,])
Redfox[is.na(Redfox)]=0
apply(Redfox,1,sum)
# Convert to proportions
Redfoxprop=prop.table(Redfox,2)
dim(Redfoxprop)
apply(Redfox,1,sum)
# No prey occurences are less than 1%
sum(Redfoxprop < .01)
sum(Redfoxprop== 0)
# Make presence/absence dataset
PA_Redfox= Redfox; PA_Redfox[Redfox>0]=1
apply(PA_Redfox,1,sum)
#POO summary
dim(PA_Redfox)# 66 samples
apply(PA_Redfox,1,sum) # Number of occurences
POO=apply(PA_Redfox,1,sum)/sum(apply(PA_Redfox,1,sum))
apply(PA_Redfox,2,sum)# how many prey per scat
#wPOO summary
dim(PA_Redfox)# 66 samples
wPOO=apply(prop.table(PA_Redfox,2),1,mean)
#RRA summary
RRA=apply(prop.table(as.matrix(Redfox),2),1,mean)
RedfoxSum <- as.data.frame(cbind(POO, RRA, wPOO))
RedfoxSum <- tibble::rownames_to_column(as.data.frame(RedfoxSum), "PreyTaxa")
RedfoxSum$Species <- rep("Redfox", nrow(RedfoxSum))
RedfoxSum$reads <- rowSums(Redfox)
RedfoxSum$final_scientific <- c("Callospermophilus lateralis",
                                "Cervus canadensis",
                                "Chordeiles minor",
                                "Dendragapus fuliginosus",
                                "Junco hyemalis",
                                "Larus sp.",
                                "Lepus americanus",
                                "Marmota flaviventris",
                                "Myodes californicus",
                                "Neotoma cinerea",
                                "Odocoileus sp.",
                                "Peromyscus maniculatus",
                                "Phenacomys intermedius",
                                "Salvelinus fontinalis",
                                "Scapanus sp.",
                                "Sorex sp.",
                                "Tamias sp.",
                                "Thomomys talpoides",
                                "Vireo sp.")

#### coyote ####
#Look at the counts of prey identified in this subset
Coyote=evalq(tapply(rep_reads,list(speciessimpleid.x,fieldid.x),sum),sealD50[Coyote,])
Coyote[is.na(Coyote)]=0
apply(Coyote,1,sum)
# Convert to proportions
Coyoteprop=prop.table(Coyote,2)
dim(Coyoteprop)
apply(Coyote,1,sum)
# No prey occurences are less than 1%
sum(Coyoteprop < .01)
sum(Coyoteprop== 0)
# Make presence/absence dataset
PA_Coyote= Coyote; PA_Coyote[Coyote>0]=1
apply(PA_Coyote,1,sum)
#POO summary
dim(PA_Coyote)# 66 samples
apply(PA_Coyote,1,sum) # Number of occurences
POO=apply(PA_Coyote,1,sum)/sum(apply(PA_Coyote,1,sum))
apply(PA_Coyote,2,sum)# how many prey per scat
#wPOO summary
dim(PA_Coyote)# 66 samples
wPOO=apply(prop.table(PA_Coyote,2),1,mean)
#RRA summary
RRA=apply(prop.table(as.matrix(Coyote),2),1,mean)
CoyoteSum <- as.data.frame(cbind(POO, RRA, wPOO))
CoyoteSum <- tibble::rownames_to_column(as.data.frame(CoyoteSum), "PreyTaxa")
CoyoteSum$Species <- rep("Coyote", nrow(CoyoteSum))
CoyoteSum$reads <- rowSums(Coyote)
CoyoteSum$final_scientific <- c("Accipiter sp.",
                                "Callospermophilus lateralis",
                                "Cervus canadensis",
                                "Chordeiles minor",
                                "Dendragapus fuliginosus",
                                "Junco hyemalis",
                                "Lepus americanus",
                                "Marmota flaviventris",
                                "Mephitis mephitis",
                                "Microtus sp.",
                                "Odocoileus sp.",
                                "Peromyscus maniculatus",
                                "Phenacomys intermedius",
                                "Scapanus sp.",
                                "Sorex sp.",
                                "Tamias sp.",
                                "Tamiasciurus douglasii",
                                "Thomomys talpoides")


library(dplyr)
AllSum <- bind_rows(CoyoteSum,
                    RedfoxSum)

#########################################################

Subset01 <- tibble::rownames_to_column(as.data.frame(Subset01), "PreyTaxa")
write.csv(Subset01, "SNRFCALADiet_GroomedRedfoxPreyMatrix_250404.csv", row.names = FALSE)
Subset02 <- tibble::rownames_to_column(as.data.frame(Subset02), "PreyTaxa")
write.csv(Subset02, "SNRFCALADiet_GroomedCoyotePreyMatrix_250404.csv", row.names = FALSE)

verttab <- ddply(data.all2[data.all2$defecator == "Vulpes vulpes",], .(phylum, class, order, family, genus, species, final_scientific), summarize, freq=length(rra), focc=length(rra)/nrow(fox.vert), totalspecreads = sum(speciesreads), finalrra_all=sum(rra)/nrow(fox.vert), finalrra_vert=sum(rra)/nrow(fox.vert))
verttab <- left_join(verttab, RedfoxSum, by = c("final_scientific"))
verttab$c <- factor(verttab$class)
verttab$final_scientific <- as.character(verttab$final_scientific)
verttab$sciname <- ifelse(verttab$final_scientific == " ", verttab$family, verttab$final_scientific)
verttab$sciname <- ifelse(verttab$sciname == "", verttab$order, verttab$sciname)
verttab$sciname <- factor(verttab$sciname, levels=verttab[order(verttab$finalrra_vert, decreasing = F),]$sciname) #order by vert RRA
write.csv(verttab, "RedFoxDiet_VertDietDataSummary_250404.csv", row.names = FALSE)
verttab <- verttab %>% filter(!sciname == "",
                              !sciname == "Passeriformes")

verttabyote <- ddply(data.all3[data.all3$defecator == "Canis latrans",], .(phylum, class, order, family, genus, species, final_scientific), summarize, freq=length(rra), focc=length(rra)/nrow(fox.vert), totalspecreads = sum(speciesreads), finalrra_all=sum(rra)/nrow(fox.vert), finalrra_vert=sum(rra)/nrow(fox.vert))
verttabyote <- left_join(verttabyote, CoyoteSum, by = c("final_scientific"))
verttabyote$c <- factor(verttabyote$class)
verttabyote$final_scientific <- as.character(verttabyote$final_scientific)
verttabyote$sciname <- ifelse(verttabyote$final_scientific == " ", verttabyote$family, verttabyote$final_scientific)
verttabyote$sciname <- ifelse(verttabyote$sciname == "", verttabyote$order, verttabyote$sciname)
verttabyote$sciname <- factor(verttabyote$sciname, levels=verttabyote[order(verttabyote$finalrra_vert, decreasing = F),]$sciname) #order by vert RRA
write.csv(verttabyote, "CoyoteDietRedFoxStudyArea_VertDietDataSummary_250404.csv", row.names = FALSE)
verttabyote <- verttabyote %>% filter(!sciname == "",
                                      !sciname == "Passeriformes")

#### load in phylopics ####
library(rphylopic)
hawk <- pick_phylopic(name = "Accipitridae", n = 8, view = 8) # pick #4 #Sharon Wegner Larsen
hare <- pick_phylopic(name = "Lepus americanus", n = 4, view = 4) # pick #2 #Margot Michaud 
squir <- pick_phylopic(name = "Urocitellus beldingi", n = 4, view = 4) # no callo, pick other ground squirrel Patricia Holroyd
mouse <- pick_phylopic(name = "Peromyscus leucopus", n = 4, view = 4) # Nina Skinner
vole <- pick_phylopic(name = "Cricetidae", n = 6, view = 6) # 5 # Callum Le Lay
mole <- pick_phylopic(name = "Talpidae", n = 6, view = 6) # 3 # Birgit lang
gull <- pick_phylopic(name = "Larus", n = 6, view = 6) # 3 # Andy wilson
vireo <- pick_phylopic(name = "Vireo", n = 6, view = 6) # caleb gordon

shrew <- pick_phylopic(name = "Sorex", n = 6, view = 6) # Becky barnes
elk <- pick_phylopic(name = "Cervus", n = 6, view = 6) #2 Steven Traver
junco <- pick_phylopic(name = "Passerellidae", n = 6, view = 6) # Ferran Sayol
chip <- pick_phylopic(name = "Tamias", n = 6, view = 6) # Chloe Schmidt
deer <- pick_phylopic(name = "Odocoileus", n = 6, view = 6) # Tracy Heath
trout <- pick_phylopic(name = "Salvelinus", n = 6, view = 6) # xgirouxb
redvole <- pick_phylopic(name = "Myodes", n = 6, view = 6) # 
grouse <- pick_phylopic(name = "Phasianinae", n = 10, view = 10) #xgirouxb
nighthawk <- pick_phylopic(name = "Caprimulgidae", n = 6, view = 6) # Ferran Sayol
woodrat <- pick_phylopic(name = "Cricetidae", n = 20, view = 6) # nsvitek
marm <- pick_phylopic(name = "Marmota", n = 6, view = 6) # Michael Keesey
doug <- pick_phylopic(name = "Tamiasciurus", n = 6, view = 6) # Chloe Schmidt
strsku <- pick_phylopic(name = "Mephitis", n = 6, view = 6) # Margot Michaud
pocket <- pick_phylopic(name = "Geomyidae", n = 6, view = 6) # Xavier Jenkins
coyote <- pick_phylopic(name = "Canis latrans", n = 6, view = 6) # Margot Michaud
fox <- pick_phylopic(name = "Vulpes", n = 6, view = 6) #1, Margot Michaud

library(dplyr)
library(forcats)
verttab <- verttab %>% dplyr::filter(!final_scientific == "Canis latrans") %>% dplyr::mutate(final_scientific = forcats::fct_reorder(final_scientific, wPOO))
verttabyote <- verttabyote %>% dplyr::mutate(final_scientific = forcats::fct_reorder(final_scientific, wPOO))

speciesfocc <- ggplot(verttab, aes(x=final_scientific,
                                  y=focc,
                                  fill=c))+ 
  geom_bar(stat="identity", width = 0.7) + 
  theme_bw(base_size=15) + #coord_flip() +
  ylim(0,0.43)+
  scale_y_reverse() +
  scale_x_discrete(name = "", position = "top") +
  coord_flip () +
  scale_fill_manual(name="Prey class", values=c("#ABDDDE",
                                                "#D69C4E",
                                                "#046C9A")) + 
  legend.params +
  ylab(" ") +
  annotate("text", x=19, y=0.34, label= "A", size=16, parse=TRUE)+
  add_phylopic(img = fox, x = 3.5, y = 0.3, ysize = 5)+
  #xlab("species name") +
  theme(axis.text.x = element_text(size=24),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.title = element_text(size=24),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(0.2,-1.6,0.2,1), 'lines'),
        legend.position = "none")
speciesfocc

speciesfoccyote <- ggplot(verttabyote, aes(x=final_scientific,
                                   y=focc,
                                   fill=c))+ 
  #ggtitle("Montane prey items")+
  geom_bar(stat="identity", width = 0.7) + 
  theme_bw(base_size=15) + #coord_flip() +
  ylim(0,0.43)+
  scale_y_reverse() +
  scale_x_discrete(name = "", position = "top") +
  coord_flip () +
  scale_fill_manual(name="Prey class", values=c("#D69C4E",
                                                "#046C9A")) + 
  legend.params +
  ylab("Frequency of occurrence") +
  annotate("text", x=18, y=0.38, label= "C", size=16, parse=TRUE)+
  add_phylopic(img = coyote, x = 3.5, y = 0.28, ysize = 6.0)+
  
  #xlab("species name") +
  theme(axis.text.x = element_text(size=24),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.title = element_text(size=24),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(0.2,-1.6,0.2,1), 'lines'),
        legend.position = "none")
speciesfoccyote

specieswpoo <- ggplot(verttab,aes(x=final_scientific, y=wPOO, fill=c)) + 
  #rra <- ggplot(verttab, aes(x=sciname, y=finalrra_all, fill=c)) + 
  geom_bar(stat="identity", width = 0.7) + theme_bw(base_size=15) + coord_flip() +
  scale_fill_manual(name="Prey class", values=c("#ABDDDE",
    "#D69C4E",
    "#046C9A")) + 
  legend.params +
  ylab(" ") +
  xlab("")+
  annotate("text", x=19, y=0.42, label= "B", size=16, parse=TRUE)+
  add_phylopic(img = squir, x = 19, y = 0.15, ysize = 1.2)+
  add_phylopic(img = mouse, x = 18, y = 0.13, ysize = 0.85)+
  add_phylopic(img = pocket, x = 17, y = 0.12, ysize = 0.85)+
  add_phylopic(img = chip, x = 16, y = 0.11, ysize = 1.0)+
  add_phylopic(img = redvole, x = 15, y = 0.09, ysize = 0.85)+
  add_phylopic(img = deer, x = 14, y = 0.09, ysize = 1)+ 
  add_phylopic(img = hare, x = 13, y = 0.08, ysize = 0.9)+ 
  add_phylopic(img = mole, x = 12, y = 0.075, ysize = 0.7)+
  add_phylopic(img = elk, x = 11, y = 0.075, ysize = 1.1)+
  add_phylopic(img = vole, x = 10, y = 0.065, ysize = 0.55)+ 
  add_phylopic(img = trout, x = 9, y = 0.06, ysize = 0.7)+ 
  add_phylopic(img = nighthawk, x = 8, y = 0.06, ysize = 0.7)+ 
  add_phylopic(img = shrew, x = 7, y = 0.05, ysize = 0.7)+ 
  add_phylopic(img = marm, x = 6, y = 0.05, ysize = 0.85)+ 
  add_phylopic(img = woodrat, x = 5, y = 0.03, ysize = 0.8)+ 
  add_phylopic(img = junco, x = 4, y = 0.03, ysize = 0.7)+ 
  add_phylopic(img = grouse, x = 3, y = 0.03, ysize = 0.8)+ 
  add_phylopic(img = vireo, x = 2, y = 0.02, ysize = 0.8)+ 
  add_phylopic(img = gull, x = 1, y = 0.02, ysize = 0.9)+ 
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(hjust = 0.5,
                                   size = 24,
                                   margin(t = 0, r = 20, b = 0, l = 0, unit = "lines")),
        axis.title = element_text(size = 24),
        plot.title = element_text(hjust = 0.5, size = 20),
        plot.margin = unit(c(0.2,1,0.2,-1.5), 'lines'))+
  ylim(0,0.43)

specieswpoo

specieswpooyote <- ggplot(verttabyote,aes(x=final_scientific, y=wPOO, fill=c)) + 
  #rra <- ggplot(verttab, aes(x=sciname, y=finalrra_all, fill=c)) + 
  geom_bar(stat="identity", width = 0.7) + theme_bw(base_size=15) + coord_flip() +
  scale_fill_manual(name="Prey class", values=c("#D69C4E",
                                                "#046C9A")) + 
  legend.params +
  ylab("Weighted percent of occurrence") +
  xlab("")+
  annotate("text", x=18, y=0.42, label= "D", size=16, parse=TRUE)+
  add_phylopic(img = hare, x = 18, y = 0.25, ysize = 1.2)+ 
  add_phylopic(img = squir, x = 17, y = 0.22, ysize = 1.0)+ 
  add_phylopic(img = deer, x = 16, y = 0.16, ysize = 0.85)+
  add_phylopic(img = chip, x = 15, y = 0.15, ysize = 1.1)+
  add_phylopic(img = marm, x = 14, y = 0.06, ysize = 0.8)+ 
  add_phylopic(img = elk, x = 13, y = 0.06, ysize = 1.1)+ 
  add_phylopic(img = pocket, x = 12, y = 0.05, ysize = 0.7)+ 
  add_phylopic(img = shrew, x = 11, y = 0.045, ysize = 0.55)+ 
  add_phylopic(img = strsku, x = 10, y = 0.045, ysize = 0.9)+
  add_phylopic(img = grouse, x = 9, y = 0.045, ysize = 0.8)+ 
  add_phylopic(img = mole, x = 8, y = 0.045, ysize = 0.7)+
  add_phylopic(img = vole, x = 7, y = 0.04, ysize = 0.55)+ 
  add_phylopic(img = hawk, x = 6, y = 0.025, ysize = 0.9)+ 
  add_phylopic(img = doug, x = 5, y = 0.025, ysize = 0.8)+
  add_phylopic(img = vole, x = 3.9, y = 0.034, ysize = 0.55)+ 
  add_phylopic(img = junco, x = 3, y = 0.03, ysize = 0.8)+ 
  add_phylopic(img = nighthawk, x = 2, y = 0.03, ysize = 0.7)+ 
  add_phylopic(img = mouse, x = 1, y = 0.03, ysize = 0.7)+ 
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(hjust = 0.5,
                                   size = 24,
                                   margin(t = 0, r = 20, b = 0, l = 0, unit = "lines")),
        axis.title = element_text(size = 24),
        plot.title = element_text(hjust = 0.5, size = 20),
        plot.margin = unit(c(0.2,1,0.2,-1.5), 'lines'))+
  ylim(0,0.43)

specieswpooyote

##### save the plots 
foccwpoo <- ggarrange(speciesfocc,specieswpoo, 
                      widths = c(1.1, 1.8))
foccwpooyote <- ggarrange(speciesfoccyote,specieswpooyote, 
                      widths = c(1.1, 1.8))

ggsave(ggarrange(foccwpoo,foccwpooyote,nrow = 2, 
                             heights = c(1, 1)), 
                   filename=paste("prelimresults/SNRFCALADiet_FoowPOO_250404.jpeg", sep=""), 
                   height=17, width=19, units="in", dpi=600)


#### sample summaries ####
#### taxa summaries ###
# class
classsum <- data.all2 %>% dplyr::group_by(fieldid.x,
                                          class) %>% dplyr::summarize(count = nrow(fieldid.x)) %>% ungroup
table(classsum$class)
#Actinopterygii       Aves       Mammalia      
#         1              5            26 


# class
classsum <- data.all3 %>% dplyr::group_by(fieldid.x,
                                          class) %>% dplyr::summarize(count = nrow(fieldid.x)) %>% ungroup
table(classsum$class)
#      Aves       Mammalia      
#        4            32 
