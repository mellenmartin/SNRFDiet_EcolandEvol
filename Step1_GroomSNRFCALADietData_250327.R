###################
##=============================
# Sympatric Sierra Nevada Red Fox and Coyote Diet
# Step 1 - grooming data, prepping for following summaries and analysis
#=============================
#### Diet data collected as part of carnivore occurrence/diet study funded by Katie Moriarty + Taal Levi and collected by Rogue Detection Teams
###################
rm(list=ls());gc() #clear the memory

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("C:/Users/yourdrivehere/DietData")

diet <- read.csv("SNRFCoyote_diet_250327.csv")
names(diet) = tolower(names(diet))
diet$best_match_pid.mean <- as.numeric(diet$best_match_pid.mean)
diet <- diet %>% pivot_longer(., cols = repa:repc,
                              names_to = "rep",
                              values_to = 'counts',
                              values_drop_na = FALSE)
diet2 <- diet %>% 
  dplyr::filter(defecator == 'Vulpes vulpes',
                best_match_pid.mean >= 98, # threshold of 98% for PID
                !final_scientific == "Vulpes vulpes") %>% # filter out defecator data
                #!fieldid == 'DESC_0190_13',
                #!fieldid == 'DESC_0190_18',
                #!fieldid == 'DESC_0190_7',
                #!fieldid == 'DESC_0682_14',
                #!fieldid == 'KLAM_0001_7') %>% 
          group_by(labid,
                   fieldid,
                   rep,
                   speciessimpleid,
                   defecator,
                   final_scientific,
                   species,
                   genus,
                   family,
                   order,
                   class,
                   final_common,
                   genus_group,
                   family_group,
                   order_group,
                   class_group) %>%                
                            dplyr::summarise(counts = sum(counts),
                                   best_match_pid.mean = mean(best_match_pid.mean))  %>% dplyr::filter(!speciessimpleid == "") %>% dplyr::filter(!speciessimpleid == "Canis latrans")
sampcount <- unique(diet2$fieldid)
diet2$phylum <- rep("Chordata") # add phylum to support phyloplots (later step)
diet2$key <- paste0(diet2$fieldid,"_",diet2$rep,diet2$speciessimpleid) # our merging key
diet2$counts[is.na(diet2$counts)] <- 0 # if NA/no data, fill in w/ 0

## Calculate the percent reads for a prey item per scat
total_reads <- diet2 %>% group_by(fieldid) %>% summarize_at(vars(counts), sum) #total reads per scat
spec_reads <- diet2 %>% group_by(fieldid,speciessimpleid) %>% summarize_at(vars(counts), sum) #total reads per species per scat
key_reads <- diet2 %>% group_by(fieldid,rep) %>% summarize_at(vars(counts), sum) # merging key

# fieldid, and rep
reads <- merge(total_reads,key_reads, by="fieldid")
reads <- merge(reads,spec_reads,by="fieldid")
head(reads)

colnames(reads) <- c("fieldid","total_reads","rep", "rep_reads", "speciessimpleid","spec_reads")
reads$p_reads <- reads$spec_reads/reads$total_reads # % reads of a given species relative to the whole scat
reads$key <- paste0(reads$fieldid,"_",reads$rep,reads$speciessimpleid) # make our key to remerge

# merge back with the larger dataframe
db <- merge(diet2,reads, by="key", all.y = F)

# write this to csv
write.csv(db,"ProcessedRedFox_DietGeoData_250404.csv", row.names = FALSE)

# d1 now has the number of times 
library(plyr)
names(db)
d1 = ddply(db,
           .(fieldid.x,
             labid,
             rep.x,
             defecator,
             speciessimpleid.x,
             final_scientific,
             species,
             genus,
             family,
             order,
             class,
             phylum,
             final_common,
             genus_group,
             family_group,
             order_group,
             class_group,
             best_match_pid.mean,
             total_reads,
             rep_reads,
             spec_reads,
             p_reads), summarise,
             replicates = ifelse(rep_reads == 0, 0, 1)) 

### filter data above by additional thresholds which are more robust/selective for analysis and data viz ####
draft.dat <- d1[which(d1$p_reads>=0.01),] # remove items with less than 1 percent of reads within a replicate
sampcount <- unique(draft.dat$fieldid.x)
write.csv(draft.dat,"ProcessedRedFox_DietGeoDataSimplified_250404.csv", row.names = FALSE)

#### repeat with coyote diet in snrf study area #####
diet3 <- diet %>% 
  dplyr::filter(defecator == 'Canis latrans',
                best_match_pid.mean >= 98, # threshold of 98% for PID
                !final_scientific == "Canis latrans") %>% # filter out defecator data
  #!fieldid == 'DESC_0190_13',
  #!fieldid == 'DESC_0190_18',
  #!fieldid == 'DESC_0190_7',
  #!fieldid == 'DESC_0682_14',
  #!fieldid == 'KLAM_0001_7') %>% 
  group_by(labid,
           fieldid,
           rep,
           speciessimpleid,
           defecator,
           final_scientific,
           species,
           genus,
           family,
           order,
           class,
           final_common,
           genus_group,
           family_group,
           order_group,
           class_group) %>%                
  dplyr::summarise(counts = sum(counts),
                   best_match_pid.mean = mean(best_match_pid.mean))%>% dplyr::filter(!speciessimpleid == "") %>% dplyr::filter(!speciessimpleid == "Canis latrans")

sampcount2 <- unique(diet3$fieldid)
diet3$phylum <- rep("Chordata")
diet3$key <- paste0(diet3$fieldid,"_",diet3$rep,diet3$speciessimpleid) # our merging key
diet3$counts[is.na(diet3$counts)] <- 0

## Calculate the percent reads for a prey item per scat
total_reads <- diet3 %>% group_by(fieldid) %>% summarize_at(vars(counts), sum) #total reads per scat
spec_reads <- diet3 %>% group_by(fieldid,speciessimpleid) %>% summarize_at(vars(counts), sum) #total reads per species per scat
key_reads <- diet3 %>% group_by(fieldid,rep) %>% summarize_at(vars(counts), sum) # merging key

# fieldid, and rep
reads <- merge(total_reads,key_reads, by="fieldid")
reads <- merge(reads,spec_reads,by="fieldid")
head(reads)

colnames(reads) <- c("fieldid","total_reads","rep", "rep_reads", "speciessimpleid","spec_reads")
reads$p_reads <- reads$spec_reads/reads$total_reads # % reads of a given species relative to the whole scat
reads$key <- paste0(reads$fieldid,"_",reads$rep,reads$speciessimpleid) # make our key to remerge

# merge back with the larger dataframe
db <- merge(diet3,reads, by="key", all.y = F)

write.csv(db,"ProcessedCoyoteRedFoxStudyArea_DietGeoData_250404.csv", row.names = FALSE)

# d1 now has the number of times 
library(plyr)
names(db)
d1 = ddply(db,
           .(fieldid.x,
             labid,
             rep.x,
             defecator,
             speciessimpleid.x,
             final_scientific,
             species,
             genus,
             family,
             order,
             class,
             phylum,
             final_common,
             genus_group,
             family_group,
             order_group,
             class_group,
             best_match_pid.mean,
             total_reads,
             rep_reads,
             spec_reads,
             p_reads), summarise,
           replicates = ifelse(rep_reads == 0, 0, 1)) 

# filter by our needs to get our final data frame. we can use this dataframe to make plots
draft.dat <- d1[which(d1$p_reads>=0.01),] 
sampcount <- unique(draft.dat$fieldid.x)
write.csv(draft.dat,"ProcessedCoyoteRedFoxStudyArea_DietGeoDataSimplified_250404.csv", row.names = FALSE)
