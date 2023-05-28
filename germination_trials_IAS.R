#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    SEED BANK ANALYSIS IN ICELANDIC RANGELANDS
#             GERMINATION TRIALS
#       Isabel C Barrio (isabel@lbhi.is)
#                July 21, 2022
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the script analyses data on seedling emergence and aboveground vegetation
# from experimental plots in the FENCES experiment
# https://tundraecology.hi.is/research-projects-2/fences-experiment/

#libraries---- 
#packages to be used
library(readxl)      #to import data directly from Excel
library(dplyr)
library(tidyverse)
library(lubridate)   # to work with dates
library(scales)      # to work with dates
library(ggplot2)     # to make plots
library(lme4)        # to build GLMMs
library(lmerTest)    # to get p-values of GLMMs
library(performance) # to check for overdispersion in Poisson GLMMs
library(vegan)       # to calculate multivariate statistics (e.g. similarity indices)
library(gridExtra)


#load datasets----
# contains dataset contains data on the number and identity of germinated seedlings 
#detailed germination data: the file contains the following variables:
    #plot_ID: random number assigned to identify each plot (numbers from 1-48)
    #plot: plot name as in the field experiment
    #species: species ID for seedlings
    #date: date when observation was recorded
    #comments: any comments or remarks
germination <- read_excel("germination_trials_AK.xlsx", sheet= "seedlings") %>%
                  #make sure date is understood as a date
                  mutate(day=as.Date(dmy(date)),
                         #make sure the species names are the same as for aboveground assessments
                         species=case_when(species == "RUMACE" ~ "RUMACETOSELLA",
                                           species == "RUMOSA" ~ "RUMACETOSA",
                                           species == "unidentified" ~ "UNID",
                                           TRUE ~ species))

write_excel_csv(germination, file = "germination_trials.csv", col_names=TRUE)
# germination <- read_csv("germination_trials.csv")

#sp.names: the file contains the following variables:
    #species: code for the species used in germination trials and aboveground assessments
    #sp_name: species name
    #life_form: perennial or annual (all perennials)
sp.names <- read_excel("germination_trials_AK.xlsx", sheet= "species.codes") 
write_excel_csv(sp.names, file = "species_names.csv", col_names=TRUE)

# sp.names <- read_csv("species_names.csv")

#plot.info: the file contains the following variables:
    #plot_ID: random number assigned to identify each plot (numbers from 1-48)
    #plot: plot name as in the field experiment
    #pair: plot pair (using NutNet codes for plots: AH1C, AH1NPK, AH2C…)
    #location	study location (Audkuluheidi or Theistareykir)
    #habitat	habitat type (H: heath or M: melur)
    #treatment	experimental treatment (F: fenced, grazing excluded; C: control, open to grazing)
plot.info <- read_excel("germination_trials_AK.xlsx", sheet= "plot.info") %>% 
                  #add full names for the habitats and treatments
                  #add variable for combined site_habitat
                  mutate(habitat=case_when(habitat == "M" ~ "melur",
                                           habitat == "H" ~ "heath"),
                         treatment=case_when(treatment == "C" ~ "control",
                                             treatment == "F" ~ "fenced"),
                         site.hab = paste(location, habitat, sep="_"),
                         site.hab.ttm = paste(location, habitat, treatment, sep="_")) %>% 
                  arrange(plot)

write_excel_csv(plot.info, file = "plot_info.csv", col_names=TRUE)
# plot.info <- read_csv("plot_info.csv")

#plot codes (point intercept data uses the original coding for the plots, but we used the NutNet code)
setwd("C:/Users/isabel/OneDrive - Menntaský/ISABEL/FENCES/vegetation")
plot.codes <- read_excel("codes_plots.xlsx", sheet= "conversion_codes") %>% 
                mutate(NutNet=case_when(NutNet=="AH2F"~"AH2CF",
                                        TRUE~NutNet),
                       plot=FENCES, plot.NutNet=NutNet) #make sure we have the same column names to merge datasets later

#point intercept data for 2021 (point intercept frames in rows, species in columns)
setwd("C:/Users/isabel/OneDrive - Menntaský/ISABEL/FENCES/vegetation")
pi.data <- read_excel("point_frame_intercept_Fence_Audkuluheidi_Theistareykir_2021.xlsx", sheet= "data_original") %>%  
                   #fix some typing issues
                   mutate(SAXCES=SAXCEP, LOIPRO=LOUPRO, SILUNI=SILMAR,
                          plot=case_when(plot == "TMOC" ~ "TM0C",
                                         TRUE ~ as.character(plot))) %>% 
                   #get the plot names for the NutNet plots
                   left_join(plot.codes, by="plot") %>%
                      select(-plot) %>% mutate(plot=plot.NutNet) %>% #keep only plot names from NutNet
                   #keep only vascular species; check with Inga: SANSP
                   select(plot, site, habitat, ttm, sitehabitat, pointframe, THYARC:POAALP, FESRIC:MINUARTIASP, SILUNI, 
                          LUZSPI:AGRVIN, MINRUB, DRANOR, TRISPI, SAXCES, THAALP, JUNTRI, CARBIG, DRYOCT, BETNAN, VACULI:EQUVAR, 
                          EQUARV, KOBMYO, AGRCAP, SALLAN, CARRUP, AGRSTO:POAPRA, SALHER, TOFPUS, HARHYP, PINVUL, BARALP, 
                          SALARC, CERFON, CERALPSSPGLA, ALCALP, POASP, ANTALP, FESRUB, DESFLE, ANTODO, LOIPRO, RANACR, CALVUL, 
                          GALVER, COEVIR, PARPAL, PLATHYP, VACMYR, ERIBOR, HIERACIUM, OMASUP, CARVAG, SALPHY) %>% 
                   group_by(plot) %>% select(-plot, -site, -habitat, -ttm, -sitehabitat, - pointframe) %>% 
                   summarise_all(sum) %>% 
                   #keep plot ID as row names
                   column_to_rownames("plot")

setwd("C:/Users/isabel/OneDrive - Menntaský/ISABEL/students/Abdubakir")

write_excel_csv(pi.data, file = "point_intercept_data.csv", col_names=TRUE)
# pi.data <- read_csv("point_intercept_data.csv")


#customised functions----
#defining the parameters for our graphs
theme_seedlings <- function(){    #create a new theme function for the style of graphs
  theme_bw()+                  #use a predefined theme as a base
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
        plot.title = element_text(size = 20, vjust = 1, hjust = 0),
        legend.text = element_text(size = 12, face = "italic"),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.9))
}


#GERMINATION DATA----
#remove plants not emerging from seeds
germination.seed <- germination %>% 
                      filter(!species %in% c("EQUVAR", "SALHER", "BISVIV"))
  
##Figure S1----
#number of emerged seedlings per day for the duration of the experiment
#and accumulated number of seedlings
#create a vector for the length of the germination trials (Jun 20 to Sep 20, 2022)
duration.trials <- data.frame(day=seq.Date(as.Date("2022-06-20"), as.Date("2022-09-20"), 1))
seedlings.day <- germination.seed %>% 
                      group_by(day) %>% 
                       summarise(N=n()) %>% 
                       mutate(accum=cumsum(N)) %>% 
                 full_join(duration.trials, by="day") %>% arrange(day) %>% 
                 replace(is.na(.), 0) #replace NAs with zeros

#after one month (Jul 20): 36 seedlings had emerged
(36/45)*100
    #when considering also plants not emerging from seeds:
    (40/50)*100

ggplot(seedlings.day) +
	geom_bar(aes(x = day, y = N), stat="identity", fill = "#4A708B") +
	scale_x_date(breaks=as.Date(c("2022-07-1", "2022-08-1", "2022-09-1")),
	             labels= c("1 Jul", "1 Aug", "1 Sep"),  #an easier way to specify this
	                                                    #would be using: date_labels="%b %d"
	                                                    #but it gives me month names in Icelandic
	             limits=as.Date(c("2022-06-20", "2022-09-20"))) +
  geom_line(aes(x=day, y=cumsum(N)), colour="#4A708B", size=0.7, linetype="dashed") +
	theme_seedlings() +
  geom_vline(xintercept=as.Date("2022-07-20"), color="orange", size = 0.8) +
  labs(x="", y="Number of seedlings") 

##Table 1----
#which species germinated, how many and where
#we report here plants from non-seed propagules as well
seedlings.spp <- germination %>% 
                        left_join(plot.info, by="plot") %>% 
                     #calculate how many seedlings for each species in each habitat
                     group_by(species, habitat) %>% 
                       summarise(N=n()) %>% ungroup() %>% 
                     #separate columns for each habitat
                     pivot_wider(names_from=habitat, values_from=N, values_fill=0) %>% 
                     mutate(total=heath+melur) %>% 
                       arrange(desc(total)) %>% 
                     left_join(sp.names, by="species") %>% 
                      mutate(species=sp_name,
                             nr_seedlings = paste(total, " (", heath, "/", melur, ")", sep="")) %>% 
                      select(species, life_form, nr_seedlings) 
                      
write.csv(seedlings.spp, "table1.csv", row.names=FALSE)

# version of table including data per site/hab/ttm
seedlings.spp <- germination %>% 
                        left_join(plot.info, by="plot") %>% 
                     #calculate how many seedlings for each species in each habitat
                     group_by(species, site.hab.ttm) %>% 
                       summarise(N=n()) %>% ungroup() %>% 
                     #separate columns for each site/hab/ttm
                     pivot_wider(names_from=site.hab.ttm, values_from=N, values_fill=0) %>% 
                     mutate(total = Audkuluheidi_heath_control + Audkuluheidi_heath_fenced +
                                             Audkuluheidi_melur_control + Audkuluheidi_melur_fenced +
                                             Theistareykir_heath_control + Theistareykir_heath_fenced + 
                                             Theistareykir_melur_control + Theistareykir_melur_fenced) %>% 
                       arrange(desc(total)) %>% 
                     left_join(sp.names, by="species") %>% 
                      mutate(species=sp_name) %>% 
                      select(species, total, Audkuluheidi_heath_control, Audkuluheidi_heath_fenced,
                                             Audkuluheidi_melur_control, Audkuluheidi_melur_fenced,
                                             Theistareykir_heath_control, Theistareykir_heath_fenced, 
                                             Theistareykir_melur_control, Theistareykir_melur_fenced) 
                      
write.csv(seedlings.spp, "table1.csv", row.names=FALSE)

#SUMMARY DATA----
#summarize per plot
seedlings <- germination.seed %>% 
                #calculate how many seedlings in each plot
                group_by(plot) %>% 
                  summarise(nr.seedlings=n()) %>% ungroup() %>%         
                full_join(plot.info, by="plot") %>% 
                replace(is.na(.), 0) %>%  #replace NAs with zeros
                select(plot_ID, plot, pair, location, habitat, treatment, site.hab, site.hab.ttm, nr.seedlings) %>% 
                mutate(seedlings.m2=(nr.seedlings/12.5)*10000)

#quick visualisation of the data
ggplot(seedlings, aes(nr.seedlings)) + 
  geom_histogram() + 
  theme_seedlings()
  #distribution of number of seedlings per plot

#A) DIFFERENCES OF THE EFFECT OF FENCING BETWEEN HABITATS AND SITES
#start by plotting; we want to plot the summary data (mean nr seedlings per plot at each 
#site/habitat/treatment combination but we also want to add information on
#species richness belowground to the same graph
spp.richness.B <- germination.seed %>% 
                    #get treatment, plot and site information
                    left_join(plot.info, by="plot") %>% 
                    select(-plot_ID.y) %>% 
                    #calculate species richness per location_habitat_treatment
                    group_by(site.hab, site.hab.ttm, species) %>% tally() %>% 
                    ungroup() %>% group_by(site.hab, site.hab.ttm) %>% 
                    summarize(nr.seedlings=sum(n), sp.rich=length(species))

summary.germination <- seedlings %>%
                       mutate(germ = ifelse(nr.seedlings>0, 1, 0)) %>% 
                       group_by(site.hab, treatment) %>% 
                       summarise(N=n(),
                                 N.germ = sum(germ),
                                 mean=mean(nr.seedlings), 
                                 sd=sd(nr.seedlings), 
                                 se=sd/sqrt(N),
                                 ci=1.96*se) %>% 
                       mutate(site.hab.ttm=paste(site.hab, treatment, sep="_")) %>% 
                       left_join(spp.richness.B, by="site.hab.ttm") %>% 
                       mutate(site.hab=site.hab.x) %>% 
                       select(-site.hab.y, -site.hab.x)

##Figure 1----
#number of emerged seedlings
ggplot(summary.germination, aes(x=site.hab, y=mean, fill=treatment, group=treatment)) +       
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(x=site.hab, ymin=mean-se, ymax=mean+se), 
                  width=.1, position=position_dodge(0.9)) +
    scale_fill_manual(values=c("#4A708B", "#6CA6CD"), labels=c("grazed", "non-grazed")) +
    geom_text(aes(x=site.hab, y=5, label=sp.rich), position=position_dodge(0.9), col="black") +
    scale_y_continuous(limits=c(-0.1,4.5)) +     #expand the limits of y-axis
    scale_x_discrete(labels=c('Auðkúluheiði heath', 'Auðkúluheiði melur', 
                              'Þeistareykir heath', 'Þeistareykir melur')) +
    labs(x="", y="Number of seedlings\n") + 
    theme_seedlings() +                 
    theme(legend.position = "right")

# with richness values
ggplot(summary.germination, aes(x=site.hab, y=mean, fill=treatment, group=treatment)) +       
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(x=site.hab, ymin=mean-se, ymax=mean+se), 
                  width=.1, position=position_dodge(0.9)) +
    scale_fill_manual(values=c("#4A708B", "#6CA6CD")) +
    geom_point(aes(x=site.hab, y=5, size=sp.rich, color=treatment), alpha=0.6, 
               position=position_dodge(0.9)) +
    scale_color_manual(values=c("#4A708B", "#6CA6CD")) +
    scale_size_continuous(range = c(7, 25))+
    geom_text(aes(x=site.hab, y=5, label=sp.rich), position=position_dodge(0.9), col="black") +
    scale_y_continuous(limits=c(-0.1,5.5)) +     #expand the limits of y-axis
    scale_x_discrete(labels=c('Auðkúluheiði heath', 'Auðkúluheiði melur', 
                              'Þeistareykir heath', 'Þeistareykir melur')) +
    labs(title="Number of emerged seedlings", x="", y="Number of seedlings\n") + 
    guides(size="none") +
    theme_seedlings() +                 
    theme(legend.position = "right")

##Figure S2----
#number of plots with emerged seedlings
#(figure not included in the end)
ggplot(summary.germination, aes(x=site.hab, y=N.germ, fill=treatment, group=treatment)) +       
    geom_bar(position="dodge", stat="identity") +
    scale_fill_manual(values=c("#4A708B", "#6CA6CD")) +
    scale_y_continuous(limits=c(0,6)) +     #expand the limits of y-axis
    scale_x_discrete(labels=c('Auðkúluheiði heath', 'Auðkúluheiði melur', 
                              'Þeistareykir heath', 'Þeistareykir melur')) +
    labs(title="Number of plots with emerged seedlings", x="", 
         y="Number of plots\n") + 
    guides(size="none") +
    theme_seedlings() +                 
    theme(legend.position = "right")



##Generalized Linear Mixed Models----
#our response variable is count data, with relatively small values, so we should use
#Generalized Linear Mixed Models (GLMMs) with a Poisson distribution

GLMM=glmer(nr.seedlings~treatment*habitat*location + (1|pair), family = poisson, data =seedlings)
  summary(GLMM) #gives a singularity, but this will be probably solved when we simplify the model

GLMM.1a=glmer(nr.seedlings~treatment*habitat*location + (1|pair), family = poisson, data =seedlings)
GLMM.1b=glmer(nr.seedlings~treatment*habitat + treatment*location + habitat*location+ (1|pair), family = poisson, data =seedlings)
  anova(GLMM.1a, GLMM.1b)         #TTM*HAB*LOCATION (p=0.8626) -- because it is not significant we should drop it 

GLMM.2a=glmer(nr.seedlings~treatment*habitat + treatment*location + habitat*location+ (1|pair), family = poisson, data =seedlings)
GLMM.2b=glmer(nr.seedlings~treatment+habitat + treatment*location + habitat*location+ (1|pair), family = poisson, data =seedlings)
  anova(GLMM.2a, GLMM.2b)         #TTM*HAB (p=0.508) -- we drop first the least significant interaction
GLMM.2c=glmer(nr.seedlings~treatment*habitat + treatment+location + habitat*location+ (1|pair), family = poisson, data =seedlings)
  anova(GLMM.2a, GLMM.2c)         #TTM*LOCATION (p=0.0085) 
GLMM.2d=glmer(nr.seedlings~treatment*habitat + treatment*location + habitat+location+ (1|pair), family = poisson, data =seedlings)
  anova(GLMM.2a, GLMM.2d)         #HAB*LOCATION (p=0.4993) 

GLMM.3a=glmer(nr.seedlings~treatment*location + habitat*location+ (1|pair), family = poisson, data =seedlings)
GLMM.3b=glmer(nr.seedlings~treatment+location + habitat*location+ (1|pair), family = poisson, data =seedlings)
  anova(GLMM.3a, GLMM.3b)         #TTM*HABITAT (p=0.01017)
GLMM.3c=glmer(nr.seedlings~treatment*location + habitat+location+ (1|pair), family = poisson, data =seedlings)
  anova(GLMM.3a, GLMM.3c)         #HABITAT*LOCATION (p=0.5985) -- we drop first the least significant interaction

GLMM.4a=glmer(nr.seedlings~treatment*location + habitat+ (1|pair), family = poisson, data =seedlings)
GLMM.4b=glmer(nr.seedlings~treatment+location + habitat+ (1|pair), family = poisson, data =seedlings)
  anova(GLMM.4a, GLMM.4b)         #TREATMENT*SITE (p=0.02758) 

#our final model:
GLMM.final=glmer(nr.seedlings~treatment*location + habitat + (1|pair), family = poisson, data =seedlings)
  summary(GLMM.final)         
plot(GLMM.final) #no funny patterns in residuals --> the assumptions of the model are fine
check_overdispersion(GLMM.final) #no overdispersion, so using Poisson is ok

audk <- seedlings %>% filter(location == "Audkuluheidi")
  GLMM.audk=glmer(nr.seedlings~treatment + habitat + (1|pair), family = poisson, data =audk)
    summary(GLMM.audk)
theis <- seedlings %>% filter(location == "Theistareykir")
  GLMM.theis=glmer(nr.seedlings~treatment + habitat + (1|pair), family = poisson, data =theis)
    summary(GLMM.theis)
    


#COMPARISON ABOVE AND BELOWGROUND----
#community data for seedbanks
plot.names <- seedlings %>% select(plot)
seed.banks <- germination.seed %>% 
              #calculate number of seedlings of each species in each plot
              group_by(plot, species) %>% tally() %>%  
              pivot_wider(names_from = species, values_from = n, values_fill = 0) %>% 
              #make sure that we have all plots
              full_join(plot.names, by="plot") %>% 
              arrange(plot) %>% 
              replace(is.na(.), 0) %>% 
              #make sure plot names of soil seed banks are identified
              mutate(plot=paste(plot,"_sb",sep="")) %>% 
              #keep plot ID as row names
              column_to_rownames("plot")

#combine dataset for seed banks and aboveground vegetation
above.belowgr <- pi.data %>% bind_rows(seed.banks) %>% 
                   replace(is.na(.), 0)
above.belowgr.presabs <- above.belowgr #create a presence/absence matrix, for the species present in each plot
    above.belowgr.presabs[above.belowgr.presabs > 0] <- 1 #convert to presence/absence
# remove rows with all zeros
above.belowgr.presabs <- above.belowgr.presabs %>% filter_all(any_vars(. != 0))
    
#calculate similarity between above and belowground species
sorensen <- vegdist(above.belowgr.presabs, binary=TRUE)
    
#extract values for above-belowground comparisons
    mSorensen <- as.matrix(sorensen)
      row_names <- grep("_sb", rownames(mSorensen), value=TRUE)
      col_names <- match(sub("_sb","",row_names),colnames(mSorensen))
    similarity <- data.frame(plot=colnames(mSorensen)[col_names],
                         similarity=(1-diag(mSorensen[row_names,col_names]))*100)

#add similarity info to seedlings dataset
seedlings <- seedlings %>% 
                left_join(similarity, by="plot") 

seedlings %>% filter(similarity == 0)
(7/20)*100 #35% of observations are zero!

max(seedlings$similarity, na.rm=T)
mean(seedlings$similarity, na.rm=T)
sd(seedlings$similarity, na.rm=T)/sqrt(20)

##Figure 2----
n.germ <- summary.germination$N.germ
ggplot(seedlings, aes(x=site.hab, y=similarity, fill=treatment)) +       
    geom_boxplot() +
    scale_fill_manual(values=c("#4A708B", "#6CA6CD"), labels=c("grazed", "non-grazed")) +
    scale_colour_manual(values=c("#4A708B", "#6CA6CD")) +
    scale_y_continuous(limits=c(-1,15)) +     #expand the limits of y-axis
    scale_x_discrete(labels=c('Auðkúluheiði heath', 'Auðkúluheiði melur', 
                              'Þeistareykir heath', 'Þeistareykir melur')) + 
        annotate("text", x= 0.8, y=-1, label="3") +
        annotate("text", x= 1.2, y=-1, label="2") +
        annotate("text", x= 1.8, y=-1, label="1") +
        annotate("text", x= 2.2, y=-1, label="3") +
        annotate("text", x= 2.8, y=-1, label="5") +
        annotate("text", x= 3.2, y=-1, label="3") +
        annotate("text", x= 3.8, y=-1, label="2") +
        annotate("text", x= 4.2, y=-1, label="1") +
    labs(x="", y="Similarity index (%)\n") + 
    theme_seedlings() +                 
    theme(legend.position = "right")

  seedlings %>% group_by(site.hab.ttm) %>% 
                       summarise(N=n())   

  
#build NMDS plots separately for above and belowground
### aboveground communities----
#non-metric multidimensional scaling (NMDS)
Z.AG <- metaMDS(pi.data)
  plot(Z.UF, display="species", type="n", las=1, main= "Above- and belowground communities")
    orditorp(Z.UF, display="sites", col="Grey", air=1)
    orditorp(Z.UF, display="species",col="DarkRed",air=1)
  #site.col <- c("#66CC00","#99CCFF","#FF6666","#006600","#FFCC99")
  #  ordiellipse(Z.UF, herbivore.diversity$site, draw="polygon", border=site.col, col=site.col, alpha=80, lty=NULL)
  #stressplot(Z.UF)

#plot with ggplot
#extract the scores for each site and species
data.scores.AG <- as.data.frame(scores(Z.AG)$sites) %>%   #use scores function to extract site scores 
                    mutate(habitat = plot.info$habitat, #create a column of patch names
                           site = plot.info$location, 
                           treatment = plot.info$treatment,
                           hab_ttm = paste0(habitat, " ", ifelse(treatment == "control", "grazed", "non-grazed"))) 
plot.AG <- ggplot() + 
  geom_point(data = data.scores.AG, aes(x = NMDS1, y = NMDS2, col = site, shape = hab_ttm), 
             alpha = 0.8, size = 2, show.legend = TRUE) + 
    scale_color_brewer(palette="Dark2") +
  annotate("text", -0.9, -0.6, label = "stress = 0.13") +
  scale_shape_manual("habitat and treatment", values = c(15, 0, 17, 2)) +
  scale_colour_manual(values=c("#4A708B", "#6CA6CD"), labels = c("Auðkúluheiði", "Þeistareykir")) + 
  coord_equal() +
  labs(title = "a. Aboveground vegetation") +
  theme_seedlings() +
  theme(legend.position = "right")

plot.AG

#species richness aboveground (to include in Table 1)
#overall
pi.data %>%
    mutate(plot = rownames(pi.data)) %>% 
    #get treatment, plot and site information
   left_join(plot.info, by="plot") %>% 
   pivot_longer(cols = THYARC:SALPHY, names_to = "species") %>% 
   group_by(location, species) %>% 
     summarize(total = sum(value)) %>% 
   mutate(value = replace(total, total>0, 1)) %>% #change to presence/absence
   ungroup() %>% 
     summarize(sp.rich.A = sum(value))

#by site, treatment and habitat
pi.data %>%    
   mutate(plot = rownames(pi.data)) %>% 
   #get treatment, plot and site information
   left_join(plot.info, by="plot") %>% 
   pivot_longer(cols = THYARC:SALPHY, names_to = "species") %>% 
   group_by(site.hab.ttm, species) %>% 
    summarize(total = sum(value)) %>% 
   mutate(value = replace(total, total>0, 1)) %>% #change to presence/absence
   group_by(site.hab.ttm) %>% 
    summarize(sp.rich.A = sum(value))

### belowground communities----
# non-metric multidimensional scaling (NMDS)
# remove rows with all zeros
plot.names <- seedlings %>% select(plot)
seed.banks <- germination.seed %>% 
              #calculate number of seedlings of each species in each plot
              group_by(plot, species) %>% tally() %>%  
              pivot_wider(names_from = species, values_from = n, values_fill = 0) %>% 
              #make sure that we have all plots
              full_join(plot.names, by="plot") %>% 
              arrange(plot) %>% 
              replace(is.na(.), 0) %>% 
              #keep plot ID as row names
              column_to_rownames("plot")
seed.banks.pres <- seed.banks %>% filter_all(any_vars(. != 0))

Z.BG <- metaMDS(seed.banks.pres, trymax = 100)

#extract the scores for each sample
data.scores.BG <- as.data.frame(scores(Z.BG)$sites) %>%   #use scores function to extract site scores 
                    rownames_to_column("plot") %>% 
                    left_join(plot.info, by = "plot") %>% 
                    mutate(hab_ttm = paste0(habitat, " ", ifelse(treatment == "control", "grazed", "non-grazed"))) 

plot.BG <- ggplot() + 
  geom_point(data = data.scores.BG, aes(x = NMDS1, y = NMDS2, col = location, shape = hab_ttm), 
             alpha = 0.8, size = 2, show.legend = TRUE) + 
    scale_color_brewer(palette="Dark2") +
  annotate("text", -8, -12, label = "stress = 0") +
  scale_shape_manual("habitat and treatment", values = c(15, 0, 17, 2)) +
  scale_colour_manual(values=c("#4A708B", "#6CA6CD"), labels = c("Auðkúluheiði", "Þeistareykir")) + 
  coord_equal() +
  labs(title = "b. Belowground vegetation") +
  theme_seedlings() +
  theme(legend.position = "right")

plot.BG

grid.arrange(plot.AG, plot.BG, ncol = 1, nrow = 2)


#species richness belowground (to include in Table 1)
#overall
germination.seed %>% 
                    #get treatment, plot and site information
                    left_join(plot.info, by="plot") %>% 
                    select(-plot_ID.y) %>% 
                    #calculate species richness
                    group_by(species) %>% tally() %>% 
                    ungroup() %>% summarize(sp.rich=length(species))
#by site, treatment and habitat
germination.seed %>% 
                    #get treatment, plot and site information
                    left_join(plot.info, by="plot") %>% 
                    select(-plot_ID.y) %>% 
                    #calculate species richness per location_habitat_treatment
                    group_by(site.hab, site.hab.ttm, species) %>% tally() %>% 
                    ungroup() %>% group_by(site.hab, site.hab.ttm) %>% 
                    summarize(nr.seedlings=sum(n), sp.rich=length(species))

