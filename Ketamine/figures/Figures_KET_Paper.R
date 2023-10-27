# Making figures for the Ketamine paper
# Polished version
# 1. Dependencies & Set up ####
# Data import and figures
library('dplyr')
library('ggplot2')
library('data.table')
library('ggh4x')
library(readxl)
library(ggrepel)

# Statistics
library(tidyverse)
library(ggpubr)
library(rstatix)

# Set path and file_name
paper_path <- "C:/Users/Paolo/dox/PhD/00_Publications/Papers/Mine/2023/Pilot_KET/NEJM_RapidReport/results"
data_LZC_allcon_path <- 'pilotKET_merged_R_V3_all.txt'
data_PerCon_path <- 'pilotKET_merged_R_V3_PerConc.txt'
data_alpha_path <-  'pilotKET_merged_R_alpha_v4.txt'
data_ECI_path <-  'ECI_table_R_bas_full.csv'
  
# 2. Figures on the Complexity ####
setwd("C:/Users/Paolo/OneDrive - Universite de Liege/Bureau/OldComputer/D/Complexit_doc/Ketamine/results")

# All trials
data_LZC_allcon <- read.csv(data_LZC_allcon_path, header = TRUE, sep = ',')
data_LZC_allcon[c('Subject', 'Electrode', 'Concentration', 'Session', 'Condition')] <- 
  lapply(data_LZC_allcon[c('Subject', 'Electrode', 'Concentration', 'Session', 'Condition')],as.factor)
levels(data_LZC_allcon$Subject) <- c('MCS-','UWS', 'MCS+')
data_LZC_allcon$Subject <- factor(data_LZC_allcon$Subject, levels =  c('UWS','MCS-','MCS+'))

data_LZC_allcon_summarize <-data_LZC_allcon %>%
  group_by(Subject, Condition) %>%
  summarise(mean_LZC = mean(LZC), sd_LZC = sd(LZC))

# Statistics
res.aov_LZC_all <- data_LZC_allcon_summarize %>% as.data.frame() %>% anova_test(dv = mean_LZC, wid = Subject, within = c(Condition))
get_anova_table(res.aov_LZC_all)
t.test(x = pull(data_LZC_allcon_summarize[data_LZC_allcon_summarize['Condition'] == 'Ketamine',3]),
       y = pull(data_LZC_allcon_summarize[data_LZC_allcon_summarize['Condition'] == 'Placebo',3]),
             paired = TRUE)

#Graph with SE
data_LZC_allcon_summarize %>% ggplot(aes(x=factor(Condition, level=c('Placebo', 'Ketamine')), y=mean_LZC, color=Subject,
                                         group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line(aes(group =Subject)) +
  geom_errorbar(aes(ymin=mean_LZC-sd_LZC/sqrt(127), ymax=mean_LZC+sd_LZC/sqrt(127)), width=.2) + 
  xlab('') + ylab('Lempel-Ziv complexity') + theme_bw() + theme(text = element_text(size = 20))

ggsave(paste0(today(),"_KET_LZC_alltrials_se.tiff"), device="tiff", width=16, height=10, units="cm", path=paper_path, dpi=300)

# Channel-wise
LZC_graph_chanwise <- data_LZC_allcon %>% 
  ggplot(aes(x=factor(Condition, level=c('Placebo', 'Ketamine')), y=LZC, 
             color=Subject, group =interaction(Condition, Electrode)))+
  geom_point(size=3)+ geom_line(aes(group = Electrode)) +
  ylab('Lempel-Ziv complexity') +xlab('')+ theme_bw() + 
  theme(text = element_text(size =10))+ 
  scale_y_continuous(limits = c(0.15, 0.61))
LZC_graph_chanwise + facet_grid( .~ Subject )
#ggsave(paste0(today(),"_KET_LZC_alltrials_chanwise.tiff"), device="tiff", width=16, height=10, units="cm", path=paper_path, dpi=300)

# 2.1.0 Whole brain LZC #
data_PerCon <- read.csv(data_PerCon_path, header = TRUE, sep = ',')
data_PerCon[c('Electrode', 'Concentration', 'Session', 'Condition')] <- 
  lapply(data_PerCon[c('Electrode', 'Concentration', 'Session', 'Condition')],as.factor)

# 2.1.1 Whole brain LZC w/ standard mean error #
data_PerCon_summarize <-data_PerCon %>%
  group_by(Subject, Condition, Concentration) %>%
  summarise(mean_LZC = mean(LZC), sd_LZC = sd(LZC))

# Statistics
res.aov_LZC <- data_PerCon_summarize %>% as.data.frame() %>% anova_test(dv = mean_LZC, wid = Subject, within = c(Concentration, Condition))
get_anova_table(res.aov_LZC)

# Standard error - all together
data_PerCon_summarize$Subject <- as.factor(data_PerCon_summarize$Subject)
levels(data_PerCon_summarize$Subject) <- c('MCS-','UWS','MCS+')
data_PerCon_summarize$Subject <- factor(data_PerCon_summarize$Subject, levels =  c('UWS','MCS-','MCS+'))

LZC_graph_summary <- data_PerCon_summarize %>% 
  ggplot(aes(x=Concentration, y=mean_LZC, shape=Condition, color=Condition,
             group =interaction(Condition, Subject)))+
  geom_point(size=3)+ geom_line() +
  geom_errorbar(aes(ymin=mean_LZC-sd_LZC/sqrt(127), ymax=mean_LZC+sd_LZC/sqrt(127)), width=.2) + 
  ylab('Lempel-Ziv complexity') +xlab('Concentration (µg/ml)')+ theme_bw() + 
  theme(text = element_text(size =10))+ 
  scale_y_continuous(limits = c(0.15, 0.61)) + geom_vline(xintercept = 1.5, linetype="dotted", 
                                                        color = "black", size=2)
LZC_graph_summary + facet_grid( .~ Subject )
ggsave(paste0(today(),"_LZC_PerCon.tiff"), device="tiff", width=16, height=10, units="cm", path=paper_path, dpi=300)



# 2.1 Alpha centrality ########
data_alpha <- read.csv(data_alpha_path, header = TRUE, sep = ',')

colnames(data_alpha)[3] <- "Concentration" # Not dose, but concentration
data_alpha[c('Concentration')] <- lapply(data_alpha[c('Concentration')],as.factor)
levels(data_alpha$Concentration) <- c("0", "0.15-0.30", "0.45-0.60", "0.75")

data_alpha$Subject <- as.factor(data_alpha$Subject)
levels(data_alpha$Subject) <- c('MCS-','UWS', 'MCS+')
data_alpha$Subject <- factor(data_alpha$Subject, levels =  c('UWS','MCS-','MCS+'))

alpha_graph <- data_alpha %>% ggplot(aes(x=Concentration, y=Alpha, shape=Condition, color=Condition,
                          group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line() + 
  ylab('Alpha Centrality') +xlab('Concentration (µg/ml)')+ theme_bw() +
  theme(text = element_text(size = 25), axis.text.x = element_text(size = 18))+ 
  geom_vline(xintercept = 1.5, linetype="dotted", color = "black", size=2)
alpha_graph + facet_grid(. ~ Subject)

ggsave(paste0(today(),"_KET_alpha_withbas.tiff"), device="tiff", width=37, height=20, units="cm", path=paper_path, dpi=300)

res.aov_alpha_centrality <- data_alpha %>% #as.data.frame() %>% 
  anova_test(dv = Alpha, wid = Subject, within = c(Concentration, Condition))
get_anova_table(res.aov_alpha_centrality)


# 2.2 ECI ####
data_ECI <- read.csv(data_ECI_path, header = TRUE, sep = ',')
data_ECI$Session <- as.factor(data_ECI$Session )
data_ECI$Subject <- as.factor(data_ECI$Subject )
levels(data_ECI$Subject) <- c('MCS-','UWS', 'MCS+')
data_ECI$Subject <- factor(data_ECI$Subject, levels =  c('UWS','MCS-','MCS+'))

ECI_aro_graph <- data_ECI %>% filter(Concentration != 'All') %>% filter(Concentration != '0')%>% 
  ggplot(aes(x=Concentration, y=ECI_Aro, color=Condition,
                                         group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line() + 
  ylab('ECI Arousal') +xlab('Concentration (µg/ml)')+ theme_bw() +
  theme(text = element_text(size = 25), axis.text.x = element_text(size = 18)) +
  ylim(c(0,1)) +geom_hline(yintercept = 0.5, linetype="dotted",color = "red", size=2)+
  geom_vline(xintercept = 1.5, linetype="dotted", color = "black", size=2)
ECI_aro_graph + facet_grid(. ~ Subject)
ggsave(paste0(today(),"_KET_ECIAro.tiff"), device="tiff", width=37, height=20, units="cm", path=paper_path, dpi=300)

ECI_awa_graph <- data_ECI %>% filter(Concentration != 'All') %>% #filter(Concentration != '0')%>% 
  ggplot(aes(x=Concentration, y=ECI_Awa, color=Condition,
                                         group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line() + 
  ylab('ECI Awareness') +xlab('Concentration (µg/ml)')+ theme_bw() +
  theme(text = element_text(size = 25), axis.text.x = element_text(size = 18)) +
  ylim(c(0,1)) +  geom_hline(yintercept = 0.5, linetype="dotted", color = "red",  size=2) + 
  geom_vline(xintercept = 1.5, linetype="dotted", color = "black", size=2)
ECI_awa_graph + facet_grid( .~ Subject )

ggsave(paste0(today(),"_KET_ECIAwa.tiff"), device="tiff", width=37, height=20, units="cm", path=paper_path, dpi=300)

# Full###
# Arousal ##
data_ECI %>% filter(Concentration == "All") %>% ggplot(aes(x=factor(Condition, level=c('Placebo', 'Ketamine')), y=ECI_Aro,
                                                           color=Subject, group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line(aes(group =Subject)) +
  xlab('') + ylab('ECI Arousal') + theme_bw() + theme(text = element_text(size = 20))+
  ylim(c(0,1)) +  geom_hline(yintercept = 0.5, linetype="dotted", color = "red",  size=2)

ggsave(paste0(today(),"_KET_ECIAro_all.tiff"), device="tiff", width=16, height=10, units="cm", path=paper_path, dpi=300)

# Awa ##
data_ECI %>% filter(Concentration == "All") %>% ggplot(aes(x=factor(Condition, level=c('Placebo', 'Ketamine')), y=ECI_Awa,
                                                           color=Subject, group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line(aes(group =Subject)) +
  xlab('') + ylab('ECI Awareness') + theme_bw() + theme(text = element_text(size = 20))+
  ylim(c(0,1)) +  geom_hline(yintercept = 0.5, linetype="dotted", color = "red",  size=2)
  
ggsave(paste0(today(),"_KET_ECIAwa_all.tiff"), device="tiff", width=16, height=10, units="cm", path=paper_path, dpi=300)


# 2. Figures on the SECONDs ####
setwd("C:/Users/Paolo/dox/PhD/01_Psychedelics/Ketamine/Pilot/results/Tables_dataset")

data_seconds <- read_excel("20230220_SECONDs_R.xlsx", 
                           col_types = c("text", "text", "text", "numeric", "text"))
data_seconds[c('Condition', 'Time', 'Diagnosis')] <- 
  lapply(data_seconds[c('Condition', 'Time', 'Diagnosis')],as.factor)
data_seconds$Time <- ordered(data_seconds$Time, levels = c("0", "30", "60", "90", "120")) 
levels(data_seconds$Time) <- c("0", "30", "60", "90", "210")

res.aov_SECONDs <- anova_test(data = data_seconds, dv = Index, wid = Subject, within = c(Time,Condition))
get_anova_table(res.aov_SECONDs)


 
#Change for visual representation
temp1 <- data_seconds %>% filter(Condition == "Placebo") %>% mutate(Index = Index -0.1)
temp2 <- data_seconds %>% filter(Condition == "Ketamine") %>% mutate(Index = Index +0.1)
data_seconds <- rbind(temp1, temp2)

seconds_graph <- data_seconds %>% ggplot(aes(x=Time, y=Index, shape=Condition, color=Condition,
                                         group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line() + 
  ylab('SECONDs score') +xlab('Time')+ theme_bw() +
  theme(text = element_text(size = 25), axis.text.x = element_text(size = 18))
seconds_graph + facet_grid(. ~ Subject)

seconds_graph_v2 <- data_seconds %>% ggplot(aes(x=Time, y=Index, shape=Diagnosis, color=Condition,
                                             group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line() + 
  ylab('SECONDs score') +xlab('Time')+ theme_bw() +
  theme(text = element_text(size = 25), axis.text.x = element_text(size = 18)) 
seconds_graph_v2 + facet_grid(. ~ Subject)
ggsave("KET_seconds_3diagn.tiff", device="tiff", width=30, height=20, units="cm", path=paper_path, dpi=300)

data_seconds_2diag <- data_seconds
levels(data_seconds_2diag$Diagnosis) <- c('MCS', 'MCS', 'UWS')

# Per sub
subject <- 'KET03'
seconds_graph_v3 <- data_seconds_2diag %>% filter(Subject == subject) %>%
  ggplot(aes(x=Time, y=Index, shape=Diagnosis, color=Condition,
                                                group =interaction(Condition, Subject)))+
  geom_point(size=15)+ geom_line() + scale_shape_manual(values=c(15,18)) +
  ylab('SECONDs score') +xlab('Time')+ theme_bw() +
  theme(text = element_text(size = 40))
seconds_graph_v3 + facet_grid(. ~ Subject) +  
  scale_y_continuous(breaks = seq(0,8,1), limits = c(-0.1, 8.1)) + 
  guides(color = FALSE)
ggsave(paste("KET_SECONDs_",".tiff",sep=subject), device="tiff", width=50, height=20, units="cm", path=paper_path, dpi=300)

# All together
data_seconds_2diag$Subject <- as.factor(data_seconds_2diag$Subject)
levels(data_seconds_2diag$Subject) <- c('MCS-','UWS','MCS+')
data_seconds_2diag$Subject <- factor(data_seconds_2diag$Subject, levels =  c('UWS','MCS-','MCS+'))

seconds_graph_v3 <- data_seconds_2diag %>% 
  ggplot(aes(x=Time, y=Index, shape=Diagnosis, color=Condition,
             group =interaction(Condition, Subject)))+
  geom_point(size=7)+ geom_line() + scale_shape_manual(values=c(15,18)) +
  ylab('SECONDs score') +xlab('Time (min)')+ theme_bw() +
  theme(text = element_text(size = 25))

seconds_graph_v3 + facet_grid(. ~ Subject) +  
  scale_y_continuous(breaks = seq(0,8,1), limits = c(-0.1, 8.1)) + geom_vline(xintercept = 4.5, linetype="dotted", 
                                                                            color = "black", size=2)

ggsave(paste0(today(),"KET_seconds_2diagn.tiff"), device="tiff", width=30, height=20, units="cm", path=paper_path, dpi=300)

# 3. Figures on the MAS ####
setwd("C:/Users/Paolo/dox/PhD/01_Psychedelics/Ketamine/Pilot/results/Tables_dataset")

data_MAS <- read_excel("20230220_Database_KET_MAS.xlsx")
data_MAS[c('Condition', 'Time', 'Side', 'Limb')] <- 
  lapply(data_MAS[c('Condition', 'Time', 'Side', 'Limb')],as.factor)

# Statistics
res.aov_MAS <- data_MAS %>% 
  filter(Subject != 'KET01')%>% anova_test(dv = MAS, wid = Subject, within = c(Time,Condition, Side, Limb))
get_anova_table(res.aov_MAS)

#Change for visual representation
data_MAS$MAS <- data_MAS$MAS + seq(-0.1, 0.1, length=8)

data_MAS$Subject <- as.factor(data_MAS$Subject)
levels(data_MAS$Subject) <- c('MCS-','UWS','MCS+')


### Figure with "0" MAS - Supplementary?
MAS_graph <- data_MAS %>% 
  filter(Subject != 'MCS-')%>% 
  ggplot(aes(x=Time, y=MAS, shape=Condition, color=Condition, label=Limb,
                                             group =interaction(Condition, Limb, Side)))+
  geom_point(size=5)+ geom_line() + 
  ylab('Modified Ashworth Scale') +xlab('Time (min)')+ theme_bw() +
  theme(text = element_text(size = 25), 
        axis.text.x = element_text(size = 18)) +
  geom_text_repel(size=4)
  
MAS_graph + theme(legend.position='none') + facet_nested(Condition ~ Subject + Side )
ggsave(paste0(today(),"_KET_MASS.tiff"), device="tiff", width=30, height=20, units="cm", path=paper_path, dpi=300)

### Figure with no "0" MAS - PAPER
data_MAS_mod <- data_MAS %>% filter(Subject!='MCS-')
data_MAS_mod <- as.data.frame(data_MAS_mod)
data_MAS_mod <- data_MAS_mod[-c(2, 10, 8, 16, 22, 30, 18, 26, 24, 32),]

MAS_graph_mod <- data_MAS_mod %>% 
  ggplot(aes(x=Time, y=MAS, shape=Condition, color=Condition, label=Limb,
             group =interaction(Condition, Limb, Side)))+
  geom_point(size=5)+ geom_line() + 
  ylab('Modified Ashworth Scale') +xlab('Time (min)')+ theme_bw() +
  theme(text = element_text(size = 25), 
        axis.text.x = element_text(size = 18)) +
  geom_text_repel(size=4)

MAS_graph_mod + theme(legend.position='none') + facet_nested(Condition ~ Subject + Side)
ggsave(paste0(today(),"_KET_MASS_mod.tiff"), device="tiff", width=30, height=20, units="cm", path=paper_path, dpi=300)
