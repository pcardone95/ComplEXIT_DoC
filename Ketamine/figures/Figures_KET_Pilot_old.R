# Making figures for the Ketamine paper

# 1. Dependencies & Set up ####
# Data import and figures
library('dplyr')
library('ggplot2')
library('data.table')
library('ggh4x')
library(readxl)
library(ggrepel)
paper_path <- "C:/Users/Paolo/dox/PhD/00_Publications/Papers/Mine/2023/Pilot_KET/NEJM_RapidReport/results"

# Statistics
library(tidyverse)
library(ggpubr)
library(rstatix)

# 2. Figures on the Complexity ####
setwd("C:/Users/Paolo/OneDrive - Universite de Liege/Bureau/OldComputer/D/Complexit_doc/Ketamine/results")

# 2.1.0 Whole brain LZC #
# data <- read.csv('pilotKET_merged_R_complete.txt', header = TRUE, sep = ',')
# data[c('Electrode', 'Dose')] <- lapply(data[c('Electrode', 'Dose')],as.factor)
data <- read.csv('LZC_KET_exp_complete.txt', header = TRUE, sep = ',')
data[c('Electrode', 'Concentration', 'Session', 'Condition')] <- 
  lapply(data[c('Electrode', 'Concentration', 'Session', 'Condition')],as.factor)
data[c('Electrode', 'Session', 'Condition')] <- 
  lapply(data[c('Electrode', 'Session', 'Condition')],as.factor)

LZC_graph <- data %>% ggplot(aes(x=Concentration, y=LZC, color=Condition, fill=Subject))+#, color=Condition, group_by(Subject))) +
  geom_line(aes(group=Subject), alpha=0.2) +
  geom_pointrange(stat = "summary_bin") +
  geom_line(stat = "summary_bin", fun.y = "mean") + # group = Subject
  ylab('Lempel-Ziv complexity') + theme_bw()

LZC_graph + facet_grid(. ~ Subject)

# Distribution data %>% ggplot(aes(x=Concentration, y=LZC, color=Condition, fill=Subject))+
#geom_point()

# 2.1.1 Whole brain LZC w/ standard mean error #
data_summarize <-data %>%
  group_by(Subject, Condition, Concentration) %>%
  summarise(mean_LZC = mean(LZC), sd_LZC = sd(LZC))

# Statistics
res.aov_LZC <- data_summarize %>% as.data.frame() %>% anova_test(dv = mean_LZC, wid = Subject, within = c(Concentration, Condition))
get_anova_table(res.aov_LZC)

# Standard error - Per sub
subject = 'KET01'
data_summarize %>% filter(Subject== subject) %>% 
  ggplot(aes(x=Concentration, y=mean_LZC, shape=Condition, color=Condition,
                                                            group =interaction(Condition, Subject)))+
  geom_point(size=10)+ geom_line() +
  geom_errorbar(aes(ymin=mean_LZC-sd_LZC/sqrt(127), ymax=mean_LZC+sd_LZC/sqrt(127)), width=.2) + 
  ylab('Lempel-Ziv complexity') +xlab('Concentration (µg/ml)')+ theme_bw() + 
  theme(text = element_text(size =45))+ 
  guides(shape = FALSE, color = FALSE) + scale_y_continuous(limits = c(0.15, 0.61))


# ggsave("KET_LZC_complete_v2.tiff", device="tiff", width=35, height=40, units="cm", path=paper_path, dpi=300)
ggsave(paste("KET_LZC_",".tiff",sep=subject), device="tiff", width=40, height=20, units="cm", path=paper_path, dpi=300)

# Standard error - all together
data_summarize$Subject <- as.factor(data_summarize$Subject)
levels(data_summarize$Subject) <- c('MCS-','UWS','MCS+')
LZC_graph_summary <- data_summarize %>% 
  ggplot(aes(x=Concentration, y=mean_LZC, shape=Condition, color=Condition,
             group =interaction(Condition, Subject)))+
  geom_point(size=3)+ geom_line() +
  geom_errorbar(aes(ymin=mean_LZC-sd_LZC/sqrt(127), ymax=mean_LZC+sd_LZC/sqrt(127)), width=.2) + 
  ylab('Lempel-Ziv complexity') +xlab('Concentration (µg/ml)')+ theme_bw() + 
  theme(text = element_text(size =10))+ 
  scale_y_continuous(limits = c(0.15, 0.61))
LZC_graph_summary + facet_grid( .~ Subject )
ggsave("KET_LZC_all_V2.tiff", device="tiff", width=16, height=10, units="cm", path=paper_path, dpi=300)

# Standard deviation
data_summarize %>% ggplot(aes(x=Concentration, y=mean_LZC, shape=Condition, color=Subject,
                              group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line() +
  geom_errorbar(aes(ymin=mean_LZC-sd_LZC, ymax=mean_LZC+sd_LZC), width=.2) + 
  ylab('LZC') + theme_bw() + theme(text = element_text(size = 20))

# Version 2 - dated 30/08/2023
data <- read.csv('pilotKET_merged_R_V2_all.txt', header = TRUE, sep = ',')
data[c('Subject', 'Electrode', 'Concentration', 'Session', 'Condition')] <- 
  lapply(data[c('Subject', 'Electrode', 'Concentration', 'Session', 'Condition')],as.factor)
levels(data$Subject) <- c('MCS-','UWS', 'MCS+')

data_summarize <-data %>%
  group_by(Subject, Condition) %>%
  summarise(mean_LZC = mean(LZC), sd_LZC = sd(LZC))
#SD or SE
data_summarize %>% ggplot(aes(x=factor(Condition, level=c('Placebo', 'Ketamine')), y=mean_LZC, color=Subject,
                              group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line(aes(group =Subject)) +
  geom_errorbar(aes(ymin=mean_LZC-sd_LZC/sqrt(127), ymax=mean_LZC+sd_LZC/sqrt(127)), width=.2) + 
  #geom_errorbar(aes(ymin=mean_LZC-sd_LZC, ymax=mean_LZC+sd_LZC), width=.2) + 
  xlab('') + ylab('Lempel-Ziv complexity') + theme_bw() + theme(text = element_text(size = 20))

ggsave("KET_LZC_alltrials_August2023_se.tiff", device="tiff", width=16, height=10, units="cm", path=paper_path, dpi=300)

# Channel-wise
LZC_graph <- data %>% 
  ggplot(aes(x=factor(Condition, level=c('Placebo', 'Ketamine')), y=LZC, 
             color=Subject, group =interaction(Condition, Electrode)))+
  geom_point(size=3)+ geom_line(aes(group = Electrode)) +
  ylab('Lempel-Ziv complexity') +xlab('')+ theme_bw() + 
  theme(text = element_text(size =10))+ 
  scale_y_continuous(limits = c(0.15, 0.61))
LZC_graph + facet_grid( .~ Subject )
ggsave("KET_LZC_alltrials_August2023_channelwise.tiff", device="tiff", width=16, height=10, units="cm", path=paper_path, dpi=300)


# 2.1 Alpha centrality ########
data_alpha <- read.csv('pilotKET_merged_R_alpha_v2.txt', header = TRUE, sep = ',')

colnames(data_alpha)[3] <- "Concentration" # Not dose, but concentration
data_alpha[c('Concentration')] <- lapply(data_alpha[c('Concentration')],as.factor)
#levels(data_alpha$Concentration) <- c(".15-.30", ".45-.60", ".75")

data_alpha$Subject <- as.factor(data_alpha$Subject)
levels(data_alpha$Subject) <- c('MCS-','UWS', 'MCS+')

alpha_graph <- data_alpha %>% ggplot(aes(x=Concentration, y=Alpha_centrality, shape=Condition, color=Condition,
                          group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line() + 
  ylab('Alpha Centrality') +xlab('Concentration (µg/ml)')+ theme_bw() +
  theme(text = element_text(size = 25), axis.text.x = element_text(size = 18))
alpha_graph + facet_grid(. ~ Subject)

ggsave("KET_alpha_v3.tiff", device="tiff", width=37, height=20, units="cm", path=paper_path, dpi=300)

res.aov_alpha_centrality <- data_alpha %>% #as.data.frame() %>% 
  anova_test(dv = Alpha_centrality, wid = Subject, within = c(Concentration, Condition))
get_anova_table(res.aov_alpha_centrality)


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

seconds_graph_v3 <- data_seconds_2diag %>% 
  ggplot(aes(x=Time, y=Index, shape=Diagnosis, color=Condition,
             group =interaction(Condition, Subject)))+
  geom_point(size=7)+ geom_line() + scale_shape_manual(values=c(15,18)) +
  ylab('SECONDs score') +xlab('Time')+ theme_bw() +
  theme(text = element_text(size = 25))

seconds_graph_v3 + facet_grid(. ~ Subject) +  
  scale_y_continuous(breaks = seq(0,8,1), limits = c(-0.1, 8.1)) + geom_vline(xintercept = 4.5, linetype="dotted", 
                                                                            color = "black", size=2)

ggsave("KET_seconds_2diagn.tiff", device="tiff", width=30, height=20, units="cm", path=paper_path, dpi=300)

# 2. Figures on the MAS ####
setwd("C:/Users/Paolo/dox/PhD/01_Psychedelics/Ketamine/Pilot/results/Tables_dataset")

data_MAS <- read_excel("20230220_Database_KET_MAS.xlsx")
data_MAS[c('Condition', 'Time', 'Side', 'Limb')] <- 
  lapply(data_MAS[c('Condition', 'Time', 'Side', 'Limb')],as.factor)

# Statistics
res.aov_MAS <- data_MAS %>% 
  filter(Subject != 'KET01')%>% anova_test(dv = MAS, wid = Subject, within = c(Time,Condition, Side, Limb))
get_anova_table(res.aov_MAS)

#Change for visual representation
#temp1 <- data_MAS %>% filter(Condition == "Placebo") %>% mutate(MAS = MAS -0.05)
#temp2 <- data_MAS %>% filter(Condition == "Ketamine") %>% mutate(MAS = MAS +0.05)
#data_MAS <- rbind(temp1, temp2)
data_MAS$MAS <- data_MAS$MAS + seq(-0.1, 0.1, length=8)

data_MAS$Subject <- as.factor(data_MAS$Subject)
levels(data_MAS$Subject) <- c('MCS-','UWS','MCS+')

MAS_graph <- data_MAS %>% 
  filter(Subject != 'MCS-')%>% 
  ggplot(aes(x=Time, y=MAS, shape=Condition, color=Condition, label=Limb,
                                             group =interaction(Condition, Limb, Side)))+
  geom_point(size=5)+ geom_line() + 
  ylab('Modified Ashworth Scale') +xlab('Time')+ theme_bw() +
  theme(text = element_text(size = 25), 
        axis.text.x = element_text(size = 18)) +
  #geom_text(aes(label = paste(Limb, "")), size =3)+
  geom_text_repel(size=4)
  
#MAS_graph + facet_grid(Subject ~ Side)

#MAS_graph + facet_nested(Subject ~ Side + Condition)
MAS_graph + facet_nested(Condition ~ Subject + Side )
ggsave("KET_MASS_v2.tiff", device="tiff", width=30, height=20, units="cm", path=paper_path, dpi=300)
#MAS_graph + facet_nested(Subject ~ Limb + Side)

# Version 2 of figure
data_MAS_mod <- data_MAS %>% filter(Subject!='MCS-')
# data_MAS_mod %>% group_by(Subject, Condition, Side, Limb) %>% 
#  summarise(MAS_sum = sum(MAS)) %>% filter(MAS_sum==0)
data_MAS_mod <- as.data.frame(data_MAS_mod)
data_MAS_mod <- data_MAS_mod[-c(2, 10, 8, 16, 22, 30, 18, 26, 24, 32),]

# Per sub
subject <- 'KET03'
MAS_graph <- data_MAS_mod %>% filter(Subject==subject)%>%
  ggplot(aes(x=Time, y=MAS, shape=Condition, color=Condition, label=Limb,
             group =interaction(Condition, Limb, Side)))+
  geom_point(size=10)+ geom_line() + 
  ylab('Modified Ashworth Scale') +xlab('Time')+ theme_bw() +
  theme(text = element_text(size = 45)) +
  #geom_text(aes(label = paste(Limb, "")), size =3)+
  geom_text_repel(size=10)

#MAS_graph + facet_grid(Subject ~ Side)

MAS_graph + facet_nested(. ~ Side + Condition) + 
  guides(shape = FALSE, color = FALSE) +  scale_y_continuous(breaks = seq(0,5,1), limits = c(-0.5, 5.5))
ggsave(paste("KET_MASS_",".tiff",sep=subject), device="tiff", width=50, height=20, units="cm", path=paper_path, dpi=300)

# all together
MAS_graph <- data_MAS_mod %>% 
  ggplot(aes(x=Time, y=MAS, shape=Condition, color=Condition, label=Limb,
             group =interaction(Condition, Limb, Side)))+
  geom_point(size=5)+ geom_line() + 
  ylab('Modified Ashworth Scale') +xlab('Time')+ theme_bw() +
  theme(text = element_text(size = 25), 
        axis.text.x = element_text(size = 18)) +
  #geom_text(aes(label = paste(Limb, "")), size =3)+
  geom_text_repel(size=4)

MAS_graph + facet_nested(Subject ~ Side + Condition)
ggsave("KET_MASS_all_V2.tiff", device="tiff", width=30, height=20, units="cm", path=paper_path, dpi=300)
