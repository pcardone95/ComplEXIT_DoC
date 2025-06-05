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
library(cowplot)
library(ggpubr)


# Statistics
library(tidyverse)
library(ggpubr)
library(rstatix)

# Set path and file_name
paper_path <- "..."
data_LZC_allcon_path <- 'pilotKET_merged_R_all.txt'
data_PerCon_path <- 'pilotKET_merged_R_PerConc.txt'
data_alpha_path <-  'pilotKET_merged_R_alpha.txt'
data_ECI_path <-  'ECI_table_R_bas_full.csv'

  
# 2. Lempel-Ziv Complexity ####
setwd("...")
data_LZC_allcon <- read.csv(data_LZC_allcon_path, header = TRUE, sep = ',')
data_LZC_allcon[c('Subject', 'Electrode', 'Concentration', 'Session', 'Condition')] <- 
  lapply(data_LZC_allcon[c('Subject', 'Electrode', 'Concentration', 'Session', 'Condition')],as.factor)
levels(data_LZC_allcon$Subject) <- c('Case 2: \nMCS-','Case 1: \nUWS', 'Case 3: \nMCS+')
data_LZC_allcon$Subject <- factor(data_LZC_allcon$Subject, levels =  c('Case 1: \nUWS','Case 2: \nMCS-','Case 3: \nMCS+'))

# 2.1 Whole-brain LZC - All trials s - Figure ####
data_LZC_allcon_summarize <-data_LZC_allcon %>%
  group_by(Subject, Condition) %>%
  summarise(mean_LZC = mean(LZC), sd_LZC = sd(LZC))

LZC_graph_grid <- data_LZC_allcon_summarize %>% ggplot(aes(x=factor(Condition, level=c('Placebo', 'Ketamine')), y=mean_LZC, color=Condition,
                                                           shape=Subject, group =interaction(Condition, Subject)))+
  geom_point(aes(fill=Condition), size=5)+ geom_line(aes(group =Subject), color = 'black') +
  geom_errorbar(aes(ymin=mean_LZC-sd_LZC/sqrt(127), ymax=mean_LZC+sd_LZC/sqrt(127)), width=.2) + 
  xlab('') + ylab('') + ggtitle('LZC')+ theme_bw() + theme(text = element_text(size = 20)) +
  scale_shape_manual(values=c(24,22,23),breaks=c('Case 1: \nUWS','Case 2: \nMCS-','Case 3: \nMCS+'),
                     labels = c('Case 1: UWS','Case 2: MCS-','Case 3: MCS+')) + guides(fill=guide_legend(title="Session"), colour=guide_legend(title="Session")) + theme(legend.position='none', plot.title = element_text(hjust = 0.5))

# 2.3  Whole-brain LZC - Per concentrations - Figure #### 
data_PerCon <- read.csv(data_PerCon_path, header = TRUE, sep = ',')
data_PerCon[c('Electrode', 'Concentration', 'Session', 'Condition')] <- 
  lapply(data_PerCon[c('Electrode', 'Concentration', 'Session', 'Condition')],as.factor)

data_PerCon_summarize <-data_PerCon %>%
  group_by(Subject, Condition, Concentration) %>%
  summarise(mean_LZC = mean(LZC), sd_LZC = sd(LZC))

data_PerCon_summarize$Subject <- as.factor(data_PerCon_summarize$Subject)
levels(data_PerCon_summarize$Subject) <- c('Case 2: \nMCS-','Case 1: \nUWS', 'Case 3: \nMCS+')
data_PerCon_summarize$Subject <- factor(data_PerCon_summarize$Subject, levels =  c('Case 1: \nUWS','Case 2: \nMCS-','Case 3: \nMCS+'))

LZC_graph_summary <- data_PerCon_summarize %>% 
  ggplot(aes(x=Concentration, y=mean_LZC,  color=Condition,
             group =interaction(Condition, Subject)))+
  geom_point(size=3)+ geom_line() +
  geom_errorbar(aes(ymin=mean_LZC-sd_LZC/sqrt(127), ymax=mean_LZC+sd_LZC/sqrt(127)), width=.2) + 
  ylab('') +xlab(bquote('Concentration ' (µg~mL^-1))) + theme_bw() + 
  theme(text = element_text(size =10))+ 
  scale_y_continuous(limits = c(0.18, 0.59)) + geom_vline(xintercept = 1.5, linetype="dashed", 
                                                        color = "black", size=1)
LZC_graph_summary_grid  <- LZC_graph_summary + facet_grid( .~ Subject ) + theme(legend.position='none', text = element_text(size = 15),
                                                                               axis.text.x = element_text(size=8, angle = 45))+ 
  guides(fill=guide_legend(title="Session"), colour=guide_legend(title="Session"))

# 3. ECI ####
data_ECI <- read.csv(data_ECI_path, header = TRUE, sep = ',')
data_ECI$Session <- as.factor(data_ECI$Session )
data_ECI$Subject <- as.factor(data_ECI$Subject )
levels(data_ECI$Subject) <- c('Case 2: \nMCS-','Case 1: \nUWS', 'Case 3: \nMCS+')
data_ECI$Subject <- factor(data_ECI$Subject, levels = c('Case 1: \nUWS', 'Case 2: \nMCS-','Case 3: \nMCS+'))


# 3.1 ECI - All trials - Figure ####
data_ECIAr_grid <- data_ECI %>% filter(Concentration == "All") %>% ggplot(aes(x=factor(Condition, level=c('Placebo', 'Ketamine')), y=ECI_Aro, 
                                                           color=Condition,
                                                           shape=Subject, group =interaction(Condition, Subject)))+
  geom_point(aes(fill=Condition),size=5)+ geom_line(aes(group =Subject), color = 'black') +
  xlab('') + ylab('') + ggtitle('ECI Arousal') +theme_bw() + theme(text = element_text(size = 20))+
  ylim(c(0,1)) +  geom_hline(yintercept = 0.5, linetype="dashed", color = "red",  size=1) +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values=c(24,22,23),breaks=c('Case 1: \nUWS','Case 2: \nMCS-','Case 3: \nMCS+'))+ guides(fill=guide_legend(title="Session"), colour=guide_legend(title="Session"))#+scale_colour_manual(values = c("blue3", "chartreuse4", "chocolate4"))

# Awa ##
data_ECIAw_grid <- data_ECI %>% filter(Concentration == "All") %>% ggplot(aes(x=factor(Condition, level=c('Placebo', 'Ketamine')), y=ECI_Awa,
                                                                              color=Condition,
                                                                              shape=Subject, group =interaction(Condition, Subject)))+
  geom_point(aes(fill=Condition),size=5)+ geom_line(aes(group =Subject), color = 'black') +
  xlab('') + ylab('') + ggtitle('ECI Awareness') +theme_bw() + theme(text = element_text(size = 20))+
  ylim(c(0,1)) +  geom_hline(yintercept = 0.5, linetype="dashed", color = "red",  size=1) +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values=c(24,22,23),breaks= c('Case 1: \nUWS','Case 2: \nMCS-','Case 3: \nMCS+'))+ guides(fill=guide_legend(title="Session"), colour=guide_legend(title="Session"))#+scale_colour_manual(values = c("blue3", "chartreuse4", "chocolate4"))


# 3.4 ECI - Per concentration- Figure ####
ECI_aro_graph <- data_ECI %>% filter(Concentration != 'All') %>% #filter(Concentration != '0')%>% 
  ggplot(aes(x=Concentration, y=ECI_Aro, color=Condition,
             group =interaction(Condition, Subject)))+
  geom_point(size=3)+ geom_line() + 
  ylab('')+ xlab(bquote('Concentration ' (µg~mL^-1))) + theme_bw() +
  theme(text = element_text(size = 10))+#, axis.text.x = element_text(size = 18)) +
  ylim(c(0,1)) +geom_hline(yintercept = 0.5, linetype="dashed",color = "red", size=1)+
  geom_vline(xintercept = 1.5, linetype="dashed", color = "black", size=1) 
ECI_aro_graph_grid <- ECI_aro_graph + facet_grid(. ~ Subject) +theme(legend.position = 'none', text = element_text(size = 15),
                                                                     axis.text.x = element_text(size=8, angle = 45)) + 
  guides(fill=guide_legend(title="Session"), colour=guide_legend(title="Session"))

ECI_awa_graph <- data_ECI %>% filter(Concentration != 'All') %>% #filter(Concentration != '0')%>% 
  ggplot(aes(x=Concentration, y=ECI_Awa, color=Condition,
             group =interaction(Condition, Subject)))+
  geom_point(size=3)+ geom_line() + 
  ylab('') +xlab(bquote('Concentration ' (µg~mL^-1)))+ theme_bw() +
  theme(text = element_text(size = 10))+#, axis.text.x = element_text(size = 18)) +
  ylim(c(0,1)) +  geom_hline(yintercept = 0.5, linetype="dashed", color = "red",  size=1) + 
  geom_vline(xintercept = 1.5, linetype="dashed", color = "black", size=1)
ECI_awa_graph_grid <- ECI_awa_graph + facet_grid( .~ Subject ) + theme(legend.position = 'none', text = element_text(size = 15),
                                                                       axis.text.x = element_text(size=8, angle = 45)) + 
  guides(fill=guide_legend(title="Session"), colour=guide_legend(title="Session"))


# 4. Figure 3 - EEG results####
# In common.legend. It is actually the first that is taken
ggarrange(LZC_graph_grid, data_ECIAr_grid, data_ECIAw_grid, 
          LZC_graph_summary_grid, ECI_aro_graph_grid, ECI_awa_graph_grid, 
          nrow = 2, ncol = 3,  labels = c('A', 'B', 'C'), font.label = list(size = 28),  
          common.legend = TRUE, legend = "right")
ggsave("PaperKET_EEGresults_iScience_r1_sub.tiff", device="tiff", width=15, height=8.25, units="in", path=paper_path, dpi=300)
ggsave("PaperKET_EEGresults_iScience_r1_sub.jpg", device="jpg", width=15, height=8.25, units="in", path=paper_path, dpi=300)



# 5 SECONDs ####
setwd("...")

data_seconds <- read_excel("20230220_SECONDs_R.xlsx", 
                           col_types = c("text", "text", "text", "numeric", "text"))
data_seconds[c('Condition', 'Time', 'Diagnosis')] <- 
  lapply(data_seconds[c('Condition', 'Time', 'Diagnosis')],as.factor)
data_seconds$Time <- ordered(data_seconds$Time, levels = c("0", "30", "60", "90", "120")) 
levels(data_seconds$Time) <- c("0", "30", "60", "90", "210")

# 5.1 SECONDs - FIGURE ####
#Change for visual representation
temp1 <- data_seconds %>% filter(Condition == "Placebo") %>% mutate(Index = Index -0.1)
temp2 <- data_seconds %>% filter(Condition == "Ketamine") %>% mutate(Index = Index +0.1)
data_seconds <- rbind(temp1, temp2)

data_seconds$Subject <- as.factor(data_seconds$Subject)
levels(data_seconds$Subject) <-  c('Case 2: \nMCS-','Case 1: \nUWS', 'Case 3: \nMCS+')
data_seconds$Subject <- factor(data_seconds$Subject, levels =  c('Case 1: \nUWS','Case 2: \nMCS-','Case 3: \nMCS+'))

seconds_graph_v2 <- data_seconds %>% ggplot(aes(x=Time, y=Index, shape=Diagnosis, color=Condition,
                                             group =interaction(Condition, Subject)))+
  geom_point(aes(fill=Condition),size=7) + geom_line() + 
  ggtitle('Awareness') +ylab('SECONDs score') +xlab('Time (min)')+ theme_bw() +
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 18))   +  
  scale_y_continuous(breaks = seq(0,8,1), limits = c(-0.1, 8.1), minor_breaks = seq(0,8,1)) +
  geom_vline(xintercept = 1.5, linetype="dashed", color = "black", size=2)+
  scale_shape_manual(values=c(24,22,23),breaks=c('UWS','MCS-','MCS+')) + 
  geom_vline(xintercept = 4.5, linetype="dashed", color = "gray", size=2)

seconds_graph_grid <- seconds_graph_v2 + facet_grid(. ~ Subject) + 
  guides(fill=guide_legend(title="Session"), colour=guide_legend(title="Session"),
         shape = guide_legend(title="Diagnosis at \nassessment"))

seconds_graph_grid <- seconds_graph_grid +  
  geom_text(aes(label=c(rep(c(" "), times=3), c(".75", "0"), rep(c(" "), times=7), 
                        c(".75", ".75", " "), rep(c("0", ".45", ".75", ".75", "0"), times=3))), 
            color='black', fontface = "bold", size = 3)


# 7. Eyes Opened - Figure ####
data_aro <- read_excel("20231117_Arousal_R.xlsx", 
                           col_types = c("text", "text", "text", "numeric", "text"))
data_aro[c('Session', 'Time', 'Diagnosis')] <- 
  lapply(data_aro[c('Session', 'Time', 'Diagnosis')],as.factor)
data_aro$Time <- ordered(data_aro$Time, levels = c("0", "30", "60", "90", "210")) 


# Value=1: 0-25; 2: 25-50; 3: 50-75; 4: 75-100
temp1 <- data_aro %>% filter(Session == "Placebo") %>% mutate(Index = Index -0.1)
temp2 <- data_aro %>% filter(Session == "Ketamine") %>% mutate(Index = Index +0.1)
data_aro <- rbind(temp1, temp2)

data_aro$Subject <- as.factor(data_aro$Subject)
levels(data_aro$Subject) <- c('Case 2: \nMCS-','Case 3: \nMCS+', 'Case 1: \nUWS')
data_aro$Subject <- factor(data_aro$Subject, levels =  c('Case 1: \nUWS','Case 2: \nMCS-','Case 3: \nMCS+'))

aro_grid <-data_aro %>% ggplot(aes(x=Time, y=Index, shape=Diagnosis, color=Session,
                                                group =interaction(Session, Subject)))+
  geom_point(aes(fill=Session),size=7) + geom_line() + 
  ggtitle('Arousal') +ylab('Time w/ EO') +xlab('Time (min)')+ theme_bw() +
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 18))   +  
  scale_y_continuous(breaks = seq(0,8,1), limits = c(0.9, 4.1)) + geom_vline(xintercept = 1.5, linetype="dashed", 
                                                                              color = "black", size=2)+
  scale_shape_manual(values=c(24,22,23),breaks=c('UWS','MCS-','MCS+')) + 
  scale_y_continuous(
    breaks = c(1,2,3,4), minor_breaks = seq(1,4,1),
    labels = c("0-25%", "25-50%", "50-75%", "75-100%")
  )+
  geom_vline(xintercept = 4.5, linetype="dashed", color = "gray", size=2) +facet_grid(. ~ Subject) + 
  geom_text(aes(label=rep(c("0", ".45", ".75", ".75", "0"), times=6)), 
             color='black', 
            fontface = "bold", size = 3) + guides(shape = guide_legend(title="Diagnosis at \nassessment"))

# 7. MAS - Figure ####
setwd("...")

data_MAS <- read_excel("20230220_Database_KET_MAS.xlsx")
data_MAS[c('Condition', 'Time', 'Side', 'Limb')] <- 
  lapply(data_MAS[c('Condition', 'Time', 'Side', 'Limb')],as.factor)

#Change for visual representation
data_MAS$MAS <- data_MAS$MAS + seq(-0.1, 0.1, length=8)

data_MAS$Subject <- as.factor(data_MAS$Subject)
levels(data_MAS$Subject) <-  c('Case 2: \nMCS-','Case 1: \nUWS', 'Case 3: \nMCS+')

### Figure with no "0" MAS
data_MAS_mod <- data_MAS %>% filter(Subject!='Case 2: \nMCS-')
data_MAS_mod <- as.data.frame(data_MAS_mod)
data_MAS_mod <- data_MAS_mod[-c(2, 10, 8, 16, 22, 30, 18, 26, 24, 32),]

MAS_graph_mod <- data_MAS_mod %>% 
  ggplot(aes(x=Time, y=MAS, color=Condition, label=Limb,
             group =interaction(Condition, Limb, Side)))+
  geom_point(size=5)+ geom_line() + #'Modified Ashworth Scale score'
  ylab('MAS score') +xlab('Time (min)')+ ggtitle('Spasticity') + theme_bw() +
  theme(text = element_text(size = 25), 
        axis.text.x = element_text(size = 18),
        plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(size=4)

MAS_graph_mod_grid <- MAS_graph_mod + theme(legend.position='none', axis.text.y = element_text(hjust = 0)) + 
  facet_nested(Condition ~ Subject + Side) + scale_y_continuous(minor_breaks = seq(0,5,1),
                                                                labels = c("0", "1", "1+", "2", "3", "4"))


# 8. - Figure 2 - Behav results ####
ggarrange(
  ggarrange(seconds_graph_grid, 
          aro_grid,
          labels = c('A', 'B'), font.label = list(size = 30), nrow = 2, 
          common.legend = TRUE, legend = "right", align = 'v'),
  MAS_graph_mod_grid, labels=c('', 'C'), font.label = list(size = 30), nrow = 2,
  heights = c(2,1))

ggsave("PaperKET_Behavresults_iScience_r1_sub.tiff", device="tiff", width=14, height=17, units="in", path=paper_path, dpi=300)
ggsave("PaperKET_Behavresults_iScience_r1_sub.jpg", device="jpg", width=14, height=17, units="in", path=paper_path, dpi=300)


# Supplementary material ########
# Alpha centrality ########
data_alpha <- read.csv(data_alpha_path, header = TRUE, sep = ',')

colnames(data_alpha)[3] <- "Concentration"
data_alpha[c('Concentration')] <- lapply(data_alpha[c('Concentration')],as.factor)
levels(data_alpha$Concentration) <- c("0", "0.15-0.30", "0.45-0.60", "0.75")

data_alpha$Subject <- as.factor(data_alpha$Subject)
levels(data_alpha$Subject) <- c('MCS-','UWS', 'MCS+')
data_alpha$Subject <- factor(data_alpha$Subject, levels =  c('UWS','MCS-','MCS+'))
levels(data_alpha$Subject) <- c('Case 1: \nUWS','Case 2: \nMCS-','Case 3: \nMCS+')

alpha_graph <- data_alpha %>% ggplot(aes(x=Concentration, y=Alpha,  color=Condition,
                                         group =interaction(Condition, Subject)))+
  geom_point(size=5)+ geom_line() + 
  ylab('Alpha Centrality') +xlab(bquote('Concentration ' (µg~mL^-1)))+ theme_bw() +
  theme(text = element_text(size = 25), axis.text.x = element_text(size = 11, angle =20))+ 
  geom_vline(xintercept = 1.5, linetype="dashed", color = "black", size=2)
alpha_graph + facet_grid(. ~ Subject) + guides(fill=guide_legend(title="Session"), colour=guide_legend(title="Session"))

ggsave("PaperKET_Alpha_iScience_r1.tiff", device="tiff", width=37, height=20, units="cm", path=paper_path, dpi=300)

# Correlation LZC and ECI ########
data_ECI_temp <- data_ECI %>% filter(Concentration != 'All')
total_EEG <- merge(data_PerCon_summarize, data_ECI_temp, by=c("Subject","Condition","Concentration"))

# Version 5
label.y.val <-  1
LZC_Awa <- total_EEG %>% ggscatter(x = "mean_LZC", y = "ECI_Awa", color="Condition", size=5, add = 'reg.line',conf.int = FALSE) +
  stat_cor(method = "spearman", cor.coef.name='rho', p.accuracy = 0.001, size =10, label.y = label.y.val) +theme_bw() +
  xlab('Whole-Brain LZC') + ylab('ECI Awareness') +theme(text = element_text(size = 20)) +
  guides(colour=guide_legend(title="Session"), fill=guide_legend(title="Session")) +
  geom_smooth(method = "lm", se = FALSE, colour = 'Black')

LZC_Aro <- total_EEG %>% ggscatter(x = "mean_LZC", y = "ECI_Aro", color="Condition", size=5,  add = 'reg.line', conf.int = FALSE) + 
  stat_cor(method = "spearman", cor.coef.name='rho', p.accuracy = 0.001, size =10, label.y = label.y.val) +theme_bw() +
  xlab('Whole-Brain LZC') + ylab('ECI Arousal') +theme(text = element_text(size = 20)) +
  guides(colour=guide_legend(title="Session")) +
  geom_smooth(method = "lm", se = FALSE, colour = 'Black')

ggarrange(LZC_Awa, LZC_Aro, nrow=1, common.legend = TRUE,legend = "right")
ggsave("LZC_ECI_corr.jpg", device="jpg", units="in", path=paper_path, dpi=300)

