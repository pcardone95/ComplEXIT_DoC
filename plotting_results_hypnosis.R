library('dplyr')
library('ggplot2')
library('data.table')
setwd("~/GitHub/ComplEXIT_DoC")

####### Whole brain ########

data <- read.csv('results_LZC_hypnosis_complete_v2.txt', header = TRUE, sep = ',')
data[c('Electrode')] <- lapply(data[c('Electrode')],as.factor)



data_summarize <-data %>% 
  group_by(Subject, Condition) %>%
  summarise(mean_LZC = mean(LZC), sd_LZC = sd(LZC))
data_summarize['Condition']<- lapply(data_summarize['Condition'],as.factor)
data_summarize['Subject']<- lapply(data_summarize['Subject'],as.factor)


####### Amy ########
#png("LZC_Hyp_Aminata_onlycond.jpg", height = 600, width = 600)

data_summarize %>%  ggplot(aes(x=Condition, y=mean_LZC, group = Subject, shape = Condition, color=Subject))+
  geom_point(size=5)+ geom_line() +
  #geom_errorbar(aes(ymin=mean_LZC-sd_LZC/sqrt(127), ymax=mean_LZC+sd_LZC/sqrt(127)), width=.2) + 
  ylab('LZC') +xlab('Condition')+ theme_bw() + theme(text = element_text(size = 30
  ), 
  plot.title = element_text(hjust = 0.5)) + guides(color=FALSE)

#dev.off()
