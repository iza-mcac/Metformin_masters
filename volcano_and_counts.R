library(tidyverse)
library(readr)

####################################################################################
#Kidney

Xie_difexpr_kidney_results_dic <- read_csv("Xie_kidney_swish.csv")

#############################################
library(ggplot2)
library(ggthemes)

Xie_kidney_results_dic_mod<- Xie_difexpr_kidney_results_dic

# add a column of NAs
Xie_kidney_results_dic_mod$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
Xie_kidney_results_dic_mod$diffexpressed[Xie_kidney_results_dic_mod$log2FC > 0.5 &
                                           Xie_kidney_results_dic_mod$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
Xie_kidney_results_dic_mod$diffexpressed[Xie_kidney_results_dic_mod$log2FC < -0.5 &
                                           Xie_kidney_results_dic_mod$pvalue < 0.05] <- "DOWN"


Xie_kidney_results_dic_mod %>% 
  filter(abs(log2FC)> 0.01) %>%
  ggplot(aes(x= log2FC, y= -log10(pvalue), col = diffexpressed))+
  geom_point()+
  scale_color_manual(values= c("blue", "grey", "red"))+
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


#####################################
#barplot

xie_kidney_DTEs_dic <- Xie_kidney_results_dic_mod %>% filter(pvalue <0.05 & abs(log2FoldChange) > 0.5)

#DTEs_lncRNA <- xie_kidney_DTEs_dic %>% filter(transcript_type == "lncRNA")

xie_kidney_DTEs_dic %>%
  filter(transcript_type %in% c("protein_coding", "retained_intron", "nonsense_mediated_decay",
                                "processed_transcript", "lncRNA")) %>%
  ggplot(aes(x=transcript_type))+
  geom_bar(aes(fill = diffexpressed))+
  scale_fill_manual(values= c("blue", "red3"))
#theme(axis.text.x = element_text(angle = 45))


########################################################
#####################

library(tidyverse)

####################################################################################
#lung

Xie_difexpr_lung_results_dic <- read_csv("Xie_lung_swish.csv")

#############################################
library(ggplot2)
library(tidyverse)
library(ggthemes)

Xie_lung_results_dic_mod<- Xie_difexpr_lung_results_dic

# add a column of NAs
Xie_lung_results_dic_mod$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Xie_lung_results_dic_mod$diffexpressed[Xie_lung_results_dic_mod$log2FC > 0.5 &
                                         Xie_lung_results_dic_mod$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Xie_lung_results_dic_mod$diffexpressed[Xie_lung_results_dic_mod$log2FC < -0.5 &
                                           Xie_lung_results_dic_mod$pvalue < 0.05] <- "DOWN"


Xie_lung_results_dic_mod %>% 
  #filter(-log10(pvalue))%>%
  filter(abs(log2FC)> 0.01) %>%
  ggplot(aes(x= log2FC, y= -log10(pvalue), col = diffexpressed))+
  geom_point()+
  scale_color_manual(values= c("blue", "grey", "red"))+
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


#####################################
#barplot

xie_lung_DTEs_dic <- Xie_lung_results_dic_mod %>% filter(pvalue <0.05 & abs(log2FoldChange) > 0.5)

#DTEs_lncRNA <- xie_kidney_DTEs_dic %>% filter(transcript_type == "lncRNA")

xie_lung_DTEs_dic %>%
  filter(transcript_type %in% c("protein_coding", "retained_intron", "nonsense_mediated_decay",
                                "processed_transcript", "lncRNA")) %>%
  ggplot(aes(x=transcript_type))+
  geom_bar(aes(fill = diffexpressed))+
  scale_fill_manual(values= c("blue", "red3"))
#theme(axis.text.x = element_text(angle = 45))

