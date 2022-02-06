library(readr)
library(tidyverse)
library(tximport)
library(stringr)

txtogenesv37 <- read_csv("C:/Users/izamamede/Desktop/Masters/R_analysis/gencode37/txtogenesv37.csv",
                         col_types = cols(X1 = col_skip()))


xie_kidney_index <- read_delim("xie_kidney_index.txt", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)


files <- file.path("C:/Users/izamamede/Desktop/Masters/Quant_files/5.Xie",
                   paste0(xie_kidney_index$X1,"_quant"), "quant.sf")


names(files) <- xie_kidney_index$X1

all(file.exists(files))

##############################################################################

transcript_id_togenename <- txtogenesv37 %>%
  dplyr::select(2,7)

tximporteddata <- tximport(files, type="salmon", txIn = TRUE, txOut = TRUE,
                           countsFromAbundance = "no")

tximporteddata$abundance

x <-as.data.frame(txi_counts)

as.data.frame(txi_counts) <- tximporteddata$counts

xie_kidney_counts <- txi_counts

write.csv(xie_kidney_counts, "xie_kidney_counts.csv")


#####################################

##########################################

xie_kidney_index_mod <- xie_kidney_index %>% dplyr::select(1,3)

library("DESeq2")

ddsTxi <- DESeqDataSetFromTximport(tximporteddata,
                                   colData = xie_kidney_index_mod,
                                   design = ~ X3)

res()

dds$Metformin

dds <- DESeq(ddsTxi)

resultsNames(dds)

test <- results(dds)

test <- as.data.frame(test)

# or to shrink log fold changes association with condition:
#10Mm Met
res_shr <- lfcShrink(dds, coef = "X3_metformin_vs_control", type="apeglm")

#####################################################

Xie_difexpr_kidney <- as.data.frame(res_shr)

Xie_difexpr_kidney_results <- rownames_to_column(Xie_difexpr_kidney)

Xie_difexpr_kidney_results <- Xie_difexpr_kidney_results %>% filter(baseMean > 0) 


Xie_difexpr_kidney_results_dic <- left_join(Xie_difexpr_kidney_results, txtogenesv37, by = c("rowname"= "transcript_id"))


Xie_difexpr_kidney_results_dic <- Xie_difexpr_kidney_results_dic %>% select(1:7, 10:13)

write.csv(Xie_difexpr_kidney_results_dic, "Xie_difexpr_kidney_results_dic.csv")



#####################
######################
#####################
#lung
##################################################################################
library(readr)
library(tidyverse)
library(tximport)
library(stringr)

txtogenesv37 <- read_csv("C:/Users/izamamede/Desktop/Masters/R_analysis/gencode37/txtogenesv37.csv",
                         col_types = cols(X1 = col_skip()))


xie_lung_index <- read_delim("xie_lung_index.txt", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)


files <- file.path("C:/Users/izamamede/Desktop/Masters/Quant_files/5.Xie",
                   paste0(xie_lung_index$X1,"_quant"), "quant.sf")


names(files) <- xie_lung_index$X1

all(file.exists(files))

##############################################################################


tximporteddata <- tximport(files, type="salmon", txIn = TRUE, txOut = TRUE,
                           countsFromAbundance = "no")

tximporteddata$abundance

txi_counts <- tximporteddata$counts

xie_lung_counts <- txi_counts

write.csv(xie_lung_counts, "xie_lung_counts.csv")


#####################################

##########################################

xie_lung_index_mod <- xie_lung_index %>% select(1,3)

library("DESeq2")

ddsTxi <- DESeqDataSetFromTximport(tximporteddata,
                                   colData = xie_lung_index_mod,
                                   design = ~ X3)

dds <- DESeq(ddsTxi)

resultsNames(dds)


# or to shrink log fold changes association with condition:
#Lung
res_shr <- lfcShrink(dds, coef = "X3_metformin_vs_control", type="apeglm")

#####################################################

Xie_difexpr_lung <- as.data.frame(res_shr)

Xie_difexpr_lung_results <- rownames_to_column(Xie_difexpr_lung)

Xie_difexpr_lung_results <- Xie_difexpr_lung_results %>% filter(baseMean > 0) 


Xie_difexpr_lung_results_dic <- left_join(Xie_difexpr_lung_results, txtogenesv37, by = c("rowname"= "transcript_id"))


Xie_difexpr_lung_results_dic <- Xie_difexpr_lung_results_dic %>% select(1:7, 10:13)

write.csv(Xie_difexpr_lung_results_dic, "Xie_difexpr_lung_results_dic.csv")



