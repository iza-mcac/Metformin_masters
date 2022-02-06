library(Hmisc)

#################
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

#
df_yue_cor_nonames <- column_to_rownames(df_yue_cor)

yue_test_cor_res <-rcorr(as.matrix(df_yue_cor_nonames), type = "spearman")

yue_test_cor_res_P <-yue_test_cor_res$P

yue_test_cor_res_r<-yue_test_cor_res$r

Yue_flatten_test <- flattenCorrMatrix(yue_test_cor_res_r, yue_test_cor_res_P)

Yue_flatten_test_filtered <- Yue_flatten_test %>% filter(p<0.05 & abs(cor)> 0.6)

saveRDS(Yue_flatten_test_filtered, "results_tables/Yue_cor_res.rds")

rm(Yue_flatten_test)
rm(yue_test_cor_res_P)
rm(yue_test_cor_res_r)
rm(yue_test_cor_res)

##################################################################################

df_laustriat_10mm_cor_nonames <- column_to_rownames(df_laustriat_10mm_cor)

laustriat_10mm_test_cor_res <-rcorr(as.matrix(df_laustriat_10mm_cor_nonames),
                                    type = "spearman")

laustriat_10mm_cor_res_P <-laustriat_10mm_test_cor_res$P

laustriat_10mm_cor_res_r<-laustriat_10mm_test_cor_res$r


laustriat_10mm_test <- flattenCorrMatrix(laustriat_10mm_cor_res_r,
                                         laustriat_10mm_cor_res_P)

laustriat_10mm_filtered <- laustriat_10mm_test %>% filter(p<0.05 & abs(cor)> 0.6)

saveRDS(laustriat_10mm_filtered, "results_tables/laustriat_10mm_cor_res.rds")

rm(laustriat_10mm_test)
rm(laustriat_10mm_test_cor_res)
rm(laustriat_10mm_cor_res_P)
rm(laustriat_10mm_cor_res_r)


####################################################################3

df_laustriat_25mm_cor_nonames <- column_to_rownames(df_laustriat_25mm_cor)

laustriat_25mm_test_cor_res <-rcorr(as.matrix(df_laustriat_25mm_cor_nonames),
                                    type = "spearman")

laustriat_25mm_cor_res_P <-laustriat_25mm_test_cor_res$P

laustriat_25mm_cor_res_r<-laustriat_25mm_test_cor_res$r


laustriat_25mm_test <- flattenCorrMatrix(laustriat_25mm_cor_res_r,
                                         laustriat_25mm_cor_res_P)

laustriat_25mm_filtered <- laustriat_25mm_test %>% filter(p<0.05 & abs(cor)> 0.6)

rm(laustriat_25mm_test)
rm(laustriat_25mm_test_cor_res)
rm(laustriat_25mm_cor_res_P)
rm(laustriat_25mm_cor_res_r)

saveRDS(laustriat_25mm_filtered, "results_tables/laustriat_25mm_cor_res.rds")

#####################################################################3
#
df_Xie_kidney_cor_nonames <- column_to_rownames(df_Xie_kidney_cor)
#
Xie_kidney_cor_res <-rcorr(as.matrix(df_Xie_kidney_cor_nonames),
                           type = "spearman")
#
Xie_kidney_cor_res_P <-Xie_kidney_cor_res$P
#
Xie_kidney_cor_res_r<-Xie_kidney_cor_res$r
#
#
Xie_kidney_test <- flattenCorrMatrix(Xie_kidney_cor_res_r,
                                     Xie_kidney_cor_res_P)
#
Xie_kidney_filtered <- Xie_kidney_test %>% filter(p<0.05 & abs(cor)> 0.6)

rm(Xie_kidney_test)
rm(Xie_kidney_cor_res_r)
rm(Xie_kidney_cor_res_P)
rm(Xie_kidney_cor_res)
#
saveRDS(Xie_kidney_filtered, "results_tables/Xie_kidney_cor_res.rds")
#
#####################################################################3
#
df_Xie_lung_cor_nonames <- column_to_rownames(df_Xie_lung_cor)

Xie_lung_cor_res <-rcorr(as.matrix(df_Xie_lung_cor_nonames),
                         type = "spearman")

Xie_lung_cor_res_P <-Xie_lung_cor_res$P

Xie_lung_cor_res_r<-Xie_lung_cor_res$r


Xie_lung_test <- flattenCorrMatrix(Xie_lung_cor_res_r,
                                   Xie_lung_cor_res_P)

Xie_lung_filtered <- Xie_lung_test %>% filter(p<0.05 & abs(cor)> 0.6)

rm(Xie_lung_test)
rm(Xie_lung_cor_res_P)
rm(Xie_lung_cor_res_r)
rm(Xie_lung_cor_res)

saveRDS(Xie_lung_filtered, "results_tables/Xie_lung_cor_res.rds")
#
#####################################################################3

df_luizon_cor_nonames <- column_to_rownames(df_luizon_cor)

Luizon_cor_res <-rcorr(as.matrix(df_luizon_cor_nonames),
                       type = "spearman")

Luizon_cor_res_P <-Luizon_cor_res$P

Luizon_cor_res_r<-Luizon_cor_res$r


Luizon_test_flatten <- flattenCorrMatrix(Luizon_cor_res_r,
                                         Luizon_cor_res_P)

Luizon_test_flatten_filtered <- Luizon_test_flatten %>% filter(p<0.05 & abs(cor)> 0.6)

rm(Luizon_cor_res)
rm(Luizon_cor_res_P)
rm(Luizon_cor_res_r)
rm(Luizon_test_flatten)

saveRDS(Luizon_test_flatten_filtered, "results_tables/Luizon_cor_res.rds")

###########################################################################
#############################################################################
############################################################################

laustriat_10mm_test <- readRDS("results_tables/laustriat_10mm_cor_res.rds")

laustriat_25mm_test <-readRDS("results_tables/laustriat_25mm_cor_res.rds")

Luizon_test_flatten <-readRDS("results_tables/Luizon_cor_res.rds")

Xie_kidney_test <-readRDS("results_tables/Xie_kidney_cor_res.rds")

Xie_lung_test <-readRDS("results_tables/Xie_lung_cor_res.rds")

Yue_flatten_test <-readRDS("results_tables/Yue_cor_res.rds")

###########################################################################
#filter results of corr analysis

only_lncRNA <-txtogenesv37 %>% filter(transcript_type == "lncRNA")

yue_test_lncRNA <- Yue_flatten_test %>% filter(row %in% only_lncRNA$transcript_name |
                            column %in% only_lncRNA$transcript_name)

Xie_lung_lncRNA <- Xie_lung_test %>% filter(row %in% only_lncRNA$transcript_name |
                                                 column %in% only_lncRNA$transcript_name)

Xie_kidney_lncRNA <- Xie_kidney_test %>% filter(row %in% only_lncRNA$transcript_name |
                                              column %in% only_lncRNA$transcript_name)

luizon_lncRNA <- Luizon_test_flatten %>% filter(row %in% only_lncRNA$transcript_name |
                                                  column %in% only_lncRNA$transcript_name)

laustriat10mm_lncRNA <- laustriat_10mm_test %>% filter(row %in% only_lncRNA$transcript_name |
                                                  column %in% only_lncRNA$transcript_name)

laustriat25mm_lncRNA <- laustriat_25mm_test %>% filter(row %in% only_lncRNA$transcript_name |
                                                         column %in% only_lncRNA$transcript_name)

####################################################

saveRDS(yue_test_lncRNA, "yue_test_lncRNA_filt.rds")
saveRDS(Xie_lung_lncRNA, "Xie_lung_lncRNA_filt.rds")
saveRDS(Xie_kidney_lncRNA, "Xie_kidney_lncRNA_filt.rds")
saveRDS(luizon_lncRNA, "luizon_lncRNA_filt.rds")
saveRDS(laustriat10mm_lncRNA, "laustriat10mm_lncRNA_filt.rds")
saveRDS(laustriat25mm_lncRNA, "laustriat25mm_lncRNA_filt.rds")

yue_test_lncRNA<-readRDS("yue_test_lncRNA_filt.rds")
Xie_lung_lncRNA<-readRDS("Xie_lung_lncRNA_filt.rds")
Xie_kidney_lncRNA<-readRDS("Xie_kidney_lncRNA_filt.rds")
luizon_lncRNA<-readRDS("luizon_lncRNA_filt.rds")
laustriat10mm_lncRNA<-readRDS("laustriat10mm_lncRNA_filt.rds")
laustriat25mm_lncRNA<-readRDS("laustriat25mm_lncRNA_filt.rds")


####################
###########################################
###########################################
##########################################

laustriat25mm_lncRNA_pairs <-laustriat25mm_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, p) %>% filter(p < 0.01)

laustriat10mm_lncRNA_pairs <-laustriat10mm_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, p) %>% filter(p < 0.01)

luizon_lncRNA_pairs <-luizon_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, p)%>% filter(p < 0.01)

Xie_kidney_lncRNA_pairs <-Xie_kidney_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, p)%>% filter(p < 0.01)

Xie_lung_lncRNA_pairs <-Xie_lung_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, p)%>% filter(p < 0.01)

yue_test_lncRNA_pairs <-yue_test_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, p)%>% filter(p < 0.01)


left_join_1 <-Xie_lung_lncRNA_pairs %>% left_join(laustriat25mm_lncRNA_pairs, by = "pairs")

left_join_2 <- left_join_1 %>% left_join(Xie_kidney_lncRNA_pairs, by = "pairs")

left_join_3 <- left_join_2 %>% left_join(yue_test_lncRNA_pairs, by = "pairs")

left_join_4 <- left_join_3 %>% left_join(luizon_lncRNA_pairs, by = "pairs")

left_join_5 <- left_join_4 %>% left_join(laustriat10mm_lncRNA_pairs, by = "pairs")

###

countNA <- function(df) apply(df, MARGIN = 1, FUN = function(x) length(x[is.na(x)]))

na_count <- countNA(left_join_5)

cor_all_series <- left_join_5[na_count < 2,]

#############################################################
##########################

laustriat25mm_lncRNA_pairs_cor <-laustriat25mm_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, cor)

laustriat10mm_lncRNA_pairs_cor <-laustriat10mm_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, cor)

luizon_lncRNA_pairs_cor <-luizon_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, cor)

Xie_kidney_lncRNA_pairs_cor <-Xie_kidney_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, cor)

Xie_lung_lncRNA_pairs_cor <-Xie_lung_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, cor)

yue_test_lncRNA_pairs_cor <-yue_test_lncRNA %>%
  mutate(pairs = paste(row, column)) %>%
  select(pairs, cor)


left_join_1 <-Xie_lung_lncRNA_pairs_cor %>% left_join(laustriat25mm_lncRNA_pairs_cor, by = "pairs")

left_join_2 <- left_join_1 %>% left_join(Xie_kidney_lncRNA_pairs_cor, by = "pairs")

left_join_3 <- left_join_2 %>% left_join(yue_test_lncRNA_pairs_cor, by = "pairs")

left_join_4 <- left_join_3 %>% left_join(luizon_lncRNA_pairs_cor, by = "pairs")

left_join_5 <- left_join_4 %>% left_join(laustriat10mm_lncRNA_pairs_cor, by = "pairs")

###

countNA <- function(df) apply(df, MARGIN = 1, FUN = function(x) length(x[is.na(x)]))

na_count <- countNA(left_join_5)

cor_all_series_cor <- left_join_5[na_count <= 2,]

#####


cor_all_series_cor_inall <- left_join_5[na_count < 2,]

#Yue_flatten_test_joinedpairs <- Yue_flatten_test %>% mutate(pairs = paste(row, column))
#
#Yue_flatten_joinedpairs<- Yue_flatten_test_joinedpairs %>% dplyr::select(pairs, cor, p)
#
#
#laustriat_10mm_joinedpairs <- laustriat_10mm_test %>%
#  mutate(pairs = paste(row, column)) %>%
#  select(pairs, p)
#
#names(laustriat_10mm_joinedpairs) <- c("pairs", "pval_laustriat10mm")
#
#
#laustriat_25mm_joinedpairs <- laustriat_25mm_test %>%
#  mutate(pairs = paste(row, column)) %>%
#  select(pairs, p)
#
#names(laustriat_25mm_joinedpairs) <- c("pairs", "pval_laustriat25mm")
#
#
#xie_kidney_joinedpairs <- Xie_kidney_test %>%
#  mutate(pairs = paste(row, column)) %>%
#  select(pairs, p)
#
#names(xie_kidney_joinedpairs) <- c("pairs", "pval_xiekidney")
#
#
#xie_lung_joinedpairs <- Xie_lung_test %>%
#  mutate(pairs = paste(row, column)) %>%
#  select(pairs, p)
#
#names(xie_lung_joinedpairs) <- c("pairs", "pval_xielung")
#
#
#luizon_joinedpairs <- Luizon_test_flatten %>%
#  mutate(pairs = paste(row, column)) %>%
#  select(pairs, p)
#
#names(luizon_joinedpairs) <- c("pairs", "pval_luizon")
#