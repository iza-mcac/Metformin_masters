library(biomaRt)
library(clusterProfiler)

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)

datasets_mart<-listDatasets(mart)

attributes_mart<-listAttributes(mart)

filters_mart<-listFilters(mart)

attributes <- c("ensembl_transcript_id_version","start_position",
                "end_position","strand","transcript_length",
                "chromosome_name","ucsc")
filters <- c("chromosome_name","start","end")

##########################################

chr_vector <-c("1","1","2","2","2","3","4","6","7","7","8","8","9","9","10","11",
               "11","11","11","12","12","13","14","16","16","17","17","17","17",
               "17","17","18","19","19","21","X")

selected_lnc_1mb_chr <-selected_lnc_1mb %>% mutate(chr_location =chr_vector)

list(selected_lnc_1mb_chr[1,])

all.genes <- list()

for (i in 1:dim(selected_lnc_1mb_chr)[1]) {
  #for(j in 1:dim(selected_lnc_1mb_chr)[2])
  values <- list(chromosome= selected_lnc_1mb_chr[i,6],
                 start= selected_lnc_1mb_chr[i,5],
                 end=selected_lnc_1mb_chr[i,4])
  all.genes[[i]] <- getBM(attributes=attributes,
                     filters=filters,
                     values=values, mart=mart)
}

aaaaaaaaaaa <-all.genes[[10]]


###################################################


#all_genes_all <-rbind(all.genes[[1]], all.genes[[2]], all.genes[[3]],
#                      all.genes[[4]], all.genes[[5]], all.genes[[6]],
#                      all.genes[[7]],all.genes[[8]], all.genes[[9]],
#                      all.genes[[10]], all.genes[[11]], all.genes[[12]],
#                      all.genes[[13]], all.genes[[14]], all.genes[[15]],
#                      all.genes[[16]], all.genes[[17]], all.genes[[18]],
#                      all.genes[[19]], all.genes[[20]], all.genes[[21]],
#                      all.genes[[22]], all.genes[[23]], all.genes[[24]],
#                      all.genes[[25]], all.genes[[26]], all.genes[[27]],
#                      all.genes[[28]], all.genes[[29]], all.genes[[30]],
#                      all.genes[[31]], all.genes[[32]], all.genes[[33]],
#                      all.genes[[34]], all.genes[[35]], all.genes[[36]])
#
#
aaaaaaaaaaaaaa <- as.data.frame(all.genes[[10]])

#############################

pairs_in_4series<-as.data.frame(str_split_fixed(cor_all_series_cor$pairs, " ", 2))

cor_all_series_cor_pairsep <- cor_all_series_cor %>% mutate(row = pairs_in_4series$V1, col = pairs_in_4series$V2)

everything_in_cor <-as.data.frame(c(pairs_in_4series$V1, pairs_in_4series$V2))

names(everything_in_cor) <- c("names")

#

all_in_cor_dic<-left_join(everything_in_cor, txtogenesv37,
                                     by = c("names"="transcript_name")) %>%
  dplyr::select(1,2,6, 8)

###################### é aqui o problema
####PROBLEMA
#filtrar um por um não tem o que fazer


cor_all_series_cor_pairsep_dic1 <- cor_all_series_cor_pairsep %>% 
  left_join(txtogenesv37, by = c("row" = "transcript_name")) %>%
  dplyr::select(1:10)

cor_all_series_cor_pairsep_dic2 <- cor_all_series_cor_pairsep_dic1 %>% 
  left_join(txtogenesv37, by = c("col" = "transcript_name")) %>%
  dplyr::select(1:11)

aaaaaaaaaaaaaaa <-all.genes[[10]]

testenoac <-aaaaaaaaaaaaaaa[-85,]


#caralhooo <-testenoac %>% left_join(txtogenesv37, by= c("ensembl_transcript_id_version"="transcript_id"))

aaaaaaaaaaaaaaaaaaaaaaaaa<-aaaaaaaaaaaaaaa %>% 
  filter(ensembl_transcript_id_version %in% cor_all_series_cor_pairsep_dic2$transcript_id.x |
           ensembl_transcript_id_version %in% cor_all_series_cor_pairsep_dic2$transcript_id.y)

aaaaaaaaaaa<-dplyr::filter(cor_all_series_cor_pairsep_dic2, 
       transcript_id.y %in% testenoac$ensembl_transcript_id_version)

######################################3

filtered_cors <- list()

for (i in 1:dim(selected_lnc_1mb_chr)[1]) {
  filtered_cors[[i]]<- filter(all.genes[[i]], 
                              ensembl_transcript_id_version %in% cor_all_series_cor_pairsep_dic2$transcript_id.x |
                                ensembl_transcript_id_version %in% cor_all_series_cor_pairsep_dic2$transcript_id.y)
}

#################################

filtered_cors[[5]]

all_cor_all <-rbind(filtered_cors[[1]], filtered_cors[[2]], filtered_cors[[3]],
                    filtered_cors[[4]], filtered_cors[[5]], filtered_cors[[6]],
                    filtered_cors[[7]],filtered_cors[[8]], filtered_cors[[9]],
                    filtered_cors[[10]], filtered_cors[[11]], filtered_cors[[12]],
                    filtered_cors[[13]], filtered_cors[[14]], filtered_cors[[15]],
                    filtered_cors[[16]], filtered_cors[[17]], filtered_cors[[18]],
                    filtered_cors[[19]], filtered_cors[[20]], filtered_cors[[21]],
                    filtered_cors[[22]], filtered_cors[[23]], filtered_cors[[24]],
                    filtered_cors[[25]], filtered_cors[[26]], filtered_cors[[27]],
                    filtered_cors[[28]], filtered_cors[[29]], filtered_cors[[30]],
                    filtered_cors[[31]], filtered_cors[[32]], filtered_cors[[33]],
                    filtered_cors[[34]], filtered_cors[[35]], filtered_cors[[36]])

filtered_cors[[36]]

finaltest <-all_cor_all %>% 
  left_join(txtogenesv37, by = c("ensembl_transcript_id_version"="transcript_id")) %>%
  dplyr::select(1:6,11,14)

finaltest_filt <- finaltest %>% filter(transcript_type != "lncRNA") 

#aaaaaaaaaaa<-dplyr::filter(cor_all_series_cor_pairsep_dic2, 
#                           !transcript_id.x %in% selected_lnc_1mb_chr$TXNAME &
#                             !transcript_id.x %in% selected_lnc_1mb_chr$TXNAME )
#

testestests <-aaaaaaaaaaa %>% filter(transcript_id.x %in% finaltest_filt$ensembl_transcript_id_version)

#cor_all_series_cor_pairsep_dic2 %>%
#  filter(transcript_id.y %in% all_cor_all$ensembl_transcript_id_version)
#

###################################################################

###############################

all_poscor_df_down <-as.data.frame(all_poscor_down_v)

all_poscor_df_down_joined<-left_join(all_poscor_df_down, txtogenesv37,
          by = c("all_poscor_down_v"="transcript_name")) %>%
  dplyr::select(1,3, 7)

teste <-all_poscor_df_down_joined %>%
  filter(transcript_id %in% all_genes_all$ensembl_transcript_id_version)


#########################################################

enrich_poscor_down <- enricher(gene = teste$gene_name, pvalueCutoff = 0.05,
                               pAdjustMethod = "BH", TERM2GENE = c2_all_terms_gmt)

enrich_poscor_down_df <-as.data.frame(enrich_poscor_down) %>% filter(pvalue < 0.05)

write.csv(enrich_poscor_down_df, "results_tables/enrich_poscor_down_gen.csv",
          row.names = FALSE)

#############################

all_poscor_df_up <-as.data.frame(all_poscor_up_v)

all_poscor_df_up_joined<-left_join(all_poscor_df_up, txtogenesv37,
                                     by = c("all_poscor_up_v"="transcript_name")) %>%
  dplyr::select(1,3,7)

teste <-all_poscor_df_up_joined %>%
  filter(transcript_id %in% all_genes_all$ensembl_transcript_id_version)

enrich_poscor_up <- enricher(gene = teste$gene_name, pvalueCutoff = 0.05,
                               pAdjustMethod = "BH", TERM2GENE = c2_all_terms_gmt)

enrich_poscor_up_df <-as.data.frame(enrich_poscor_up) %>% filter(pvalue < 0.05)

write.csv(enrich_poscor_up_df, "results_tables/enrich_poscor_up_gen.csv",
          row.names = FALSE)

#############################

all_negcor_df_up <-as.data.frame(all_negcor_up_v)

all_negcor_df_up_joined<-left_join(all_negcor_df_up, txtogenesv37,
                                   by = c("all_negcor_up_v"="transcript_name")) %>%
  dplyr::select(1,3,7)

teste <-all_negcor_df_up_joined %>%
  filter(transcript_id %in% all_genes_all$ensembl_transcript_id_version)

enrich_negcor_up <- enricher(gene = teste$gene_name, pvalueCutoff = 0.05,
                             pAdjustMethod = "BH", TERM2GENE = c2_all_terms_gmt)

enrich_negcor_up_df <-as.data.frame(enrich_negcor_up) %>% filter(pvalue < 0.05)

write.csv(enrich_negcor_up_df, "results_tables/enrich_negcor_up_gen.csv",
          row.names = FALSE)

#############################

all_negcor_df_down <-as.data.frame(all_negcor_down_v)

all_negcor_df_down_joined<-left_join(all_negcor_df_down, txtogenesv37,
                                   by = c("all_negcor_down_v"="transcript_name")) %>%
  dplyr::select(1,3,7)

teste <-all_negcor_df_down_joined %>%
  filter(transcript_id %in% all_genes_all$ensembl_transcript_id_version)

enrich_negcor_down <- enricher(gene = teste$gene_name, pvalueCutoff = 0.05,
                             pAdjustMethod = "BH", TERM2GENE = c2_all_terms_gmt)

enrich_negcor_down_df <-as.data.frame(enrich_negcor_down) %>% filter(pvalue < 0.05)

write.csv(enrich_negcor_down_df, "results_tables/enrich_negcor_down_gen.csv",
          row.names = FALSE)

