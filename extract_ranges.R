

BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
library(tidyverse)

download.file("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz",
              "gencode.v37.comp_assembly.annotation.gtf.gz")

txdbv37 <- makeTxDbFromGFF('gencode.v37.comp_assembly.annotation.gtf.gz')

transcripts_v37 <- transcripts(txdbv37)

heatmap_isoforms <- read.csv("references/heatmap_isoforms.csv")

txtogenesv37_heatmap_only <-txtogenesv37 %>% 
  filter(transcript_name %in% 
           heatmap_isoforms$rownames.swish_results_only_6s_heatmap_reorderd.)


vector_of_transcripts <-txtogenesv37_heatmap_only$transcript_id

#########################

#filter a txdb by a list of transcripts I want on final table

columns(txdbv37)

keys <- vector_of_transcripts

cols_i_want <- c("TXSTART", "TXEND", "TXSTRAND", "TXCHROM")

lncRNAisos_start_and_end <-biomaRt::select(txdbv37, keys=keys, columns=cols_i_want, keytype="TXNAME")

#############################
#################
#teste com 1 megabase

lncRNAisos_start_and_end_1mb <-lncRNAisos_start_and_end %>%
  mutate(added_value =1000000)

lncRNAisos_start_and_end_1mb_s <- lncRNAisos_start_and_end_1mb %>% 
  dplyr::select(-5) %>%
  mutate(new_start = c(TXSTART-added_value)) %>% 
  mutate(TXEND = lncRNAisos_start_and_end_1mb$TXEND)

lncRNAisos_start_and_end_1mb_e <- lncRNAisos_start_and_end_1mb_s %>% 
  dplyr::select(-6) %>%
  mutate(new_end = c(TXEND+added_value)) %>%
  mutate(new_start = lncRNAisos_start_and_end_1mb_s$new_start)

selected_lnc_1mb <- lncRNAisos_start_and_end_1mb_e %>% dplyr::select(1:3, 7,8)


selected_granges_1mb <-makeGRangesFromDataFrame(selected_lnc_1mb,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("TXCHROM"),
                         start.field="new_start",
                         end.field="new_end",
                         strand.field="TXSTRAND",
                         starts.in.df.are.0based=FALSE)


###############################################
library(GenomicRanges)
library(IRanges)
library(data.table)

#transcripts_v37_df <-as.data.frame(transcripts_v37)
#selected_granges_1mb_df <- as.data.frame(selected_granges_1mb)

#transcripts_v37_df %>% select()

test_2 <-findOverlaps(selected_granges_1mb, transcripts_v37,
                      type=c("within"),
                      select=c("all"))
test_2_ <-countOverlaps(selected_granges_1mb, transcripts_v37)


vector_lnc_names <-selected_lnc_1mb %>%
  left_join(txtogenesv37, by = c("TXNAME"="transcript_id")) %>%
  select(9)



barplot(test_2_, names.arg = vector_lnc_names$transcript_name,
        main="Possible targets in 1mb",ylab="Count",las=2)


#test_4 <- as.data.frame(transcripts_v37[subjectHits(test_2)])


#aaaaaaaaaaaa <-as.data.frame(subsetByOverlaps(transcripts(txdbv37),
                                             # selected_granges_1mb))


#

#overlaps <- pintersect(transcripts_v37[subjectHits(test_2)],
                       #selected_granges_1mb[queryHits(test_2)])


#overlaps_test<-as.data.frame(overlaps)


test_res <-left_join(test_4, txtogenesv37, by = c("tx_name"="transcript_id"))


#test <-intersect(selected_granges_1mb,transcripts_v37)

##################################

xie_lung_lncRNA_heatmap_poscor_up %>% 
  filter(row %in% test_res$transcript_name | column %in% test_res$transcript_name)





