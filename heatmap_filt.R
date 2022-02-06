library(pheatmap)
library(viridis)
library(RColorBrewer)
library(readr)


#
colors_list_filt <- list(
  Cell_type = c('PANC-1'="#1b9e77", HEsC ="#d95f02", kidney_tissue ="#7570b3",
                Lung_tissue ="#e7298a", Primary_hep ="#e6ab02"),
  Met_concentration = c('2.5 mM'="#edf8fb", '5 mM'="#ccece6",
                        '10 mM'="#99d8c9",'25 mM'="#66c2a4", '32.27 mM'="#2ca25f",
                        '69.97 mM'="#006d2c"),
  Read_length = c('100' = "#41b6c4", '150' = "#253494"),
  Hours_post_treatment = c('8' = "#efedf5", '48' = "#bcbddc",
                           '72' = "#756bb1"))





#

all_4plot_5_series<- all_joined_present_in_5_series_at_least %>% 
  select(-2)

all_4plot_5_series <- as.data.frame(all_4plot_5_series)

rownames(all_4plot_5_series) <- all_joined_present_in_5_series_at_least$transcript_name

all_4plot_5_series[is.na(all_4plot_5_series)] <- 0


pheatmap_output <- pheatmap(mat = all_4plot_5_series, color = viridis(3),
                            #annotation_col = lib_index_filt,
                            show_rownames = TRUE, breaks = c(5, 0.5, 0, -0.5, -5),
                            annotation_colors = colors_list)


#####################################################

lib_index_filt_2 <- read_delim("lib_index_filt.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

lib_index_filt_2 <- column_to_rownames(lib_index_filt_2, var = "X1")


lib_index_filt_3<- lib_index_filt_2 %>% mutate(order = c(2,6,5,1,3,4))

lib_index_filt_2_reord <- lib_index_filt_3 %>% arrange(order)

lib_index_no_series_type <- lib_index_filt_2_reord %>% select(-6, -1)


test <- brewer.pal(n = 11, name = 'RdBu')



test <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7",
         "#F4A582", "#D6604D", "#B2182B", "#67001F")



all_4plot_5_series_reorderd<- data.frame(all_4plot_5_series$log2FC_Luizon,
                                            all_4plot_5_series$log2FC_Yue,
                                            all_4plot_5_series$log2FC_10Mm_Laustriat,
                                            all_4plot_5_series$log2FC_25Mm_Laustriat,
                                            all_4plot_5_series$log2FC_lung_Xie,
                                            all_4plot_5_series$log2FC_kidney_Xie)

rownames(all_4plot_5_series_reorderd) <- rownames(all_4plot_5_series)

names(all_4plot_5_series_reorderd)<- c("log2FC_Luizon", "log2FC_Yue",
                                       "log2FC_10Mm_Laustriat", "log2FC_25Mm_Laustriat",
                                       "log2FC_lung_Xie", "log2FC_kidney_Xie")

pheatmap_output <- pheatmap(mat = all_4plot_5_series_reorderd, color = test,
                            annotation_col = lib_index_no_series_type,
                            annotation_colors = colors_list_filt,
                            show_rownames = TRUE,
                            breaks = c(-5,-4,-3,-2,-1, -0.5, 0, 0.5,1,2,3,4, 5),
                            border_color = "grey", cluster_cols = FALSE,
                            clustering_distance_rows = "canberra")


####################################################################################

all_4plot_4_series_filt <- all_joined_present_in_4_series_at_least_filt %>% 
  select(-2)

all_4plot_4_series_filt <- as.data.frame(all_4plot_4_series_filt)

rownames(all_4plot_4_series_filt) <- all_joined_present_in_4_series_at_least_filt$transcript_name

all_4plot_4_series_filt[is.na(all_4plot_4_series_filt)] <- 0


all_4plot_4_series_filt_reorderd<- data.frame(all_4plot_4_series_filt$log2FC_Luizon,
                                         all_4plot_4_series_filt$log2FC_Yue,
                                         all_4plot_4_series_filt$log2FC_10Mm_Laustriat,
                                         all_4plot_4_series_filt$log2FC_25Mm_Laustriat,
                                         all_4plot_4_series_filt$log2FC_lung_Xie,
                                         all_4plot_4_series_filt$log2FC_kidney_Xie)

rownames(all_4plot_4_series_filt_reorderd) <- rownames(all_4plot_4_series_filt)

names(all_4plot_4_series_filt_reorderd)<- c("log2FC_Luizon", "log2FC_Yue",
                                       "log2FC_10Mm_Laustriat", "log2FC_25Mm_Laustriat",
                                       "log2FC_lung_Xie", "log2FC_kidney_Xie")




pheatmap_output <- pheatmap(mat = all_4plot_4_series_filt_reorderd, color = test,
                            annotation_col = lib_index_no_series_type,
                            show_rownames = TRUE,
                            annotation_colors = colors_list_filt, border_color = "grey8",
                            breaks = c(-10,-4,-3,-2,-1, -0.5, 0, 0.5,1,2,3,4, 10),
                            cluster_cols = FALSE,
                            clustering_distance_rows = "canberra")

