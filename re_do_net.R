library(ggraph)
library(tidygraph)

cor_all_in_heatmap_ffilt <-cor_all_in_heatmap %>% dplyr::select(2:7)

rownames(cor_all_in_heatmap_ffilt) <-cor_all_in_heatmap[,1]

class(cor_all_in_heatmap_ffilt$cor.y.y.y)

cor_all_in_heatmap_ffilt[is.na(cor_all_in_heatmap_ffilt)] <- 0

pos_cor_heatmap_partial <-cor_all_in_heatmap_ffilt %>% 
  dplyr::filter(cor.x >=0 & cor.y >=0 & cor.x.x >=0 &cor.y.y >= 0 & cor.x.x.x >=0 & cor.y.y.y>=0)

neg_cor_heatmap_partial <-cor_all_in_heatmap_ffilt %>% 
  dplyr::filter(cor.x <=0 & cor.y <=0 & cor.x.x <=0 &cor.y.y <= 0 & cor.x.x.x <=0 & cor.y.y.y<=0)

###

cor_all_in_heatmap_pos <-cor_all_in_heatmap %>% filter(pairs %in% rownames(pos_cor_heatmap_partial)) %>%
  mutate(cor_type = "pos_cor") %>% dplyr:: select(1,8,9,10)

cor_all_in_heatmap_neg <-cor_all_in_heatmap %>% filter(pairs %in% rownames(neg_cor_heatmap_partial)) %>%
  mutate(cor_type = "neg_cor") %>% dplyr:: select(1,8,9,10)

cor_all_nodic<- rbind(cor_all_in_heatmap_pos,cor_all_in_heatmap_neg) %>% dplyr::select(-1)

##

cor_all_in_heatmap_neg_dic<- cor_all_in_heatmap_neg %>% left_join(txtogenesv37, by = c("row"="transcript_name")) %>% 
  select(1:4, 9)

cor_all_in_heatmap_neg_dic<- cor_all_in_heatmap_neg_dic %>% left_join(txtogenesv37, by = c("col"="transcript_name")) %>% 
  select(1:5, 10)

##

cor_all_in_heatmap_pos_dic<- cor_all_in_heatmap_pos %>% left_join(txtogenesv37, by = c("row"="transcript_name")) %>% 
  select(1:4, 9)

cor_all_in_heatmap_pos_dic<- cor_all_in_heatmap_pos_dic %>% left_join(txtogenesv37, by = c("col"="transcript_name")) %>% 
  select(1:5, 10)


cor_all_in_heatmap_all_dic<-rbind(cor_all_in_heatmap_neg_dic, cor_all_in_heatmap_pos_dic)

#######################

library(ggraph)
library(tidygraph)

names(cor_all_nodic) <- c("From", "To", "cor_type")

cor_all_tbl_graph <- as_tbl_graph(cor_all_nodic)

cor_all_tbl_graph %>%
  ggraph() +
  geom_node_point() +
  geom_edge_link()


nodes_mod <- cor_all_tbl_graph %>% as_tibble() %>% 
  left_join(txtogenesv37, by = c("name"="transcript_name")) %>%
  dplyr::select(1,8)

edges_mod<-cor_all_tbl_graph %>% activate(edges) %>% as_tibble()

write.csv(nodes_mod, "tables_fornet/nodes_forgraph.csv")
write.csv(edges_mod, "tables_fornet/edges_forgraph.csv")

cor_all_tbl_graph %>%
  mutate(trans_type =nodes_mod$transcript_type) %>%
  ggraph(layout = 'kk') +
  geom_node_point() +
  geom_edge_link( alpha =0.4)

  
############################################

gene_targets <- c("PCNP-207", "HNRNPD-201", "HNRNPD-203", "HNRNPDL-206", "COQ2-206",
                  "MRPS18C-204", "SYNCRIP-216", "SYNCRIP-208", "DDX56-201", "POLR2J4-203",
                  "SYNCRIP-203", "MELK-201", "MELK-203", "EIF4G2-214", "EIF4G2-225",
                  "POLA2-201", "SF1-204", "EHBP1L1-201", "SIPA1-209", "MUS81-201",
                  "PELI3-201", "POLA2-201", "SF1-204", "EHBP1L1-201", "SIPA1-209",
                  "MUS81-201", "PELI3-201", "FOXM1-202", "DMAC2L-211", "NIN-204",
                  "PKMYT1-201", "HCFC1R1-201", "THOC6-202", "CREBBP-206", "PIGL-214",
                  "PIGL-214", "SRSF2-202", "H3-3B-206", "UNK-202", "EXOC7-216",
                  "POFUT2-201", "HAUS7-203")

lncRNAs_source <-c("ZBTB11-AS1-201","THAP9-AS1-204","THAP9-AS1-204","THAP9-AS1-204",
                   "THAP9-AS1-204","THAP9-AS1-204","SNHG5-201","SNHG5-201","SNHG15-214",
                   "SNHG15-214","SNHG5-201","EBLN3P-204","EBLN3P-204","ZBED5-AS1-207",
                   "ZBED5-AS1-207","MALAT1-215","MALAT1-215","MALAT1-215","MALAT1-215",
                   "MALAT1-215","MALAT1-215","MALAT1-215","MALAT1-215","MALAT1-215",
                   "MALAT1-215","MALAT1-215","MALAT1-215","AC125807.2-201","LINC01588-205",
                   "LINC01588-205","ZNF213-AS1-203","ZNF213-AS1-203","ZNF213-AS1-203",
                   "ZNF213-AS1-203","SNHG29-243","SNHG29-209","LINC00511-243","LINC00511-243",
                   "LINC00511-243","LINC00511-243","MCM3AP-AS1-214","CH17-340M24.3-201")
####################################################
df_fornet <-data.frame(gene_targets,lncRNAs_source)

cor_net_fnal  <- as_tbl_graph(df_fornet)

cor_net_fnal %>%
  ggraph() +
  geom_node_point() +
  geom_edge_link()


nodes_mod <- cor_net_fnal %>% as_tibble() %>% 
  left_join(txtogenesv37, by = c("name"="transcript_name")) %>%
  dplyr::select(1,8)



cor_net_fnal %>%
  mutate(trans_type =nodes_mod$transcript_type) %>%
  ggraph(layout = 'sugiyama') +
  geom_node_point() +
  geom_edge_diagonal2(alpha =0.4)+
  geom_node_text(aes(label = name))+
  coord_flip()+
  theme_minimal()

edges_mod<-cor_net_fnal %>% activate(edges) %>% as_tibble()



########################





