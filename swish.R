library(tximeta)

library(fishpond)


files_df <- as.data.frame(files)

files_df <- rownames_to_column(files_df)

files_df <- left_join(files_df, xie_kidney_index, by = c("rowname"="X1"))

files_df <- files_df %>% select(1,2,4)

names(files_df) <- c("names", "files", "condition")

se <- tximeta(files_df)

se$condition <- factor(se$condition, levels = c("control","metformin"))

y<-se

y$condition
## Running Swish at the transcript level

# Squash Technical Reps
y <- fishpond::scaleInfReps(y)
# Filter rows
y <- fishpond::labelKeep(y, minCount = 10, minN = 3)
y <- y[mcols(y)$keep,]
# set.seed(1)
y <- fishpond::swish(y, x = "condition")


##############################

DET_results <- S4Vectors::mcols(y) %>%
  base::as.data.frame()  %>%
  tibble::rownames_to_column("transcript_id")  %>%
  tibble::as_tibble()


#####

test <- left_join(DET_results, txtogenesv37, by = "transcript_id")

Xie_kidney_swish<- test %>% dplyr::select(1,8,9,11,15,16,18)

write.csv(Xie_kidney_swish, "Xie_kidney_swish.csv")
##############

#lung swish

files_df <- as.data.frame(files)

files_df <- rownames_to_column(files_df)

files_df <- left_join(files_df, xie_lung_index, by = c("rowname"="X1"))

files_df <- files_df %>% dplyr::select(1,2,4)

names(files_df) <- c("names", "files", "condition")

se <- tximeta(files_df)

se$condition <- factor(se$condition, levels = c("control","metformin"))

y<-se

y$condition
## Running Swish at the transcript level

# Squash Technical Reps
y <- fishpond::scaleInfReps(y)
# Filter rows
y <- fishpond::labelKeep(y, minCount = 10, minN = 3)
y <- y[mcols(y)$keep,]
# set.seed(1)
y <- fishpond::swish(y, x = "condition")


##############################

DET_results <- S4Vectors::mcols(y) %>%
  base::as.data.frame()  %>%
  tibble::rownames_to_column("transcript_id")  %>%
  tibble::as_tibble()


#####

test <- left_join(DET_results, txtogenesv37, by = "transcript_id")

Xie_lung_swish<- test %>% dplyr::select(1,8,9,11,15,16,18)

write.csv(Xie_lung_swish, "Xie_lung_swish.csv")
