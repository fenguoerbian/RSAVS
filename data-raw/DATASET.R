## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)
por_df <- read.csv(file = "student-por.csv", sep = ";")
mat_df <- read.csv(file = "student-mat.csv", sep = ";")
full_df <- rbind(por_df, mat_df)