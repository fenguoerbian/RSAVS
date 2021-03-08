por_df <- read.csv(file = "data-raw/student-por.csv", sep = ";")
mat_df <- read.csv(file = "data-raw/student-mat.csv", sep = ";")
full_df <- rbind(por_df, mat_df)

usethis::use_data(por_df, mat_df, full_df, overwrite = TRUE)