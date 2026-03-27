## code to prepare `mydataset` dataset goes here
Left_Kidney <- load("data-raw/Left_Kidney.rda")
D_sim <- load("data-raw/simulated_data.rda")
usethis::use_data(Left_Kidney, D_sim, overwrite = TRUE)
