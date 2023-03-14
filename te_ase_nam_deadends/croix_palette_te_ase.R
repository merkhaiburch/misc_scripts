if (!require("devtools")) install.packages("devtools")
devtools::install_github("btmonier/croix")


# 
crPal <- croix::croix_palette(name = "mov_edward_scissorhands")
plot(crPal)

# continious
crPal <- croix::croix_palette(
  name = "mov_edward_scissorhands", 
  n = 29, 
  type = "continuous"
)
plot(crPal)


crPal <- croix::croix_palette(
  name = "mov_edward_scissorhands", 
  n = 5, 
  type = "continuous"
)
plot(crPal)
$