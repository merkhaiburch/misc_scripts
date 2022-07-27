# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-09-21 
# Updated... 2022-05-26
# Description 
#   - Make a QR code using R with any link (for free!) 
# ---------------------------------------------------------------

install.packages("qrcode")
library(qrcode)

# Maize genetics website link
png("maize_genetics_MKB.png")
qr_code("https://www.maizegenetics.net/merrittkhaipho-burch")
dev.off()

# rick roll
png("~/Downloads/rick.png")
qr_code("https://youtu.be/dQw4w9WgXcQ")
dev.off()

# Google scholar profile
png("~/Box Sync/Cornell_PhD/labProjects/presentations/2022_PEQG/google_scholar_mbkb.png")
qr <- qrcode::qr_code("https://scholar.google.com/citations?user=lvFMbpwAAAAJ&hl=en")
plot(qr)
dev.off()