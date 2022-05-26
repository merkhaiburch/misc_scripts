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

png("maize_genetics_MKB.png")
qrcode_gen("https://www.maizegenetics.net/merrittkhaipho-burch")
dev.off()

# rick roll
png("~/Downloads/rick.png")
qr_code("https://youtu.be/dQw4w9WgXcQ")
dev.off()
