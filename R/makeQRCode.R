# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-09-21 
#
# Description 
#   - Make a QR code using R with any link (for free!) 
# ---------------------------------------------------------------

install.packages("qrcode")
library(qrcode)

png("maize_genetics_MKB.png")
qrcode_gen("https://www.maizegenetics.net/merrittkhaipho-burch")
dev.off()
