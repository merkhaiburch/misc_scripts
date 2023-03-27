# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-02-20
# Updated... 2023-02-20
#
# Description:
# Make QR codes for maize meeting poster
# ------------------------------------------------------------------------------

# Load package
library(qrcode)

# Maize genetics website
png("~/Downloads/maizeGeneticsGame.png")
qr <- qrcode::qr_code("https://www.maizegenetics.net/game")
plot(qr)
dev.off()

# Bitbucket repo with all of the materials
png("~/Downloads/bitbucketGame.png")
qr <- qrcode::qr_code("https://bitbucket.org/bucklerlab/seed_dispersal_game/src/master/")
plot(qr)
dev.off()
