# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-28
# Updated... 2023-06-28
#
# Description:
# Download Kotlin and run it
# ------------------------------------------------------------------------------

# Suggested to download into my bioinformatics folder
cd /home/mbb262/bioinformatics

wget https://github.com/JetBrains/kotlin/releases/download/v1.8.22/kotlin-compiler-1.8.22.zip

unzip kotlin-compiler-1.8.22.zip

export PATH=/workdir/mbb262/kotlin-compiler-1.8.22/bin:$PATH
export _JAVA_OPTIONS=-Xmx300g