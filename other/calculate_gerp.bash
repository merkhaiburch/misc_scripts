#!/bin/bash

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-07-15 
#
# Description 
#   - calculate GERP scores from Kistler data
#   - https://datadryad.org/stash/dataset/doi:10.5061/dryad.70t85k2
#   - readme: http://mendel.stanford.edu/SidowLab/downloads/gerp/Readme.txt
# ---------------------------------------------------------------

# Export path
export PATH=/programs/GERPplusplus:$PATH

# gerpcol already calculated by Kistler

# Get help
gerpelem -h

# Generate scores
gerpelem \
    -f Zea_mays.allChr.rates \
    -x kistler_gerp2 \
    -v
