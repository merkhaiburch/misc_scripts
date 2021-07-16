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
gerpelem -f Zea_mays.allChr.rates -v


# log file
gpp-gerpelem
Thu Jul 15 13:24:38 EDT 2021
processing file Zea_mays.allChr.rates containing 398381696 position scores
median neutral rate is 5.41
0 non-border shallow positions excluded
terminate called after throwing an instance of 'std::bad_alloc'
  what():  std::bad_alloc
Aborted (core dumped)