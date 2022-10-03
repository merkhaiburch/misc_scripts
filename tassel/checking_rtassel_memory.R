# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-09- 
# Updated... 2022-09-
#
# Description:
# 
# ------------------------------------------------------------------------------

# Set memory options
options(java.parameters = "-Xmx500g")
library(rJava)
.jinit()
memUsageAfter <- rJava::.jcall(rJava::.jnew("java/lang/Runtime"), "J", "maxMemory") / 1024^3
message("Your max memory usage after changing it is: ", memUsageAfter, "GB")

# Check java jvm
rJava::.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
