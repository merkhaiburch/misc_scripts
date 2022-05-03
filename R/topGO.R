# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-04-07 
# Updated... 2022-04-07

# Description 
# Run topGO on the example data set using quantitative values
# ---------------------------------------------------------------

# Load in packages
library(topGO)
library(ALL)
library(hgu95av2.db)

# Load required data
data(ALL)
data(geneList)

# microarray used in the experiment is the hgu95av2 from Affymetrix, as we can see from the affyLib object.
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)


# -------------------------------------
# 3.1. Data preparation
# ------------------------------------

# count the number of genes with a p-value less than 0.01
# there are 50 genes with a raw p-value less than 0.01 out of a total of 323 genes
sum(topDiffGenes(geneList))

# build topGOdata object for Biological process terms using all 323 genes
sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, #calling all 323 genes
                    geneSel = topDiffGenes, # select the 50 genes with p-value less than 0.01
                    nodeSize = 10, # only select GO terms with 10 observations
                    annot = annFUN.db, # used to extract the gene-to-GO mappings from the affyLib object
                    affyLib = affyLib) # human genome GO terms

# look at summary of 
sampleGOdata


# -------------------------------------
# 3.2 Performing the enrichment tests
# -------------------------------------

#  perform a classical enrichment analysis by testing the over-representation 
# of GO terms within the group of differentially expressed genes. For the 
# method classic each GO category is tested independently.
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

# runTest returns an object of class topGOresult.
resultFisher

# test the enrichment using the Kolmogorov-Smirnov test. We will use the 
# both the classic and the elim method.
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")


# -------------------------------------
# 3.3. Analysis of results
# -------------------------------------

# GenTable is an easy to use function for analysing the most significant 
# GO terms and the corresponding pvalues. In the following example, we 
# list the top 10 significant GO terms identified by the elim method. At
# the same time we also compare the ranks and the p-values of these GO 
# terms with the ones obtained by the classic method.
# GenTable function returns a data frame containing the top topNodes GO terms identified by the elim algorithm
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
allRes$temp <- rep("feature", nrow(allRes))
allRes$prop <- allRes$Significant/allRes$Annotated
allRes$expression <- rnorm(nrow(allRes))

# Plot results using ggplot2
ggplot(allRes, aes(x = temp, y = GO.ID, size = prop, color = as.numeric(as.character(classicFisher)))) + 
  geom_point() +
  scale_size_area()

# make heatmap
ggplot(allRes, aes(x = temp, y = GO.ID, fill= expression)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") 


# Access GO terms p-values using score function
pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
gCol <- colMap(gstat$Significant)

# Plot differnces in two methods
plot(pValue.classic, pValue.elim, 
     xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)

# look at how significant GO terms are distributed over the GO graph
# significant nodes are rectangles
library(Rgraphviz)
showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')


# -------------------------------------
# Loading genes and annotations data
# -------------------------------------

# get all biological process terms
BPterms <- ls(GOBPTerm)
head(BPterms)

# -------------------------------------
# Using the genes score
# -------------------------------------

y <- as.integer(sapply(eset$BT, function(x) return(substr(x, 1, 1) == 'T')))
table(y)

geneList <- getPvalues(exprs(eset), classlabel = y, alternative = "greater")

topDiffGenes <- function(allScore) {
  return(allScore < 0.01)
}
x <- topDiffGenes(geneList)
sum(x) ## the number of selected genes

# build top go object
GOdata <- new("topGOdata",
              description = "GO analysis of ALL data; B-cell vs T-cell",
              ontology = "BP",
              allGenes = geneList,
              geneSel = topDiffGenes,
              annot = annFUN.db,
              nodeSize = 5, # nodeSize argument the user can control the size of the GO terms used in the analysis.
              affyLib = affyLib)