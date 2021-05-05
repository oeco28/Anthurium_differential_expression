library(WGCNA);
library(plyr);
library("DESeq2");

stage.group <- rep(c("early","late"),times=3,each=3)
colour.group <- rep(c("white","red","purple"),times=1,each=6)

myGroups <- data.frame(cbind(stage.group,colour.group))
colnames(myGroups) <- c("stage","Colour")
rownames(myGroups) <- colnames(data_ec_2)
#rownames(datExpr0) <- colnames(data_ec_2)

select = order(rowMeans(data_ec_2), decreasing =T)[1:30000]
x.common = data_ec_2[select,]
dim(x.common)
x.round = round(x.common)
ddConsensus = DESeqDataSetFromMatrix(countData = x.round, colData = myGroups, design = ~Colour)

vsdConsensus <- varianceStabilizingTransformation(ddConsensus, blind=T)
transformed.consensus= assay(vsdConsensus)
datExpr = t(assay(vsdConsensus))
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, blockSize=5000)



par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
abline(h=0.80,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

################################################

adjacencyAll = adjacency(datExpr,power=6,type="signed");
diag(adjacencyAll)=0
dissTOMAll   = 1-TOMsimilarity(adjacencyAll, TOMType="signed")
geneTreeAll  = hclust(as.dist(dissTOMAll), method="average")
collectGarbage()


stage.group.num <- rep(c(1,2),times=3,each=3)
colour.group.num <- rep(c(1,2,3),times=1,each=6)

datTrait <- data.frame(cbind(stage.group.num,colour.group.num))
colnames(datTrait) <- c("stage","Colour")
rownames(datTrait) <- colnames(data_ec_2)



################################
#  Another plot
##################################

plot(geneTreeAll,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (combined dataset)", labels=FALSE,hang=0.04);

# Plot relationship among samples 
sampleTree = hclust(dist(datExpr), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Add trait heatmap to sample tree plot 
plotDendroAndColors(sampleTree, datTraits,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

# Set color mapping to NULL
mColorh=NULL


# Run cuttreeHybrid for four different minimum cluster sizes and deep split values 
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTreeAll, pamStage=FALSE,
                      minClusterSize = (30-3*ds), cutHeight = 0.99, 
                      deepSplit = ds, distM = dissTOMAll)
  mColorh=cbind(mColorh,labels2colors(tree$labels));
}

# Plot gene tree with different split values and clusters plotted below
plotDendroAndColors(geneTreeAll, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE);

# Assign color to each value in trait matrix
traitColors = numbers2colors(datTrait, signed = FALSE);

# Choose cutoff of first split (ds=0)
moduleColors = mColorh[,1]
modulesAll = mColorh[,1]

# Plot gene tree with chosen module 
plotDendroAndColors(geneTreeAll, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")


#############################################

PCs1A    = moduleEigengenes(datExpr,  colors=modulesAll) 
ME_1A    = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a") 
MDS_1A   = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(modulesAll))


#####################################
# It is time to plot more
#######################################


# Plot Module Eigengene values, the below is best plotted into a pdf all together
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="",sub="")
plot(MDS_1A, col= colorsA1,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTreeAll$order
plotMat(scale(log(t(datExpr)[ordergenes,])) , rlabels= modulesAll[ordergenes], clabels= colnames(t(datExpr)), rcols=modulesAll[ordergenes])

for (which.module in names(table(modulesAll))){
  ME = ME_1A[, paste("ME",which.module, sep="")] 
  barplot(ME, col=which.module, main="", cex.main=2, 
          ylab="eigengene expression",xlab="array sample") 
}; 
# End of plot function would be here
# dev.off();


#=====================================================================================
#
#  Association with trait
#
#=====================================================================================

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
datTrait = data.frame(cbind(stage.group,colour.group))
moduleTraitCor = cor(MEs, datTrait, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTrait),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# All below is pulling out information about the different modules. 
modNames = substring(names(MEs), 3)
module = "magenta"
column = match(module, modNames);
moduleGenes = moduleColors==module;


# Define variable weight containing the weight column of datTrait
flower.colour = as.data.frame(datTrait$colour.group);
names(flower.colour) = "white-red-purple"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, flower.colour, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(flower.colour), sep="");
names(GSPvalue) = paste("p.GS.", names(flower.colour), sep="");


sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for fresh-sulfidic",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
colnames(datExpr)
colnames(datExpr)[moduleColors=="magenta"]

magentaNames = colnames(datExpr)[moduleColors=="magenta"]
greenNames = colnames(datExpr)[moduleColors=="green"]

