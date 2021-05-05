library(devtools)
#install_github('sinhrks/ggfortify')
library(ggfortify); library(ggplot2); library(edgeR); library(gplots); library(LDheatmap); library(RColorBrewer)

setwd("/media/venezuela3/Anthurium/Adapter_Trimmed_fastq")
purple_early_1 <- data.frame(read.table("Purple_Early_1/rsem/3_run/pe1/pe1.genes.results",header=TRUE))
purple_early_2 <- data.frame(read.table("Purple_Early_2/rsem/3_run/pe2/pe2.genes.results",header=TRUE))
purple_early_3 <- data.frame(read.table("Purple_Early_3/rsem/3_run/pe3/pe3.genes.results",header=TRUE))
purple_late_1 <- data.frame(read.table("Purple_Late_1/rsem/3_run/pl1/pl1.genes.results",header=TRUE))
purple_late_2 <- data.frame(read.table("Purple_Late_2/rsem/3_run/pl2/pl2.genes.results",header=TRUE))
purple_late_3 <- data.frame(read.table("Purple_Late_3/rsem/3_run/pl3/pl3.genes.results",header=TRUE))
red_early_1 <- data.frame(read.table("Red_Early_1/rsem/3_run/re1/re1.genes.results",header=TRUE))
red_early_2 <- data.frame(read.table("Red_Early_2/rsem/3_run/re2/re2.genes.results",header=TRUE))
red_early_3 <- data.frame(read.table("Red_Early_3/rsem/3_run/re3/re3.genes.results",header=TRUE))
red_late_1 <- data.frame(read.table("Red_Late_1/rsem/3_run/rl1/rl1.genes.results",header=TRUE))
red_late_2 <- data.frame(read.table("Red_Late_2/rsem/3_run/rl2/rl2.genes.results",header=TRUE))
red_late_3 <- data.frame(read.table("Red_Late_3/rsem/3_run/rl3/rl3.genes.results",header=TRUE))
white_early_1 <- data.frame(read.table("White_Early_1/rsem/3_run/we1/we1.genes.results",header=TRUE))
white_early_2 <- data.frame(read.table("White_Early_2/rsem/3_run/we2/we2.genes.results",header=TRUE))
white_early_3 <- data.frame(read.table("White_Early_3/rsem/3_run/we3/we3.genes.results",header=TRUE))
white_late_1 <- data.frame(read.table("White_Late_1/rsem/3_run/wl1/wl1.genes.results",header=TRUE))
white_late_2 <- data.frame(read.table("White_Late_2/rsem/3_run/wl2/wl2.genes.results",header=TRUE))
white_late_3 <- data.frame(read.table("White_Late_3/rsem/3_run/wl3/wl3.genes.results",header=TRUE))

purple_early_1 <- purple_early_1[order(purple_early_1$gene_id),]
purple_early_2 <- purple_early_2[order(purple_early_2$gene_id),]
purple_early_3 <- purple_early_3[order(purple_early_3$gene_id),]
purple_late_1 <- purple_late_1[order(purple_late_1$gene_id),]
purple_late_2 <- purple_late_2[order(purple_late_2$gene_id),]
purple_late_3 <- purple_late_3[order(purple_late_3$gene_id),]
red_early_1 <- red_early_1[order(red_early_1$gene_id),]
red_early_2 <- red_early_2[order(red_early_2$gene_id),]
red_early_3 <- red_early_3[order(red_early_3$gene_id),]
red_late_1 <- red_late_1[order(red_late_1$gene_id),]
red_late_2 <- red_late_2[order(red_late_2$gene_id),]
red_late_3 <- red_late_3[order(red_late_3$gene_id),]
white_early_1 <- white_early_1[order(white_early_1$gene_id),]
white_early_2 <- white_early_2[order(white_early_2$gene_id),]
white_early_3 <- white_early_3[order(white_early_3$gene_id),]
white_late_1 <- white_late_1[order(white_late_1$gene_id),]
white_late_2 <- white_late_2[order(white_late_2$gene_id),]
white_late_3 <- white_late_3[order(white_late_3$gene_id),]

data_ec <- cbind(white_early_1$expected_count,white_early_2$expected_count,white_early_3$expected_count,white_late_1$expected_count,white_late_2$expected_count,white_late_3$expected_count,red_early_1$expected_count,red_early_2$expected_count,red_early_3$expected_count,red_late_1$expected_count,red_late_2$expected_count,red_late_3$expected_count,purple_early_1$expected_count,purple_early_2$expected_count,purple_early_3$expected_count,purple_late_1$expected_count,purple_late_2$expected_count,purple_late_3$expected_count);
data_FPKM <- cbind(white_early_1$FPKM,white_early_2$FPKM,white_early_3$FPKM,white_late_1$FPKM,white_late_2$FPKM,white_late_3$FPKM,red_early_1$FPKM,red_early_2$FPKM,red_early_3$FPKM,red_late_1$FPKM,red_late_2$FPKM,red_late_3$FPKM,purple_early_1$FPKM,purple_early_2$FPKM,purple_early_3$FPKM,purple_late_1$FPKM,purple_late_2$FPKM,purple_late_3$FPKM);
data_TPM <- cbind(white_early_1$TPM,white_early_2$TPM,white_early_3$TPM,white_late_1$TPM,white_late_2$TPM,white_late_3$TPM,red_early_1$TPM,red_early_2$TPM,red_early_3$TPM,red_late_1$TPM,red_late_2$TPM,red_late_3$TPM,purple_early_1$TPM,purple_early_2$TPM,purple_early_3$TPM,purple_late_1$TPM,purple_late_2$TPM,purple_late_3$TPM);

colnames(data_ec) <- c("we1","we2","we3","wl1","wl2","wl3","re1","re2","re3","rl1","rl2","rl3","pe1","pe2","pe3","pl1","pl2","pl3")
colnames(data_FPKM) <- c("we1","we2","we3","wl1","wl2","wl3","re1","re2","re3","rl1","rl2","rl3","pe1","pe2","pe3","pl1","pl2","pl3")
colnames(data_TPM) <- c("we1","we2","we3","wl1","wl2","wl3","re1","re2","re3","rl1","rl2","rl3","pe1","pe2","pe3","pl1","pl2","pl3")
rownames(data_ec) <- white_late_3$gene_id
rownames(data_FPKM) <- white_late_3$gene_id
rownames(data_TPM) <- white_late_3$gene_id

data_TPM_2 <- data_TPM[rowSums(data_TPM) >= 18,]  # Exclude genes where most expression is below 1 TPM
data_FPKM_2 <- data_FPKM[rowSums(data_TPM) >= 18,] # Exclude genes where most expression is below 1 TPM
data_ec_2 <- data_ec[rowSums(data_TPM) >= 18,] # Exclude genes where most expression is below 1 TPM

data_ec_3 <- data_ec_2[,-8]
targets2 <- data.frame(matrix(NA,17,2))
targets2[,1] <- c(rep("early",3),rep("late",3),rep("early",2),rep("late",3),rep("early",3),rep("late",3))
targets2[,2] <- c(rep("White",6),rep("Red",5),rep("Purple",6))
rownames(targets2) <- c("we1","we2","we3","wl1","wl2","wl3","re1","re3","rl1","rl2","rl3","pe1","pe2","pe3","pl1","pl2","pl3")
colnames(targets2) <- c("stage","color")
Group2 <- factor(paste(targets2$stage,targets2$color,sep="."))
samples2 <- c("we1","we2","we3","wl1","wl2","wl3","re1","re3","rl1","rl2","rl3","pe1","pe2","pe3","pl1","pl2","pl3")
samples2 <- data.frame(cbind(samples2, as.character(Group2)))
colnames(samples2)<- c("samples","group")
m2 = DGEList(counts = data_ec_3, group = samples2$Group2)
m2 = calcNormFactors(m2)
design2 <- model.matrix(~0+Group2)
colnames(design2) <- levels(Group2)
m2 <- estimateCommonDisp(m2)
m2 <- estimateGLMTrendedDisp(m2, design2)
m2 <- estimateGLMTagwiseDisp(m2, design2)
m2plot <- calcNormFactors(m2)

mds <- plotMDS(m2plot,top=10000,gene.selection="common")
mds2 <- plotMDS(mds, top=10000,gene.selection="common", col=c(rep(colors()[342],6), rep(colors()[35],5), rep(colors()[99],6)), labels= c(rep("White",6), rep("Red",5), rep("Purple",6)))
mds2 <- plotMDS(mds, top=10000,gene.selection="common", col=c(rep(colors()[342],6), rep(colors()[35],5), rep(colors()[99],6)), pch=c(rep(16,3),rep(17,3),rep(16,2),rep(17,3),rep(16,3),rep(17,3)), cex=2.5)
pdf("working_analysis2/MDS_plot.pdf")
mds2 <- plotMDS(mds, top=10000,gene.selection="common", col=c(rep(colors()[342],6), rep(colors()[35],5), rep(colors()[99],6)), pch=c(rep(16,3),rep(17,3),rep(16,2),rep(17,3),rep(16,3),rep(17,3)), cex=2.5)
dev.off()

fit <- glmQLFit(m2,design2)
Stage <- c(rep("early",3),rep("late",3),rep("early",2),rep("late",3),rep("early",3),rep("late",3))
Color <- c(rep("White",6),rep("Red",5),rep("Purple",6))
#A <- cbind(mds2$x,mds2$y,Stage,Color)


my.contrasts <- makeContrasts(
  stage_white = late.White - early.White,
  stage_red = late.Red - early.Red,
  stage_purple = late.Purple - early.Purple,
  white_red = late.Red - late.White,
  white_purple = late.Purple - late.White,
  purple_red = late.Purple - late.Red,
  ewhite_red = early.Red - early.White,
  ewhite_purple = early.Purple - early.White,
  epurple_red = early.Purple - early.Red,
  levels = design2
)

SW_qlrt = glmQLFTest(fit, contrast = my.contrasts[,"stage_white"])
SR_qlrt = glmQLFTest(fit, contrast = my.contrasts[,"stage_red"])
SP_qlrt = glmQLFTest(fit, contrast = my.contrasts[,"stage_purple"])
WR_qlrt = glmQLFTest(fit, contrast = my.contrasts[,"white_red"])
WP_qlrt = glmQLFTest(fit, contrast = my.contrasts[,"white_purple"])
RP_qlrt = glmQLFTest(fit, contrast = my.contrasts[,"purple_red"])
eWR_qlrt = glmQLFTest(fit, contrast = my.contrasts[,"ewhite_red"])
eWP_qlrt = glmQLFTest(fit, contrast = my.contrasts[,"ewhite_purple"])
eRP_qlrt = glmQLFTest(fit, contrast = my.contrasts[,"epurple_red"])

qlSW = topTags(SW_qlrt, n = 50000)[[1]]
qlSR = topTags(SR_qlrt, n = 50000)[[1]]
qlSP = topTags(SP_qlrt, n = 50000)[[1]]
qlWR = topTags(WR_qlrt, n = 50000)[[1]]
qlWP = topTags(WP_qlrt, n = 50000)[[1]]
qlRP = topTags(RP_qlrt, n = 50000)[[1]]
qleWR = topTags(eWR_qlrt, n = 50000)[[1]]
qleWP = topTags(eWP_qlrt, n = 50000)[[1]]
qleRP = topTags(eRP_qlrt, n = 50000)[[1]]

qlSW.sub <- qlSW[which(qlSW$FDR <= 0.0005 & abs(qlSW$logFC) >= 2),] 
qlSR.sub <- qlSR[which(qlSR$FDR <= 0.0005 & abs(qlSW$logFC) >= 2),] 
qlSP.sub <- qlSP[which(qlSP$FDR <= 0.0005 & abs(qlSW$logFC) >= 2),] 
qlWR.sub <- qlWR[which(qlWR$FDR <= 0.0005 & abs(qlSW$logFC) >= 2),] 
qlWP.sub <- qlWP[which(qlWP$FDR <= 0.0005 & abs(qlSW$logFC) >= 2),] 
qlRP.sub <- qlRP[which(qlRP$FDR <= 0.0005 & abs(qlSW$logFC) >= 2),] 
qleWR.sub <- qleWR[which(qleWR$FDR <= 0.0005 & abs(qlSW$logFC) >= 2),]
qleWP.sub <- qleWP[which(qleWP$FDR <= 0.0005 & abs(qlSW$logFC) >= 2),]
qleRP.sub <- qleRP[which(qleRP$FDR <= 0.0005 & abs(qlSW$logFC) >= 2),]

write.csv(qlSW, file = "working_analysis2021/sw.csv")
write.csv(qlSR, file = "working_analysis2021/sr.csv")
write.csv(qlSP, file = "working_analysis2021/sp.csv")
write.csv(qlWR, file = "working_analysis2021/late_wr.csv")
write.csv(qlWP, file = "working_analysis2021/late_wp.csv")
write.csv(qlRP, file = "working_analysis2021/late_rp.csv")
write.csv(qleWR, file = "working_analysis2021/early_wr.csv")
write.csv(qleWP, file = "working_analysis2021/early_wp.csv")
write.csv(qleRP, file = "working_analysis2021/early_rp.csv")

###########################################
#################
############################################

#library(gplots)
#library(LDheatmap)
#library(RColorBrewer)

data_FPKM_3 <- data_FPKM_2[,-8]

#heatmap.2(data_FPKM_3[rownames(qlSW.sub),c(1:6)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( rev(brewer.pal(9, "Spectral")) )(255), labRow="", Colv=FALSE, main="Differentially expressed genes SW")
pdf("working_analysis2021/SW_heatmap.pdf",width=11,height=9)
heatmap.2(data_FPKM_3[rownames(qlSW.sub),c(1:6)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk1","cornsilk2","cornsilk","cornsilk3","cornsilk4"))(100), labRow="", Colv=FALSE, main="DE genes stage White") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()

pdf("working_analysis2021/SW_heatmap2.pdf",width=11,height=9)
heatmap.2(data_ec_3[rownames(qlSW.sub),c(1:6)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk1","cornsilk2","cornsilk","cornsilk3","cornsilk4"))(100), labRow="", Colv=FALSE, main="DE genes stage White") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()
#heatmap.2(data_FPKM_3[rownames(SW),c(1:6)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk1","cornsilk2","cornsilk","cornsilk3","cornsilk4"))(100), labRow="", Colv=FALSE, main="DE genes stage White")

#heatmap.2(data_FPKM_3[rownames(qlSR.sub),c(7:11)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), labRow="", Colv=FALSE, main="Differentially expressed genes SR")
pdf("working_analysis2021/SR_heatmap.pdf",width=11,height=9)
heatmap.2(data_FPKM_3[rownames(qlSR.sub),c(7:11)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("pink","indianred1","firebrick2","firebrick3","firebrick4"))(100), labRow="", Colv=FALSE, main="DE genes stage Red") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()

pdf("working_analysis2021/SR_heatmap2.pdf",width=11,height=9)
heatmap.2(data_ec_3[rownames(qlSR.sub),c(7:11)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("pink","indianred1","firebrick2","firebrick3","firebrick4"))(100), labRow="", Colv=FALSE, main="DE genes stage Red") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()
#heatmap.2(data_FPKM_3[rownames(SR),c(7:11)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("pink","indianred1","firebrick2","firebrick3","firebrick4"))(100), labRow="", Colv=FALSE, main="DE genes stage Red")

#heatmap.2(data_FPKM_3[rownames(qlSP.sub),c(12:17)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), labRow="", Colv=FALSE, main="Differentially expressed genes SP")
pdf("working_analysis2021/SP_heatmap.pdf",width=11,height=9)
heatmap.2(data_FPKM_3[rownames(qlSP.sub),c(12:17)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("thistle1","thistle2","slateblue","slateblue3","slateblue4"))(255), labRow="", Colv=FALSE, main="DE genes stage Purple") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()

pdf("working_analysis2021/SP_heatmap2.pdf",width=11,height=9)
heatmap.2(data_ec_3[rownames(qlSP.sub),c(12:17)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("thistle1","thistle2","slateblue","slateblue3","slateblue4"))(255), labRow="", Colv=FALSE, main="DE genes stage Purple") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()
#heatmap.2(data_FPKM_3[rownames(SP),c(12:17)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("thistle1","thistle2","slateblue","slateblue3","slateblue4"))(255), labRow="", Colv=FALSE, main="DE genes stage Purple")


#heatmap.2(data_FPKM_3[rownames(qlWR.sub),c(4:6,10:12)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), labRow="", Colv=FALSE, main="Differentially expressed genes WR")
pdf("working_analysis2021/WR_heatmap.pdf",width=11,height=9)
heatmap.2(data_FPKM_3[rownames(qlWR.sub),c(4:6,9:11)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","cornsilk","indianred1","firebrick3","firebrick4"))(100), labRow="", Colv=FALSE, main="DE genes late WR") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()

pdf("working_analysis2021/WR_heatmap2.pdf",width=11,height=9)
heatmap.2(data_ec_3[rownames(qlWR.sub),c(4:6,9:11)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","cornsilk","indianred1","firebrick3","firebrick4"))(100), labRow="", Colv=FALSE, main="DE genes late WR") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()
#heatmap.2(data_FPKM_3[rownames(WR),c(4:6,9:11)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","cornsilk","indianred1","firebrick3","firebrick4"))(100), labRow="", Colv=FALSE, main="DE genes late WR")


#heatmap.2(data_FPKM_3[rownames(qlWP.sub),c(4:6,16:18)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), labRow="", Colv=FALSE, main="Differentially expressed genes WP")
pdf("working_analysis2021/WP_heatmap.pdf",width=11,height=9)
heatmap.2(data_FPKM_3[rownames(qlWP.sub),c(4:6,15:17)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","thistle1","slateblue3","slateblue4"))(100), labRow="", Colv=FALSE, main="DE genes late WP") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()

pdf("working_analysis2021/WP_heatmap2.pdf",width=11,height=9)
heatmap.2(data_ec_3[rownames(qlWP.sub),c(4:6,15:17)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","thistle1","slateblue3","slateblue4"))(100), labRow="", Colv=FALSE, main="DE genes late WP") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()
#heatmap.2(data_FPKM_3[rownames(WP),c(4:6,15:17)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","thistle1","slateblue3","slateblue4"))(100), labRow="", Colv=FALSE, main="DE genes late WP")

#heatmap.2(data_FPKM_3[rownames(qlRP.sub),c(10:12,16:18)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), labRow="", Colv=FALSE, main="Differentially expressed genes RP")
pdf("working_analysis2021/RP_heatmap.pdf",width=11,height=9)
heatmap.2(data_FPKM_3[rownames(qlRP.sub),c(9:11,15:17)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( c("firebrick4","firebrick3","firebrick2","indianred1","pink","thistle1","thistle2","slateblue","slateblue3","slateblue4") )(100), labRow="", Colv=FALSE, main="DE genes late RP") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()

pdf("working_analysis2021/RP_heatmap2.pdf",width=11,height=9)
heatmap.2(data_ec_3[rownames(qlRP.sub),c(9:11,15:17)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( c("firebrick4","firebrick3","firebrick2","indianred1","pink","thistle1","thistle2","slateblue","slateblue3","slateblue4") )(100), labRow="", Colv=FALSE, main="DE genes late RP") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()
#heatmap.2(data_FPKM_3[rownames(RP),c(9:11,15:17)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( c("firebrick4","firebrick3","firebrick2","indianred1","pink","thistle1","thistle2","slateblue","slateblue3","slateblue4") )(100), labRow="", Colv=FALSE, main="DE genes late RP")


#early_colors
pdf("working_analysis2021/eWR_heatmap.pdf",width=11,height=9)
heatmap.2(data_FPKM_3[rownames(qleWR.sub),c(1:3,7:8)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","cornsilk","indianred1","firebrick3","firebrick4"))(100), labRow="", Colv=FALSE, main="DE genes early WR") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()

pdf("working_analysis2021/eWR_heatmap2.pdf",width=11,height=9)
heatmap.2(data_ec_3[rownames(qleWR.sub),c(1:3,7:8)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","cornsilk","indianred1","firebrick3","firebrick4"))(100), labRow="", Colv=FALSE, main="DE genes early WR") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()
#heatmap.2(data_FPKM_3[rownames(eWR),c(1:3,7:8)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","cornsilk","indianred1","firebrick3","firebrick4"))(100), labRow="", Colv=FALSE, main="DE genes early WR")

pdf("working_analysis2021/eWP_heatmap.pdf",width=11,height=9)
heatmap.2(data_FPKM_3[rownames(qleWP.sub),c(1:3,12:14)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","thistle1","slateblue3","slateblue4"))(100), labRow="", Colv=FALSE, main="DE genes early WP") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()

pdf("working_analysis2021/eWP_heatmap2.pdf",width=11,height=9)
heatmap.2(data_ec_3[rownames(qleWP.sub),c(1:3,12:14)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","thistle1","slateblue3","slateblue4"))(100), labRow="", Colv=FALSE, main="DE genes early WP") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()
#heatmap.2(data_FPKM_3[rownames(eWP),c(1:3,12:14)],trace="none", dendrogram="row",scale="row",col = colorRampPalette(c("cornsilk4","cornsilk2","thistle1","slateblue3","slateblue4"))(100), labRow="", Colv=FALSE, main="DE genes early WP")

pdf("working_analysis2021/eRP_heatmap.pdf",width=11,height=9)
heatmap.2(data_FPKM_3[rownames(qleRP.sub),c(7:8,12:14)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( c("firebrick4","firebrick3","firebrick2","indianred1","pink","thistle1","thistle2","slateblue","slateblue3","slateblue4") )(100), labRow="", Colv=FALSE, main="DE genes early RP") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()

pdf("working_analysis2021/eRP_heatmap2.pdf",width=11,height=9)
heatmap.2(data_ec_3[rownames(qleRP.sub),c(7:8,12:14)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( c("firebrick4","firebrick3","firebrick2","indianred1","pink","thistle1","thistle2","slateblue","slateblue3","slateblue4") )(100), labRow="", Colv=FALSE, main="DE genes early RP") + theme(axis.text.x= element_text(size = rel(1.4), angle=90)) 
dev.off()

#heatmap.2(data_FPKM_3[rownames(eRP),c(7:8,12:14)],trace="none", dendrogram="row",scale="row",col = colorRampPalette( c("firebrick4","firebrick3","firebrick2","indianred1","pink","thistle1","thistle2","slateblue","slateblue3","slateblue4") )(100), labRow="", Colv=FALSE, main="DE genes early RP")

#SW_t = topTags(SW_lrt, n = nrow(m))
#SR_t = topTags(SR_lrt, n = nrow(m))
#SP_t = topTags(SP_lrt, n = nrow(m))
#WR_t = topTags(WR_lrt, n = nrow(m))
#WP_t = topTags(WP_lrt, n = nrow(m))
#RP_t = topTags(RP_lrt, n = nrow(m))

#WR_tail <- tail(WR_t$table, n=1000)
#WP_tail <- tail(WP_t$table, n=1000)
#RP_tail <- tail(RP_t$table, n=1000)

#logFC_baseline <- c(WR_tail$logFC,WP_tail$logFC,RP_tail$logFC)
#logFC_baseline <- t(logFC_baseline)
#logFC_baseline <- t(logFC_baseline)
#logFC_baseline <- cbind(logFC_baseline,c(rep("WR",1000),rep("WP",1000),rep("RP",1000)))
#colnames(logFC_baseline) <- c("logFC","Comparison")

#logFC_baseline <- data.frame(logFC_baseline)
#library(ggplot2)
#ggplot(logFC_baseline, aes(x=Comparison, y=logFC)) + geom_boxplot()

library(reshape)
library(ggh4x)

A <- data_FPKM[c("c857_g857","c109497_g109497","c84141_g84141","c110775_g110775","c100209_g100209","c51578_g51578","c3479_g3479","c60106_g60106","c7734_g7734"),]
rownames(A) <- c("CHI","CHS","DFR","F3H","F3'5'H/F3'H","LDOX","UFGT","Cytochrome P450","MYB_domain")
A_cs <- A


for(i in 1:dim(A)[1]){
  for(j in 1:dim(A)[2]){
      A_cs[i,j] <- (A[i,j] - median(A[i,]))/sd(A[i,])
  }
  }


A_wr <- A_cs[,c(1:12)]
A_wr1 <- A_wr
for(i in 1:dim(A_wr)[1]){
  for(j in 1:dim(A_wr)[2]){
      A_wr1[i,j] <- (A_wr[i,j] - median(A_wr[i,]))/sd(A_wr[i,])
  }
  }

A_wp <- A_cs[,c(1:6,13:18)]
A_wp1 <- A_wp
for(i in 1:dim(A_wp)[1]){
  for(j in 1:dim(A_wp)[2]){
      A_wp1[i,j] <- (A_wp[i,j] - median(A_wp[i,]))/sd(A_wp[i,])
  }
  }

A_rp <- A_cs[,c(7:18)]
A_rp1 <- A_rp
for(i in 1:dim(A_rp)[1]){
  for(j in 1:dim(A_rp)[2]){
      A_rp1[i,j] <- (A_rp[i,j] - median(A_rp[i,]))/sd(A_rp[i,])
  }
  }


#B <- melt(data_FPKM_2[c("c857_g857","c109497_g109497","c84141_g84141","c110775_g110775","c100209_g100209","c110775_g110775","c3479_g3479"),])

B <- melt(A_cs)
colnames(B) <- c("gene","sample","FPKM")
B$treatment <- c(rep("we",27),rep("wl",27),rep("re",27),rep("rl",27),rep("pe",27),rep("pl",27))
B$treatment <- factor(B$treatment, levels=c("we","wl","re","rl","pe","pl"))

B1 <- B[which(B$gene == "c857_g857",]

pdf("dotplots4.3.pdf")
ggplot(data=B,aes(x=treatment, y=FPKM)) + geom_jitter(colour="darkred", size=2, width=0.1, height = 0.1) + theme(axis.text.x= element_text(size = rel(0.7), angle=90)) + scale_x_discrete(labels=c("we" = "WE","wl" = "WL","re" = "RE","rl" = "RL","pe" = "PE","pl" = "PL")) + ylab("Normalized expression") + labs(colour='Treatment') + facet_wrap(~ gene)
dev.off()




B1 <- cbind(melt(A_wr1),rep("WR",108))
B1$treatment <- c(rep("we",27),rep("wl",27),rep("re",27),rep("rl",27))
B1$treatment <- factor(B1$treatment, levels=c("we","wl","re","rl"))
B2 <- cbind(melt(A_wp1),rep("WP",108))
B2$treatment <- c(rep("we",27),rep("wl",27),rep("pe",27),rep("pl",27))
B2$treatment <- factor(B2$treatment, levels=c("we","wl","pe","pl"))
B3 <- cbind(melt(A_rp1),rep("RP",108))
B3$treatment <- c(rep("re",27),rep("rl",27),rep("pe",27),rep("pl",27))
B3$treatment <- factor(B3$treatment, levels=c("re","rl","pe","pl"))


colnames(B1) <- c("gene","sample","FPKM","comparison","treatment")
colnames(B2) <- c("gene","sample","FPKM","comparison","treatment")
colnames(B3) <- c("gene","sample","FPKM","comparison","treatment")

B <- rbind(B1,B2,B3)
B$comparison <- factor(B$comparison, levels=c("WR","WP","RP"))


pdf("dotplots5.pdf")
ggplot(data=B,aes(x=treatment, y=FPKM)) + geom_jitter(colour="darkred", size=2, width=0.1, height = 0.1) + theme(axis.text.x= element_text(size = rel(0.7), angle=90)) + scale_x_discrete(labels=c("we" = "WE","wl" = "WL","re" = "RE","rl" = "RL","pe" = "PE","pl" = "PL")) + ylab("Normalized expression") + labs(colour='Treatment') + facet_nested(gene ~ comparison) 
dev.off()

pdf("dotplots5.2.pdf")
ggplot(data=B,aes(x=treatment, y=FPKM)) + geom_jitter(colour="darkred", size=2, width=0.1, height = 0.1) + theme(axis.text.x= element_text(size = rel(0.7), angle=90)) + scale_x_discrete(labels=c("we" = "WE","wl" = "WL","re" = "RE","rl" = "RL","pe" = "PE","pl" = "PL")) + ylab("Normalized expression") + labs(colour='Treatment') + facet_nested(comparison ~ gene) 
dev.off()

pdf("dotplots5.3.pdf")
ggplot(data=B,aes(x=treatment, y=FPKM)) + geom_jitter(colour="darkred", size=2, width=0.1, height = 0.1) + theme(axis.text.x= element_text(size = rel(0.7), angle=90)) + scale_x_discrete(labels=c("we" = "WE","wl" = "WL","re" = "RE","rl" = "RL","pe" = "PE","pl" = "PL")) + ylab("Normalized expression") + labs(colour='Treatment') + facet_wrap(~ gene) 
dev.off()


