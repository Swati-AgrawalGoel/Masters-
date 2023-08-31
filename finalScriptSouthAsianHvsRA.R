library(minfi)
library(tibble)
library(dplyr)
library(EnhancedVolcano)
library("Gviz")
library("FlowSorted.Blood.EPIC")
library("RColorBrewer")
library(minfi)
library(wateRmelon)
library(ChAMP)
#BiocManager::install("DMRcate")
library("DMRcate")
#BiocManager::install("DMRcatedata")
library("DMRcatedata")
#BiocManager::install("stringr")
library("stringr")
library("mCSEA")
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library("mCSEA")
library("limma")
library(ggplot2)
library(karyoploteR)
library(MEAL)

# reading sample data file: just for checking
sam_target = read.csv("Desktop/meth_data/Sampleinfo_rasa.csv", 
                      strip.white = TRUE, stringsAsFactors = FALSE)

table(sam_target$Sample_Group) # to see number of sample gropus

# pointing R  towards the directory in whch data is saved.
datadir <- "Desktop/meth_data/"
# seeing the list of the list of all files: recursive open all the subfiles
list.files(datadir, recursive = TRUE)

#------------read in the sample sheet for the experiment----------
targets <- read.metharray.sheet(datadir, pattern="Sampleinfo_rasa.csv")
str(targets)

# --------subsetting the SA health and RA patients only.------------
tarHvR<- subset(targets, targets$Sample_Group == "SA_RA"|targets$Sample_Group == "SA_HC")
# arranging the sample group in an order
tarHvR <- tarHvR[order(tarHvR$Sample_Group),]

#------- making RGchannel object--------------
RGsetHvR <- read.metharray.exp(targets=tarHvR) # takes time be patient
RGsetHvR 

# -----------making RG channel extended object-----------------
RGsetHvR2 <- read.metharray.exp(targets=tarHvR, extended = TRUE)

# -----------1st filter  using watermelon package-----------------------
Filter_Mset <- pfilter(RGsetHvR2, pnthresh = 0.01)
#for p value = 0.01
#0 samples having 1 % of sites with a detection p-value greater than 0.01 were removed 
#Samples removed:  
# 8291 sites were removed as beadcount <3 in 5 % of samples 
#3520 sites having 1 % of samples with a detection p-value greater than 0.01 were removed 

#---------- checking for pdata and manifest files-------------------
# checking the metadata for Rgset1
pData(RGsetHvR)
# checking information regarging the probes used in EPIC array
getManifest(RGsetHvR)
# getting manifest
manifest <- getManifest(RGsetHvR)
head(getProbeInfo(manifest))
annotation(RGsetHvR) # "ilm10b4.hg19" 

#------------- methyl set from RG channel --------------
MSetHvR <- preprocessRaw(RGsetHvR)

# ----------getting ratio convert---------------
ratioSetHvR <- ratioConvert(MSetHvR , what = "both", keepCN = TRUE)
# ---------getting genomic ratio set by adding genomic coordinates to ratioset object----- 
gSetHvR <- mapToGenome(ratioSetHvR)
gSetHvR

#---------------sanity and quality control for getsex using genomic ratio set
#----------------and estimate sex from watermelon pkg--------------------------------
beta_val_gsetnoNorm <- getBeta(gSetHvR)
sex_pred <- estimateSex(beta_val_gsetnoNorm, do_plot=TRUE)

#------------getting beta and density plots from methy set 
phenoDataHvR <- pData(MSetHvR)
densityPlot(MSetHvR, sampGroups = tarHvR$Sample_Group)

# -----------quality control before normalization--------
qc <- getQC(MSetHvR)
plotQC(qc)
controlStripPlot(RGsetHvR, controls="BISULFITE CONVERSION II")
#------------- mds plot befor normalization-------------
mdsPlot(RGsetHvR, sampGroups =tarHvR$Sample_Group)
mdsPlot(RGsetHvR, sampGroups =tarHvR$Slide)


#----methyset from RG channel set extended object with p filter-----------
MSet_pft <- preprocessNoob(Filter_Mset)
MSet_pft
#---------------getting ratio convert for the methy set------------
ratioSetHvR2 <- ratioConvert(MSet_pft , what = "both", keepCN = TRUE)
# getting genomic ratio set by adding genomic coordinates to ratioset object 
gSetHvR2 <- mapToGenome(ratioSetHvR2)
gSetHvR2

#---------getting beta,copy number and M values
b_AfterfiltandNoob <- getBeta(gSetHvR2)
nrow(b_AfterfiltandNoob) # with pvalue <0.01 854565 cpg probes r there

#---------checking for outlier again using watermelon package
#------not needed. Needed only for sanity check---------------
#outliers2 <- outlyx(MSet_pft, plot=TRUE) # takes time
#print(outliers2)

#----------filering cross reactive probes:output of 
#-----------rmSNPandCH is matrix with rows after filteration--------------
# 
b_val_filt_noob <- rmSNPandCH(b_AfterfiltandNoob, dist = 2, 
                         mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY=TRUE)
nrow(b_val_filt_noob) #778181
dim(b_val_filt_noob)  # 778181 x25


#---------------------visualization of density plot of beta of all probes--------------
#-----------------comparing normalize with unnormalized------------------------
par(mfrow=c(1,2))
# Plot distributions prior to normalization for sample 1
plotBetasByType(MSetHvR[,1],main="Raw")
plotBetasByType(MSet_pft[,1],main="pfiltered in RG extended set")
dev.off()

#-------------qc by champ-------------
dev.off()
QC_champ <- champ.QC(b_val_filt_noob,pheno=tarHvR$Sample_Group,
                     mdsPlot=TRUE,
                     densityPlot=TRUE,
                     dendrogram=TRUE,
                     PDFplot=TRUE,
                     Rplot=TRUE,
                     Feature.sel="None")
QC.GUI(b_val_filt_noob, "EPIC") # do not run this one

#------ champ BMIQ normalization -----------

bval_noobBmiq <- champ.norm(beta=b_val_filt_noob,arraytype="EPIC",
                            method="BMIQ", plotBMIQ=TRUE,cores=3)

densityPlot(bval_noobBmiq , sampGroups = tarHvR$Sample_Group,
            main="Noob and BMIQ normalization")

#------- getting the m VAUES FROM BETA MATRIX---------
Mval_noobBmiq<- logit2(bval_noobBmiq) 

#------------visualization of mds plot-------------
par(mfrow=c(2,2))
pal <- brewer.pal(8,"Dark2")
plotMDS(getM(MSetHvR), top=1000, gene.selection="common", pch= 19,
        col=pal[factor(tarHvR$Sample_Group)], dim =c(1,2), cex=0.8, 
        main= "raw")
plotMDS(getM(MSet_pft), top=1000, gene.selection="common", pch= 19,
        col=pal[factor(tarHvR$Sample_Group)], dim =c(1,2), cex=0.8, 
        main= "before filtering noob normalization")
plotMDS(Mval_noobBmiq, top=5000, gene.selection="common", pch= 19,
        col=pal[factor(tarHvR$Sample_Group)], dim =c(1,2), cex=0.8, 
        main= "M val of noob+BMIQ- top 5000")
legend("topleft", legend=levels(factor(tarHvR$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
plotMDS(Mval_noobBmiq, top=10000, gene.selection="common", pch= 19,
        col=pal[factor(tarHvR$Sample_Group)], dim =c(1,2), cex=0.8, 
        main= "M val of noob+BMIQ- top 10000")
legend("topleft", legend=levels(factor(tarHvR$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
dev.off()

mdsPlot(Mval_noobBmiq, sampGroups=as.factor(tarHvR$Slide), pch=19)

#------------ performing SVD to check what is causing batch effect--------
bval_noobBmiq <- as.data.frame(bval_noobBmiq)
svd_champ<-champ.SVD(bval_noobBmiq,pd=tarHvR, RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages/")
dev.off()
#----------estimation of the blood cell components---------
cellCounts <- estimateCellCounts2(RGsetHvR, compositeCellType = "Blood",
                                  referencePlatform= "IlluminaHumanMethylationEPIC",
                                  processMethod = "preprocessNoob",
                                  probeSelect = "IDOL",
                                  meanPlot= TRUE)
head(cellCounts)
tarHvR <- transform(tarHvR, sampleID = paste(Slide, Array, sep = "_"))
cell.data<- as.data.frame(cellCounts[["prop"]])
str(cell.data)
all(rownames(cell.data) == tarHvR$sampleID)
cell.data$sample_group <- tarHvR$Sample_Group
head(cell.data)

melted_celldata<-  melt(cell.data, id.vars = "sample_group", 
                measure.vars = c("CD8T","CD4T","NK","Bcell","Mono","Neu"))
head(melted_celldata)

plot <- ggplot(melted_celldata, aes(x=variable,y= value, 
                        fill =sample_group ))+geom_boxplot()+
  labs(title="Blood cell proportion in different sample groups",
       x="Different blood cells", y = "Cell proportion")+ theme_classic()

#--- calculating the difference between two groups by wilcox test-------
#-------ggpubr package was used for this----------
dev.off()

compare_means(CD8T ~ sample_group, data = cell.data)
b <- ggboxplot(cell.data, x = "sample_group", y = "CD8T",color="sample_group",
               palette = "jco")+stat_compare_means(label = "p.signif")
compare_means(CD4T ~ sample_group, data = cell.data)
c <- ggboxplot(cell.data, x = "sample_group", y = "CD4T",color="sample_group",
          palette = "jco")+stat_compare_means(label = "p.signif")
compare_means(NK ~ sample_group, data = cell.data)
d<- ggboxplot(cell.data, x = "sample_group", y = "NK",color="sample_group",
          palette = "jco")+stat_compare_means(label = "p.signif")
compare_means(Bcell ~ sample_group, data = cell.data)
e<- ggboxplot(cell.data, x = "sample_group", y = "Bcell",color="sample_group",
          palette = "jco")+stat_compare_means()
compare_means(Mono ~ sample_group, data = cell.data)
f<- ggboxplot(cell.data, x = "sample_group", y = "Mono",color="sample_group",
          palette = "jco")+stat_compare_means()
compare_means(Neu ~ sample_group, data = cell.data)
g<- ggboxplot(cell.data, x = "sample_group", y = "Neu",color="sample_group",
          palette = "jco")+stat_compare_means(label = "p.signif")
# par does not work for ggpubr so use ggarrange
ggarrange(b,c,d,e,f,g+ rremove("x.text"), 
          labels = c("A", "B", "C","D","E","F"),
          ncol = 3, nrow = 2)
#--------- batch correcting by combact(champ package)--------
Mval_noobBmiq<- logit2(bval_noobBmiq)
bval_noobBmiq_bc <- champ.runCombat(bval_noobBmiq,
                                    pd=tarHvR,batchname=c("Slide"),
                                    logitTrans=TRUE)
dim(bval_noobBmiq_bc)

# checking again the SVD
bval_noobBmiq_bc <- as.data.frame(bval_noobBmiq_bc)
svd_champ1<-champ.SVD(bval_noobBmiq_bc ,pd=tarHvR, RGEffect=FALSE,
                     PDFplot=TRUE,
                     Rplot=TRUE,
                     resultsDir="./CHAMP_SVDimages/")

tar_cell <- cbind(tarHvR,cell.data)
tar_cell <- tar_cell [, -c(12,19)]
svd_champ2<-champ.SVD(bval_noobBmiq_bc ,pd=tar_cell, RGEffect=FALSE,
                      PDFplot=TRUE,
                      Rplot=TRUE,
                      resultsDir="./CHAMP_SVDimages/asians_cellproportion")
dev.off()

Mval_noobBmiq_bc<- logit2(bval_noobBmiq_bc)


#---------- dmp finder by champ package--------------------------------

champ_dmp <- champ.DMP(beta = bval_noobBmiq_bc,pheno = tarHvR$Sample_Group,
                       compare.group = NULL,
                       adjPVal = 0.045,
                       adjust.method = "BH",
                       arraytype = "EPIC") # 85463 significantdmp withP-value below 0.05


DMP.GUI(DMP=champ_dmp[[1]],beta=bval_noobBmiq_bc,pheno=tarHvR$Sample_Group)
#--------- making DMP as dataframe and data visualization------------------------------------------
dmp_champ_dataframe<- as.data.frame(champ_dmp[[1]])
dim(dmp_champ_dataframe) #85463    20
head(dmp_champ_dataframe,2)
dmp_champ_dataframe[order(dmp_champ_dataframe$adj.P.Val, decreasing = FALSE),]  
summary(dmp_champ_dataframe$adj.P.Val)
ggplot(dmp_champ_dataframe, aes(x = adj.P.Val)) + 
  geom_histogram(color = "black", fill = "#1FA187") +
  labs(title = "histogram id Adjusted pVAL of DMPS in South Asian population",
       x =  NULL) +
  theme_minimal() 

#----------plot to see DMP in different chromosomes-------------
DMP_df2$Methylation_status <- "None"
# manhattan plot:very bad
DMP_df2 <- DMP_df2 %>% arrange(chr)
ggplot(DMP_df2, aes(x = chr, y= logFC, fill= Methylation_status)) + 
  geom_col(position = "identity",color = "black",width = 0.4) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme_minimal()+
  labs(title = "DMPs of South Asian RA in all chromosomes" , y= "Log Fold Change",
       x= "Chromosomes")

#---------making column in the DMP file if the CpG is hypo or hyper methylated 
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
DMP_df2$Methylation_status[DMP_df2$logFC > 0.05 & DMP_df2$adj.P.Val < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
DMP_df2$Methylation_status[DMP_df2$logFC < -0.05 & DMP_df2$adj.P.Val < 0.01] <- "Down"
DMP_df2$labels<- NA
DMP_df2$labels[DMP_df2$Methylation_status != "None"] <- DMP_df2$gene[DMP_df2$Methylation_status != "None"]
write.xlsx(DMP_df2, "dmp_sa.xlsx", col.names = TRUE, row.names = TRUE, append = FALSE)

#--------trying to make a volcano plot------------
mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("Down", "UP", "None")

p <- ggplot(data=DMP_df2, aes(x=logFC, y=-log10(adj.P.Val), label=gene)) + 
  geom_point() +theme_minimal()+
  geom_vline(xintercept=c(-0.05, 0.05), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")+
  geom_text_repel(data=DMP_df2[1:40,],aes(x = logFC, y = -log10(adj.P.Val),label=gene))+
  scale_color_manual(values=c('gray','red',"green"))+
  labs(title = "Volcano Plot with DMPs in South Asian RA patients" ,
       x= "Log Fold Change", y= "Log 10 Adj P value")
p+scale_color_manual(values=c('gray','red',"green"))
  p+geom_text_repel(data=DMP_df2[1:20,],aes(x = logFC, y = -log10(adj.P.Val),label=gene))+
  geom_point(data=DMP_df2[1:20,],aes(x = logFC, y = -log10(adj.P.Val),color="red", size=1.5,alpha=0.4))+
  scale_color_manual(values=c('gray','red',"green"))

DMP_df3 <- subset(DMP_df2, DMP_df2$adj.P.Val <0.01 | DMP_df2$adj.P.Val == 0.01)
ggplot(data=DMP_df3, aes(x=logFC, y=-log10(adj.P.Val), col=feature, label=DMP_df3$gene)) + 
  geom_point() +theme_minimal()+
  scale_colour_manual(values = mycolors)+
  labs(title = "Volcano Plot with DMPs in South Asian RA patients" ,
       x= "Log Fold Change", y= "Log 10 Adj P value")

# ---------making volcano plot with enhanced volcano (best way for this data)---------------
DMP_df4 <- DMP_df2[,c("gene","logFC","adj.P.Val")]
DMP_df4 <-colnames("ID","logFC","adj.P.Val")
EnhancedVolcano(DMP_df4,
                lab = rownames(DMP_df4),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 10e-3,
                title = 'Volcano Plot with DMPs in South Asian RA patients',
                FCcutoff = 0.05,xlim = c(-0.2, 0.2),legendPosition = 'right',
                ylim = c(0,7.5),
                legendLabSize = 10,
                legendIconSize = 2.0,
                )

#--------- visulaization of log fc in all chrmosomes)

dmp_df2_hyper <- subset(DMP_df2, DMP_df2$Methylation_status =="Hyper")

ggplot(dmp_champ_dataframe, aes(x=chr, y=))
summary(dmp_champ_dataframe$logFC)
ggplot(dmp_champ_dataframe, aes(x = logFC)) + 
  geom_histogram(color = "black", fill = "#1FA187") +
  labs(title = "histogram of logFC of DMPS in South Asian population",
       x =  NULL) +
  theme_minimal() 
bval_noobBmiq_bc <- as.matrix(bval_noobBmiq_bc)

# ---------- susetting adjusted  p<0.01 dmp ---------
sigdmp <- filter(dmp_champ_dataframe, dmp_champ_dataframe$adj.P.Val< 10e-3) #10e-3 means 0.010
dim(sigdmp)#25895  cps
#----- visulization spefic DMPs-------------
hla_sa <- dmp_champ_dataframe%>%filter(grepl("HLA", gene))
remove_samples <- c("HHLA2", "SCHLAP1", "HHLA3","HHLA1")
hla_sa <-  hla_sa[!(hla_sa$gene %in% remove_samples), ]
ggplot(hla_sa, aes(x = gene, y = logFC, fill = feature)) +
   geom_col(position = "identity",color = "black",width = 0.4)+coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y = "Log Fold Change (logFC)", x = "HLA Genes")

mirna_sa <- sigdmp%>%filter(grepl("MIR", gene))
mirna_sa_tss200  <-mirna_sa%>%filter(grepl("TSS200", feature))
mirna_sa_tss200  <-mirna_sa_tss200 %>% arrange(mirna_sa_tss200$logFC)
ggplot(mirna_sa_tss200, aes(x = gene, y = logFC, fill = cgi)) +
  geom_col(position = "identity",color = "black",width = 0.4)+coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
labs(y = "Log Fold Change (logFC)", x = "MicroRNA (MiRNA) Genes")+ theme_bw()

ggplot(mirna_sa, aes(x = gene, y = logFC, fill = cgi)) +
  geom_col(position = "identity",color = "black",width = 0.4)+coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y = "Log Fold Change (logFC)", x = "MicroRNA (MiRNA) Genes")+ theme_bw()



#-------- checking for DMR-----------------

champ_dmr <- champ.DMR(beta=bval_noobBmiq_bc,pheno=tarHvR$status,
                       arraytype="EPIC",method="Bumphunter",maxGap = 500)
head(champ_dmr)
DMR.GUI(DMR=champ_dmr, beta=bval_noobBmiq_bc, pheno=tarHvR$status,
        arraytype="EPIC")

champ_dmr_df <- as.data.frame(champ_dmr[["BumphunterDMR"]])
write.xlsx(champ_dmr_df, "dmR_sa.xlsx", colNames = TRUE, rowNames = TRUE, append = FALSE)


#------------trying to do dm carte so that we can make tracks-----------
statusH <- factor(tarHvR$status)
# use the above to create a design matrix
design1 <- model.matrix(~1+statusH, data=tarHvR)
#colnames(design1) <- levels(statusH)[-1]
Mval_noobBmiq_bc<- logit2(bval_noobBmiq_bc)
Mval_noobBmiq_bc <- as.matrix(Mval_noobBmiq_bc)
myAnnotation <- cpg.annotate(object = Mval_noobBmiq_bc, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design1, 
                             coef= "statusHRA",arraytype = "EPIC",fdr = 0.01)
str(myAnnotation)
DMRs_SA <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs_SA)
results.ranges
pal <- brewer.pal(8,"Dark2")
groups <- pal[1:length(unique(tarHvR$status))]
names(groups) <- levels(factor(tarHvR$status))
cols <- groups[as.character(factor(tarHvR$status))]
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 2, CpGs = bval_noobBmiq_bc, phen.col = cols, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")
write.csv(results.ranges,"DMR_sa_crate", sep=",")

# trying to visulaize the results

#------ finding fuctional regions in the data using mCSEA (PORMOTER REGION)---------------------

myRank <- dmp_champ_dataframe$logFC
names(myRank) <- rownames(dmp_champ_dataframe)

# Reshape the phenotype data to a format suitable for mCSEA
pData(gSetHvR2)
pheno <- as.data.frame(pData(gSetHvR2))
pheno <- pheno[,"Sample_Group", drop=FALSE]

# Run the mCSEA
mcsea_res  <- mCSEATest(myRank,
                        bval_noobBmiq_bc,
                       pheno,
                       regionsTypes = "promoters",
                       platform = "EPIC")
head(mcsea_res$promoters)

mCSEAresData <- as.data.frame(mcsea_res[["promoters"]])
sigDSEA_data <- mCSEAresData[mCSEAresData$padj<0.01 ,]
mCSEAPlot(mcsea_res,
          regionType = "promoters",
          dmrName = "RNF19A",leadingEdge = TRUE,CGI = TRUE,
          transcriptAnnotation = "symbol",
          makePDF = FALSE)
write.csv(sigDSEA_data,"sig_CSEAres_southAsians.csv",sep=",")

#--------------GSE ANALYSIS-------------------------------
#---------------------
myGSEA <- champ.GSEA(beta=bval_noobBmiq_bc,DMP=champ_dmp[[1]], 
                     DMR=champ_dmr, pheno=tarHvR$status,
                     arraytype="EPIC",adjPval=0.05, method="gometh", Rplot=TRUE)

myGSEA_dmr <- champ.GSEA(beta=bval_noobBmiq_bc,DMP=champ_dmp[[1]], 
                     DMR=champ_dmr, pheno=tarHvR$status,
                     arraytype="EPIC",adjPval=0.05, method="gometh", Rplot=TRUE)
head(myGSEA$DMP)
tail(myGSEA$DMP)

gsea_dmp <- as.data.frame(myGSEA[["DMP"]])
gsea_dmp <- subset(gsea_dmp, gsea_dmp$FDR <0.01)
str(gsea_dmp)
table(gsea_dmp$ONTOLOGY)

gsea_dmp_MF <- subset(gsea_dmp, gsea_dmp$ONTOLOGY == "MF")

head(gsea_dmp_MF)
summary(gsea_dmp_MF )

gsea_dmp_MF_fdr<- subset(gsea_dmp_MF, gsea_dmp_MF$FDR < 0.01 | gsea_dmp_MF$FDR == 0.01)
gsea_dmp_MF_fdr1 <- gsea_dmp_MF_fdr %>% arrange(desc(DE))
summary(gsea_dmp_MF_fdr1)
MF_plot <- ggplot(gsea_dmp_MF_fdr,aes(x= TERM, y=FDR))+
  geom_bar(stat='identity',fill="steelblue")+theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))

gsea_dmp_BP <- subset(gsea_dmp, gsea_dmp$ONTOLOGY == "BP")
gsea_dmp_BP_fdr <- subset(gsea_dmp_BP, gsea_dmp_BP$FDR < 0.01)
gsea_dmp_BP_fdr1 <- gsea_dmp_BP_fdr %>% arrange(FDR)

gseaSA_fdr1<- gsea_dmp_BP_fdr1[1:15,]
ggplot(gseaSA_fdr1,aes(x= TERM, y=-log10(FDR), size= DE))+
  geom_point(stat='identity',fill="blue", color="darkred")+theme_bw()+coord_flip()
theme(axis.text.x = element_text(angle = 40,hjust=2))


write.xlsx(gsea_dmp_BP_fdr1, "gseaDMPsaBP.xlsx", colNames = TRUE,rowNames=TRUE)







