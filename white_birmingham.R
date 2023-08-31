library(minfi)
library(tibble)
library(dplyr)
library("Gviz")
library("RColorBrewer")
library(minfi)
library(ChAMP)
#BiocManager::install("DMRcate")
library("DMRcate")
#BiocManager::install("DMRcatedata")
library("DMRcatedata")
#BiocManager::install("stringr")
library("stringr")
library("mCSEA")
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", force= TRUE)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#BiocManager::install("mCSEA")
library("mCSEA")
library("limma")
library(ggplot2)

#-------------- data directory------------
# pointing R  towards the directory in which data is saved.
datadir <- "Desktop/meth_data/"
# seeing the list of the list of all files: recursive open all the subfiles
list.files(datadir, recursive = TRUE)

#------------read in the sample sheet for the experiment----------
Bir_white <- read.metharray.sheet(datadir, pattern="SampleInfo_all_White.csv")
str(Bir_white) 
#--------subsetting netherland polulation----------
Bir_white <- subset(Bir_white, Bir_white$Subgroup=="Birminghan")
#------ sanity check alert and cleaning a bit----------
head(Bir_white)
Bir_tar <- Bir_white[,-c(2,3,10,11)]

#------- making RGchannel object--------------
rgchan_bir <- read.metharray.exp(targets=Bir_tar) # takes time be patient
# probes 1051943 37

# -----------making RG channel extended object-----------------
rgchan_birEx <- read.metharray.exp(targets=Bir_tar, extended = TRUE)

## -----------1st filter  using watermelon package-----------------------
Fil_Msetbir <- pfilter(rgchan_birEx, pnthresh = 0.01)
#0 samples having 1 % of sites with a detection p-value greater than 0.01 were removed 
#Samples removed:  
# 14759 sites were removed as beadcount <3 in 5 % of samples 
#3279 sites having 1 % of samples with a detection p-value greater than 0.01 were removed 
#---------- checking for pdata and manifest files-------------------
# checking the metadata for Rgset1
pData(rgchan_bir)
# checking information regarging the probes used in EPIC array
getManifest(rgchan_bir)
# getting manifest
manifest_white <- getManifest(rgchan_bir)
head(getProbeInfo(rgchan_bir))
annotation(rgchan_bir) # "ilm10b4.hg19" 

#------------- methyl set from RG channel --------------
MSet_bir <- preprocessRaw(rgchan_bir)

# ----------getting ratio convert---------------
ratioSetbir <- ratioConvert(MSet_bir  , what = "both", keepCN = TRUE)
# ---------getting genomic ratio set by adding genomic coordinates to ratioset object----- 
gSetbir <- mapToGenome(ratioSetbir)
gSetbir

#---------------sanity and quality control for getsex using genomic ratio set
#----------------and estimate sex from watermelon pkg--------------------------------
bval_noNormBir <- getBeta(gSetbir)
sex_pred <- estimateSex(bval_noNormBir, do_plot=TRUE)

#------------getting beta and density plots from methy set 
phenoDatabir <- pData(MSet_bir)
densityPlot(MSet_bir, sampGroups = Bir_tar$status)

# -----------quality control before normalization--------
qc <- getQC(MSet_bir)
plotQC(qc)
controlStripPlot(rgchan_bir, controls="BISULFITE CONVERSION II")
#------------- mds plot befor normalization-------------
mdsPlot(rgchan_bir, sampGroups =Bir_tar$status)
mdsPlot(rgchan_bir, sampGroups =Bir_tar$Slide)


#----methyset from RG channel set extended object with p filter-----------
MSet_pbir <- preprocessNoob(rgchan_birEx)
MSet_pbir
#---------------getting ratio convert for the methy set------------
ratioSetpB <- ratioConvert(MSet_pbir , what = "both", keepCN = TRUE)
# getting genomic ratio set by adding genomic coordinates to ratioset object 
gSetpBir <- mapToGenome(ratioSetpB)
gSetpBir

#---------getting beta,copy number and M values
b_fN <- getBeta(gSetpBir)
nrow(b_fN) # with pvalue <0.01 865859 cpg probes r there

#---------checking for outlier again using watermelon package
#----- Needed only for sanity check---------------
outliers2 <- outlyx(MSet_pbir, plot=TRUE) # takes time
print(outliers2)

#----------filering cross reactive probes:output of 
#-----------rmSNPandCH is matrix with rows after filteration--------------
# 
b_BM_bir <- rmSNPandCH(b_fN , dist = 2, 
                              mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY=TRUE)
nrow(b_BM_bir) #787087
dim(b_BM_bir)  # 787087x37


#---------------------visualization of density plot of beta of all probes--------------
#-----------------comparing normalize with unnormalized------------------------
par(mfrow=c(1,2))
# Plot distributions prior to normalization for sample 1
plotBetasByType(MSet_bir[,1],main="Raw")
plotBetasByType(MSet_pbir[,1],main="pfiltered in RG extended set")
dev.off()

#-------------qc by champ-------------
dev.off()
QC_champ_bir <- champ.QC(b_BM_bir,pheno=Bir_tar$status,
                     mdsPlot=TRUE,
                     densityPlot=TRUE,
                     dendrogram=TRUE,
                     PDFplot=TRUE,
                     Rplot=TRUE,
                     Feature.sel="None")
QC.GUI(b_BM_bir, "EPIC") # do not run this one

#------ champ BMIQ normalization -----------

bval_NB_bir <- champ.norm(beta=b_BM_bir,arraytype="EPIC",
                            method="BMIQ", plotBMIQ=TRUE,cores=3)

densityPlot(bval_NB_bir , sampGroups = Bir_tar$status,
            main="Noob and BMIQ normalization")

#------- getting the m VAUES FROM BETA MATRIX---------
Mval_NB_bir<- logit2(bval_NB_bir) 

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
bval_NB_bir <- as.data.frame(bval_NB_bir)
svd_champ_bir<-champ.SVD(bval_NB_bir,pd=Bir_tar, RGEffect=FALSE,
                     PDFplot=TRUE,
                     Rplot=TRUE,
                     resultsDir="./CHAMP_SVDimages/")
dev.off()
#----------estimation of the blood cell components---------
cellCounts_bir <- estimateCellCounts2(rgchan_bir, compositeCellType = "Blood",
                                  referencePlatform= "IlluminaHumanMethylationEPIC",
                                  processMethod = "preprocessNoob",
                                  probeSelect = "IDOL",
                                  meanPlot= TRUE)
head(cellCounts_bir)
Bir_tar <- transform(Bir_tar, sampleID = paste(Slide, Array, sep = "_"))
cell.data_bir<- as.data.frame(cellCounts_bir[["prop"]])
str(cell.data_bir)
all(rownames(cell.data_bir) == Bir_tar$sampleID)
cell.data_bir$sample_group <- Bir_tar$status
head(cell.data_bir)

melted_celldata_b<-  melt(cell.data_bir, id.vars = "sample_group", 
                        measure.vars = c("CD8T","CD4T","NK","Bcell","Mono","Neu"))
head(melted_celldata_b)

plot_b <- ggplot(melted_celldata_b, aes(x=variable,y= value, 
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

#--------- batch correcting by combact(champ package) --------
Mval_NB_bir<- logit2(bval_NB_bir)

# checking again the SVD 
bval_noobBmiq_bc <- as.data.frame(bval_noobBmiq_bc)
svd_champ1<-champ.SVD(bval_NB_bir ,pd=Bir_tar, RGEffect=FALSE,
                      PDFplot=TRUE,
                      Rplot=TRUE,
                      resultsDir="./CHAMP_SVDimages/")
# ------ running combat to remove batch effect-----------
bval_NB_bcB <- champ.runCombat(bval_NB_bir,variablename = "sex",pd=Bir_tar,
                               batchname=c("Slide"),
                               logitTrans=TRUE)
dim(bval_noobBmiq_bc)
#------ failed to remove batch effect due to too many covariate-------

#---------- dmp finder by champ package--------------------------------

champ_dmp_Bir <- champ.DMP(beta = bval_NB_bir,pheno = Bir_tar$status,
                       compare.group = NULL,
                       adjPVal = 0.045,
                       adjust.method = "BH",
                       arraytype = "EPIC") # 87235 significantdmp withP-value below 0.05
DMP.GUI(DMP=champ_dmp[[1]],beta=bval_noobBmiq_bc,pheno=tarHvR$Sample_Group)
#-----------------sanity checks and data visualization --------
head(champ_dmp_Bir[[1]])

dmp_chp_dfBir<- as.data.frame(champ_dmp_Bir[[1]])
dim(dmp_chp_dfBir) #87235    20
head(dmp_chp_dfBir,2)
dmp_chp_dfBir[order(dmp_chp_dfBir$adj.P.Val, decreasing = FALSE),]  
summary(dmp_chp_dfBir$adj.P.Val)
ggplot(dmp_champ_dataframe, aes(x = adj.P.Val)) + 
  geom_histogram(color = "black", fill = "#1FA187") +
  labs(title = "histogram id Adjusted pVAL of DMPS in South Asian population",
       x =  NULL) +
  theme_minimal() 
summary(dmp_champ_dataframe$logFC)
ggplot(dmp_champ_dataframe, aes(x = logFC)) + 
  geom_histogram(color = "black", fill = "#1FA187") +
  labs(title = "histogram of logFC of DMPS in South Asian population",
       x =  NULL) +
  theme_minimal() 
bval_noobBmiq_bc <- as.matrix(bval_noobBmiq_bc)
#--------- subsetting adj p val <0.01 dmp -----------------

sigdmp_Bir <- filter(dmp_chp_dfBir, dmp_chp_dfBir$adj.P.Val< 10e-3) #10e-3 means 0.010
dim(sigdmp_Bir)#46227  cps
#----- visulization spefic DMPs-------------
hla_sa_Bir <- dmp_chp_dfBir%>%filter(grepl("HLA", gene))

hla_sa_Bir1<-hla_sa_Bir[!(hla_sa_Bir$gene=="HHLA2"),]
hla_sa_Bir1<-hla_sa_Bir1[!(hla_sa_Bir1$gene=="SCHLAP1"),] 
hla_sa_Bir1<-hla_sa_Bir1[!(hla_sa_Bir1$gene=="HHLA3"),] 


ggplot(hla_sa_Bir1, aes(x = gene, y = logFC, fill = feature)) +
  geom_col(position = "identity",color = "black",width = 0.4)+coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x ="HLA genes",y="Log Fold Change")

mirna_Bir <- sigdmp_Bir%>%filter(grepl("MIR", gene))
dim(mirna_Bir)

mirna_Bir$Methylation_status <- "NA"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
mirna_Bir$Methylation_status[mirna_Bir$logFC > 0.02] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
mirna_Bir$Methylation_status[mirna_Bir$logFC < -0.02] <- "Down"

mirna_down <- subset(mirna_Bir, mirna_Bir$Methylation_status =="Down")
mirna_up <- subset(mirna_Bir, mirna_Bir$Methylation_status =="UP")

mirna_up <- mirna_up %>% arrange(desc(logFC))
mirna_down <- mirna_down %>% arrange(logFC)
mirna_up1 <- mirna_up[1:10,]
mirna_down1 <- mirna_down [1:10,]
mirna_fltw<- rbind(mirna_up1,mirna_down1 )
ggplot(mirna_fltw, aes(x = gene, y = logFC, fill = feature)) +
  geom_col(position = "identity",color = "black",width = 0.4)+coord_flip()+
  theme(axis.text.y = element_text(size =2, angle = 90, vjust = 0.4, hjust=1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust=1))+theme_minimal()+
  labs(x ="MicroRNA  genes",y="Log Fold Change")





#-----finding DMR in champ-----------
bval_NB_bir <- as.data.frame(bval_NB_bir)
champ_dmr_bir <- champ.DMR(beta=bval_NB_bir,pheno=Bir_tar$status,
                       arraytype="EPIC",method="Bumphunter")

DMR.GUI(DMR=champ_dmr_bir , beta=bval_noobBmiq, pheno=tarHvR$status,
        arraytype="EPIC")


champDMR_Bir <- as.data.frame(champ_dmr_bir[["BumphunterDMR"]])
write.xlsx(champDMR_Bir, "DMR_WHITE.xlsx", 
           colNames = TRUE, rowNames = TRUE, append = FALSE)

########-------- GSEA------------

myGSEA_WHITE <- champ.GSEA(beta=bval_NB_bir,DMP=champ_dmp_Bir[[1]], 
                     DMR=champ_dmr_bir, pheno=Bir_tar$status,
                     arraytype="EPIC",adjPval=0.05, method="gometh", Rplot=TRUE)

head(myGSEA_WHITE$DMP)

gsea_dmp_white <- as.data.frame(myGSEA_WHITE[["DMP"]])
str(gsea_dmp_white)

table(gsea_dmp_white $ONTOLOGY)
gseawhitMF<- subset(gsea_dmp_white, gsea_dmp_white$ONTOLOGY == "MF")

head(gsea_dmp_MF)
summary(gsea_dmp_MF )

gseadmpMFW_fdr<- subset(gseawhitMF, gseawhitMF$FDR < 0.01)
gseadmpMFW_fdr<- gseadmpMFW_fdr %>% arrange(DE)
gseadmpMFW_fdr1<- gseadmpMFW_fdr[1:15,]
 ggplot(gseadmpMFW_fdr1,aes(x= TERM, y=-log10(FDR)))+
  geom_bar(stat='identity',fill="blue")+theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))

gseaBPwhite<- subset(gsea_dmp_white, gsea_dmp_white$ONTOLOGY == "BP")
gseaBPW_fdr <- subset(gseaBPwhite, gseaBPwhite$FDR < 0.01 )
gseaBPW_fdr <- gseaBPW_fdr %>% arrange(DE)


gseaBPW_fdr1<- gseaBPW_fdr[1:15,]
ggplot(gseaBPW_fdr1,aes(x= TERM, y=-log10(FDR), size= DE))+
  geom_point(stat='identity',fill="blue", color="darkred")+theme_bw()+coord_flip()
  theme(axis.text.x = element_text(angle = 40,hjust=2))
  
  ggplot(gseaBPW_fdr1,aes(x= TERM, y=-log10(FDR), size= DE))+
  geom_point(alpha=0.5) +coord_fip()
    scale_size(range = c(.1, 24)) +
    theme(legend.position="bottom") 
    
  

#--------------csea promoter------------
myRank1 <- dmp_chp_dfBir$logFC
names(myRank1) <- rownames(dmp_chp_dfBir)

# Reshape the phenotype data to a format suitable for mCSEA
pData(gSetpBir)
phenoB <- as.data.frame(pData(gSetpBir))
phenoB <- phenoB[,"status", drop=FALSE]

# Run the mCSEA
mcsea_resB  <- mCSEATest(myRank1,
                        bval_NB_bir,
                        phenoB,
                        regionsTypes = "promoters",
                        platform = "EPIC")
head(mcsea_resB$promoters)

mCSEAres_Bir <- as.data.frame(mcsea_resB [["promoters"]])
sigDSEA_bir <- mCSEAres_Bir[mCSEAres_Bir$padj<0.01 ,]
mCSEAPlot(mcsea_resB,
          regionType = "promoters",
          dmrName = "BRD2",leadingEdge = TRUE,CGI = TRUE,
          transcriptAnnotation = "symbol",
          makePDF = FALSE)

Sa_hypermethylatedRA <- subset(sigdmp, sigdmp$logFC > 0 | sigdmp$logFC ==  "0")
Sa_hypormethylatedRA <- subset(sigdmp, sigdmp$logFC < 0 )

saHypo<- as.character(Sa_hypormethylatedRA$gene)
saHyper<- as.character(Sa_hypermethylatedRA$gene)

White_hyperRA <- subset(sigdmp_Bir, sigdmp_Bir$logFC > 0 | sigdmp_Bir$logFC ==  "0")
white_hypoRA <- subset(sigdmp_Bir, sigdmp_Bir$logFC < 0 )

whiteHypo<- as.character(white_hypoRA$gene)
whiteHyper<- as.character(White_hyperRA$gene)


x <- list(SA_RA_Hypo= saHypo,
          SA_RA_hyper = saHyper,
          White_hyper = whiteHyper,
           White_hypo = whiteHypo)
myCol <- brewer.pal(3, "Pastel2")
ggVennDiagram(x, color = "black", lwd = 0.8, lty = 1)+
scale_fill_gradient(low="lightblue",high = "pink")

venn.diagram(x,
  category.names = c("South Asian RA: hypomethylated" ,"South Asian RA: hypermethylated",
                     "White RA Hypermethylated", "White RA Hypomethylated" ),
  filename = 'GENES_venn_diagramm.png',
  output=TRUE,lwd = 2,
  lty = 'blank',
  fill = myCol)





ggVennDiagram(x, label_alpha = 0)
 








