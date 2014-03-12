/data/customers/Musculoskeletal/Yaobo/Osteoarthritic_Hip_Cartilage/R_Analysis

## cummeRbund analysis of zebrafish data
## ************* Based on STAR ALIGNMENT ***************

biocLite("cummeRbund")
library("cummeRbund")

# working directory should be the directory containing cuffdiff output 
setwd("/data/customers/Musculoskeletal/Yaobo/Osteoarthritic_Hip_Cartilage/R_Analysis")
setwd("/data/customers/Corbett/STAR_cuffdiff/cuffdiff_out_reduced/")

cuff<-readCufflinks(gtfFile='merged.gtf',genome='hg19')
cuff = readCufflinks()

#dend.rep.isoforms <- csDendro(isoforms(cuff),replicates=T)
dend.rep.genes <- csDendro(genes(cuff),replicates=T)

gg1 <- ggdendrogram(dend.rep.genes, rotate = TRUE, size = 4, theme_dendro = FALSE) +
  theme(panel.background = element_blank()) +
  labs(title = "Dendogram of NOF and OA RNA Seq Samples (genes)")

#gg <- ggdendrogram(dend.rep.isoforms, rotate = TRUE, size = 4, theme_dendro = FALSE) +
#  theme(panel.background = element_blank()) +
#  labs(title = "Dendogram of NOF and OA RNA Seq Samples (isoforms)")
#print(gg)

cb_pallette <- c()
cb_pallette[1:16] <- 'black'
genes.PCA = MDSplot(genes(cuff),replicates=T) +
  theme(panel.background = element_blank()) +
  theme(legend.position="none") + 
  scale_color_manual(values = cb_pallette) +
  labs(title = "")
  #scale_x_continuous(limits = c(-.5, .5)) +
  #scale_y_continuous(limits = c(-.5, .5))
genes.PCA

# Boxplot
cb_pallette <- c()
cb_pallette[1:16] <- 'grey'
bplot = csBoxplot(genes(cuff), replicates = TRUE) +
  theme(panel.background = element_blank()) +
  theme(legend.position="none") + 
  scale_fill_manual(values = cb_pallette) +
  labs(title = "")


v<-csVolcano(genes(cuff),showSignificant=TRUE, "OA","NOF", alpha=0.05) +
  geom_hline(size=.5,aes(yintercept=2.5)) +
  theme(legend.position="none") + 
  theme(panel.background = element_blank()) +
  labs(title = "") 
v

pdf("RNA_Seq_Images.pdf")
print(gg1)
print(genes.PCA)
print(bplot)
print(v)
dev.off()

gene_fpkm_data = read.csv("heatmap.csv")
##'HEATMAP

diff_exp <- as.matrix(gene_fpkm_data)
gene_fpkm_data <- diff_exp[,3:4]
rownames(gene_fpkm_data) <- diff_exp[,1]
gene_fpkm_data <- gene_fpkm_data[,3:4]
class(gene_fpkm_data) <- "numeric"
my_palette <- colorRampPalette(c("white", "grey", "black"))(n = 1000)
par(cex.main=0.8)
pdf("Diff_Exp_Genes_FPKM_Heatmap_BW.pdf")
heatmap.2(t(gene_fpkm_data_log), col=my_palette, density.info="none",
          trace="none", cexRow=0.2, cexCol=0.8,
          main="NOF Vs OA Differentially Expressed Genes (-log10(FPKM))")
dev.off()

heatmap.2(gene_fpkm_data, col="redblue", density.info="none",
          trace="none", cexRow=0.2, cexCol=0.8, 
          main=paste(unlist(strsplit(gene_lists_1[5], "_"))[3:4], collapse=" Vs "))


##'GENE SETS
##'DIFFERENTIALLY EXPRESSED GENES

myGeneIDs <- read.table("gene_names.csv", header=FALSE, stringsAsFactors=FALSE)
myGenes   <- getGenes(cuff,geneIDs)
head(fpkm(myGenes))

geneIDs <- c()
for(i in 1:length(myGeneIDs))
{
  geneIDs[i] <- myGeneIDs[i]  
}

h <- csHeatmap(myGenes,cluster="both",replicates=T)
h
b <- expressionBarplot(myGenes)
b  
s <- csScatter(myGenes,"NOF","OA"),smooth=T)
s
v <- csVolcano(myGenes,"NOF","OA")
v

myGeneID <- "XLOC_015821"
cummeRbund:getGene
myGene   <- getGene(cuff,myGeneID)

myGeneId<-"ACAN"
cummeRbund:getGene
myGene<-getGene(cuff,myGeneId)
myGene


gl<-expressionPlot(myGene)
gl
gl.rep<-expressionPlot(myGene,replicates=TRUE)
gl.rep

gl.cds.rep<-expressionPlot(CDS(myGene),replicates=T)
gl.cds.rep

gb<-expressionBarplot(myGene)
gb


source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("Gviz")
biocLite("GenomeGraphs")
library('Gviz')
library('GenomeGraphs')
######################

trackList <- list()
myStart <- min(features(myGene)$start)
myEnd <- max(features(myGene)$end)
myChr <- unique(features(myGene)$seqnames)
genome <- 'hg19'
ideoTrack <- IdeogramTrack(genome = genome, chromosome = myChr)
trackList <- c(trackList,ideoTrack)
axtrack <- GenomeAxisTrack()
trackList <- c(trackList,axtrack)
genetrack <- makeGeneRegionTrack(myGene)
genetrack

trackList <- c(trackList,genetrack)
biomTrack <- BiomartGeneRegionTrack(genome=genome,chromosome=as.character(myChr),
                                      start=myStart,end=myEnd,name="ENSEMBL",showId=T)
trackList<-c(trackList,biomTrack)
conservation <- UcscTrack(genome = genome, chromosome = myChr,
                  track = "Conservation", table = "phyloP100wayAll",
                  from = myStart-2000, to = myEnd+2000, trackType = "DataTrack",
                  start = "start", end = "end", data = "score",
                  type = "hist", window = "auto", col.histogram = "darkblue",
                  fill.histogram = "darkblue", ylim = c(-3.7, 4),
                  name = "Conservation")
trackList <- c(trackList,conservation)
plotTracks(trackList,from=myStart-8000,to=myEnd+2000)
#####################

mySigMat<-sigMatrix(cuff,level='genes',alpha=0.05)
mySigGeneIds<-getSig(cuff,alpha=0.05,level='genes')
NOF_vs_OA.sigIsoformIds<-getSig(cuff,x='NOF' ,y='OA' ,alpha=0.05,level='isoforms')
mySigGenes<-getGenes(cuff,mySigGeneIds)
mySigTable<-getSigTable(cuff,alpha=0.01,level='genes')
myDistHeat<-csDistHeat(genes(cuff),replicates=T)
myDistHeat
ic<-csCluster(myGene,k=4)
icp<-csClusterPlot(ic)
icp

mySimilar<-findSimilar(cuff,"ACAN",n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)


#####################
pdf('ACAN_Plots.pdf')
gb.rep<-expressionPlot(myGene,replicates=T)
print(gb.rep)
gl.iso.rep<-expressionPlot(isoforms(myGene),replicates=T) +
  labs(title = element_text(size = 2))
print(gl.iso.rep)
##'CHANGE WIDTH
genetrack<-makeGeneRegionTrack(myGene)
print(plotTracks(genetrack))
print(plotTracks(trackList,from=myStart-2000,to=myEnd+2000))
print(plotTracks(trackList,from=myStart-13000,to=myEnd+2000))
dev.off()
###################

fpkm_reps_1 <- c()
tmp_ <- as.vector(fpkm_reps[,3])
tmp_ <- as.numeric(tmp_)
for(i in seq(from=1, to=length(tmp_), by=16))
{
  j <- i+15
  fpkm_reps_1 <- rbind(fpkm_reps_1, tmp_[i:j])
}
colnames(fpkm_reps_1) <- as.vector(fpkm_reps[1:16,2])
rownames(fpkm_reps_1) <- as.vector(unique(fpkm_reps[,1]))
fpkm_reps_1<--log10(fpkm_reps_1)
fpkm_reps_1[fpkm_reps_1 == 'Inf'] <- 0
heatmap.2(t(fpkm_reps_1))

fpkm_reps <- cbind(fpkm_reps, myGenes@repFpkm$gene_id)
fpkm_reps <- cbind(fpkm_reps, myGenes@repFpkm$rep_name)
fpkm_reps <- cbind(fpkm_reps, myGenes@repFpkm$fpkm)

my_palette <- colorRampPalette(c("white", "grey", "black"))(n = 1000)
par(cex=1.5)
heatmap.2(t(fpkm_reps_1), col=my_palette, density.info="none",
          trace="none", cexRow=0.8, cexCol=.2,
          main="NOF Vs OA Differentially Expressed Genes (-log10(FPKM))")



