####PACKAGES
install.packages("gplots")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("lumi")
biocLite("arrayQualityMetrics")
biocLite("lumiHumanAll.db")
biocLite("annotate")
biocLite("limma")
biocLite("GOstats")
biocLite("RamiGO")
biocLite("pathview")

library(lumi)
library(gplots)
library(arrayQualityMetrics)
library(lumiHumanAll.db)
library(annotate)
library(limma) 
library(ggplot2)
library("GOstats")
library("RamiGO")
library("GOstats")
library("RamiGO")
library('pathview')
####END PACKAGES

###RAW DATASET & NORMALISATION
#lumib for multiple file reads
raw_data            <- lumiR("Raw_data.txt")
vst_data            <- lumiT(raw_data, method="vst")
rsn_data            <- lumiN(vst_data, method = "rsn")
analysis_ready_data <- lumiQ(rsn_data)
##quantile analysis sucks, don't use it 
##quan_data           <- lumiN(vst_data, method = "quantile")
###END RAW DATASET & NORMALISATION

###RAW DATA VISUALISATION
pdf("Raw_Dataset.pdf")
  plot(raw_data,  main="Density Plot of Raw Dataset")
  boxplot(raw_data,  main="Boxplot of Raw Dataset")
  plot(vst_data,  main="Density Plot of Raw Dataset - Post VST")
  boxplot(vst_data,  main="Boxplot of Raw Dataset - Post VST")
  plot(rsn_data,  main="Density Plot of Raw Dataset - Post RSN")
  boxplot(rsn_data,  main="Boxplot of Raw Dataset - Post RSN")
  #plot(quan_data, main="Density Plot of Raw Dataset - Post Quantile Analysis - No RSN")
  #boxplot(quan_data,  main="Boxplot of Raw Dataset - Post Quantile Analysis - No RSN")
  plot(raw_data[, c(1, 2, 3, 4)], what="MAplot", main=paste(sampleNames(analysis_ready_data)[1:4], collapse=", "), cex=0.5)
  plot(analysis_ready_data[, c(31:36)], cex=0.9, what="MAplot", main=paste(sampleNames(analysis_ready_data)[31:36], collapse=", "))
  plot(analysis_ready_data, what="sampleRelation", method="mds")
dev.off()
###END RAW DATA VISUALISATION

###REMOVAL OF TROUBLESOME ARRAYS
removal_elements <- c(1,2,4,19,24,31:36)
removal          <- sampleNames(analysis_ready_data)[removal_elements]
raw_data_post    <- raw_data[, !(sampleNames(raw_data) %in% removal)]
sampleNames(raw_data_post)
###END REMOVAL OF TROUBLESOME ARRAYS

###OPTIMAL DATASET
vst_data_post            <- lumiT(raw_data_post, method="vst")
rsn_data_post            <- lumiN(vst_data_post, method="rsn")
analysis_ready_data_post <- lumiQ(rsn_data_post)

pdf("Optimal_Dataset.pdf")
  plot(raw_data_post,  main="Density Plot of Optimal Dataset")
  boxplot(raw_data_post,  main="Boxplot of Optimal Dataset")
  plot(vst_data_post,  main="Density Plot of Optimal Dataset - Post VST")
  boxplot(vst_data_post,  main="Boxplot of Optimal Dataset - Post VST")
  plot(rsn_data_post,  main="Density Plot of Optimal Dataset - Post RSN")
  boxplot(rsn_data_post,  main="Boxplot of Optimal Dataset - Post RSN")
  plot(analysis_ready_data_post, main="Density Plot of Optimal Dataset - Post QC")
  boxplot(analysis_ready_data_post,  main="Boxplot of Optimal Dataset - Post QC")
  plot(raw_data_post[, c(1, 2, 3, 4)], what="MAplot", main=paste(sampleNames(analysis_ready_data_post)[1:4], collapse=", "))
  plot(analysis_ready_data_post[, c(26:30)], cex=0.9, what="MAplot", main=paste(sampleNames(analysis_ready_data_post)[26:30], collapse=", "))
  plot(analysis_ready_data_post, what="sampleRelation", method="mds")
dev.off()
###END OF OPTIMAL DATASET

###ARRAY QUALITY METRICS
arrayQualityMetrics(expressionset=analysis_ready_data, outdir="QC_INITIAL_DATASET")
arrayQualityMetrics(expressionset=analysis_ready_data_post, outdir="QC_OPTIMAL_DATASET")
###END OF ARRAY QUALITY METRICS

###DETECTION CALL THRESHOLD REMOVAL
exprs_data                   <- exprs(analysis_ready_data_post)
present_count                <- detectionCall(analysis_ready_data_post) 
filtered_analysis_ready_data <- exprs_data[present_count > 0, ]
###PERCENTAGE OF PROBES REMOVED []
(nrow(filtered_analysis_ready_data)/nrow(analysis_ready_data_post))*100
###END OF DETECTION CALL THRESHOLD REMOVAL

###ARRAY ANNOTATION
probe_list <- rownames(filtered_analysis_ready_data)
nuIDs      <- probeID2nuID(probe_list)[, "nuID"]
symbol     <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name       <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df    <- data.frame(ID = nuIDs, probe_list, symbol, name)
###END OF ARRAY ANNOTATION

###DIFFERENTIAL EXPRESSION ANALYSIS 
treatments  <- c("NOF", "OA")
array_names <- c("NOF","NOF","NOF","OA","OA","OA","NOF","NOF","OA","OA","NOF",
                 "OA","OA","NOF","NOF","OA","NOF","NOF","NOF","OA","OA","OA",
                 "OA","NOF","NOF")
design <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design) <- treatments
num_parameters <- ncol(design)
fit <- lmFit(analysis_ready_data_post, design)
cont_mat <- makeContrasts(NOF-OA, levels=treatments)
fit2 <- contrasts.fit(fit, contrasts=cont_mat)
fit2 <- eBayes(fit2)
#fit2$genes = anno_df

## Filter by fold change (2x) and p-value (0.01) cutoffs
## Adjusted using Benjamini-Hochberg False Discovery Rate
topTable(fit2, coef="NOF - OA", p.value=0.01, lfc=log2(2))
gene_list   = topTable(fit2, coef = "NOF - OA", number = nrow(anno_df))

pdf("Differential_Expression_Plots.pdf")
###GGPLOT2
gene_list$threshold = as.factor(abs(gene_list$logFC) > 1 & gene_list$P.Value < 0.01/nrow(gene_list))#/nrow(anno_df))
g = ggplot(data=gene_list, aes(x=gene_list$logFC, y=-log10(gene_list$adj.P.Val), colour=gene_list$threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-5, 5)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle("Volcano Plot of Fold Change vs P-Value (0.01 Colour Change)")
g
###END GGPLOT2
#dev.off()

###VENN
##Use for Comparison of Multiple array classesw
results = classifyTestsP(fit2, p.value = 0.01, method = "fdr")
vennDiagram(results[, c("NOF - OA", "OA - NOF")])
###END VENN

gene_list_complete = topTable(fit2, coef = "NOF - OA", p.value=0.01, lfc=log2(2), number = nrow(anno_df))
###END OF DIFFERENTIAL EXPRESSION ANALYSIS

###FUNCTIONAL ANALYSIS

##HIERARCHICAL CLUSTERING
fa_data          <- exprs(analysis_ready_data_post)
dist_measure     <- dist(t(exprs(analysis_ready_data_post)), method="euclidian")
clust_by_samples <- hclust(dist_measure)
plot(clust_by_samples)

kmeans_by_samples <- kmeans(t(exprs(analysis_ready_data_post)), centers=8)
print(kmeans_by_samples$cluster)

##EXTRACT NOF-OA DIFFERENTIALLY EXPRESSED GENE DATA
nof_oa <- topTable(fit2, coef="NOF - OA", p.value=0.01,
                   lfc=log2(2), number=nrow(anno_df))
write.csv(file="Differentially_Expressed_Genes.csv", nof_oa)

##HEATMAP OF TOP 100 GENES
diff_exp <- fa_data[rownames(nof_oa),1:length(colnames(fa_data))]
heatmap.2(diff_exp)
heatmap.2(diff_exp, col="redblue", density.info="none",
          trace="none", cexRow=0.2)
dev.off()
###END OF FUNCTIONAL ANALYSIS

###RAMI_GO
sig_probes      <- as.character(nof_oa$PROBE_ID)
entrez<- unique(unlist(lookUp(nuIDs[sig_probes],
                              "lumiHumanAll.db", "ENTREZID")))
entrez          <- as.character(entrez[!is.na(entrez)])
entrez_universe <- unique(unlist(lookUp(nuIDs,"lumiHumanAll.db", "ENTREZID")))
entrez_universe <- as.character(entrez_universe[!is.na(entrez_universe)])

params <- new("GOHyperGParams", geneIds=entrez,universeGeneIds=entrez_universe,annotation="lumiHumanAll.db",ontology="BP",
              pvalueCutoff= 0.01, conditional=FALSE, testDirection="over")
hyperg_result <- hyperGTest(params)
print(hyperg_result)
str(summary(hyperg_result))

pval_go     <- pvalues(hyperg_result)
go_fdr      <- p.adjust(pvalues(hyperg_result), method="fdr")
sig_go_id   <- names(go_fdr[go_fdr < 0.01])
sig_go_term <- getGOTerm(sig_go_id)[["BP"]]

top_5_go = sig_go_id[1:5]
top_5_p = go_fdr[go_fdr < 0.01][1:5]
amigo_tree = getAmigoTree(top_5_go, pvalues=top_5_p, filename='RamiGo.png')
###END RAMI_GO

###KEGG
kegg_params <- new("KEGGHyperGParams",
                   geneIds=entrez,
                   universeGeneIds=entrez_universe,
                   annotation="lumiHumanAll.db",
                   pvalueCutoff= 0.01,
                   testDirection="over")
kegg_hyperg_result <- hyperGTest(kegg_params)
print(kegg_hyperg_result)

sig_logFCs <- nof_oa$logFC
fc_entrez <- unlist(lookUp(nuIDs[sig_probes],
                   "lumiHumanAll.db", "ENTREZID"))
names(sig_logFCs) <- fc_entrez
sig_logFCs <- sig_logFCs[!is.na(entrez)]

for(i in 1:length(sigCategories(kegg_hyperg_result)))
{
  pv <- pathview(gene.data=sig_logFCs,
                 pathway.id=sigCategories(kegg_hyperg_result)[i],
                 species="hsa",
                 limit=list(gene=c(round(min(sig_logFCs)),
                                   round(max(sig_logFCs))), cpd=1))
}
###END KEGG