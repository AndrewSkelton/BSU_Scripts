setwd("/data/customers/Musculoskeletal/Williams/Leptin_IL1/Raw_data")
list.files()

##'RAW DATASET NORMALISATION
raw_data            <- lumiR("raw_data.txt")
vst_data            <- lumiT(raw_data, method = "vst")
rsn_data            <- lumiN(vst_data, method = "rsn")
analysis_ready_data <- lumiQ(rsn_data)
arrayQualityMetrics(expressionset=analysis_ready_data, outdir="QC_INITIAL_DATASET")

##'VISUALISATION OF NORMALISATION
pdf("Normalisation_Plots.pdf")
plot(raw_data,  main="Density Plot of Raw Dataset")
boxplot(raw_data,  main="Boxplot of Raw Dataset")
plot(vst_data,  main="Density Plot of Raw Dataset - Post VST")
boxplot(vst_data,  main="Boxplot of Raw Dataset - Post VST")
plot(rsn_data,  main="Density Plot of Raw Dataset - Post RSN")
boxplot(rsn_data,  main="Boxplot of Raw Dataset - Post RSN")
dev.off()

analysis_ready_data_post <- analysis_ready_data
##'END OF NORMALISATION

###DETECTION CALL THRESHOLD REMOVAL
exprs_data                   <- exprs(analysis_ready_data_post)
present_count                <- detectionCall(analysis_ready_data_post) 
filtered_analysis_ready_data <- exprs_data[present_count > 0, ]
###PERCENTAGE OF PROBES REMOVED [47.04142%]
(nrow(filtered_analysis_ready_data)/nrow(analysis_ready_data_post))*100
###END OF DETECTION CALL THRESHOLD REMOVAL

###ARRAY ANNOTATION
probe_list <- rownames(analysis_ready_data)
nuIDs      <- probeID2nuID(probe_list)[, "nuID"]
symbol     <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name       <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df    <- data.frame(ID = nuIDs, probe_list, symbol, name)
###END OF ARRAY ANNOTATION

##'BATCH CORRECTION
pdf("PCA_Plots.pdf")
plot(analysis_ready_data, what="sampleRelation", method="mds", main="2D PCA of Samples")
batches <- c(1,2,3,1,2,3,1,2,3,1,2,3)
pheno            <- data.frame(sample = c(1:length(batches)), outcome = array_names, batch = batches)
rownames(pheno)  <- colnames(filtered_analysis_ready_data)
batch = pheno$batch
mod = model.matrix(~as.factor(outcome), data=pheno)
combat_edata = ComBat(dat=filtered_analysis_ready_data, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
plotSampleRelation(combat_edata, method= 'mds' , cv.Th=0, main="2D PCA of Samples - Batch Corrected")
dev.off()
##'

###EXPERIMENTAL SETUP
treatments  <- c("Control", "Leptin", "IL1", "Leptin_IL1")
array_names <- c("Control","Control","Control","Leptin","Leptin","Leptin",
                 "IL1","IL1","IL1","Leptin_IL1","Leptin_IL1","Leptin_IL1")

design <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design) <- treatments
num_parameters <- ncol(design)
fit <- lmFit(combat_edata, design)

cont_mat <- makeContrasts(Control-Leptin, Control-IL1, Leptin_IL1-Leptin, 
                          Leptin_IL1-IL1, Leptin_IL1-Control, levels=treatments)
fit2 <- contrasts.fit(fit, contrasts=cont_mat)
fit2 <- eBayes(fit2)
fit2$genes = anno_df
##'

##'FILTERING
##'Filter by fold change (2x) and p-value (<0.01) cutoffs
##'Adjusted using Benjamini-Hochberg False Discovery Rate
all_genes <- topTable(fit2, p.value=0.01, lfc=log2(2), adjust.method="BH")

#CONTROL VS LEPTIN
topTable(fit2, coef="Control - Leptin", p.value=0.01, lfc=log2(2), adjust.method="BH")
gene_list_Control_Leptin_unfiltered         <- topTable(fit2, coef = "Control - Leptin", number = nrow(anno_df), adjust.method="BH")
gene_list_Control_Leptin                    <- topTable(fit2, coef = "Control - Leptin", p.value=0.01, lfc=log2(2), number = nrow(anno_df), 
                                                        adjust.method="BH")

#CONTROL VS IL-1
topTable(fit2, coef="Control - IL1", p.value=0.01, lfc=log2(2), adjust.method="BH")
gene_list_Control_IL1_unfiltered         <- topTable(fit2, coef = "Control - IL1", number = nrow(anno_df), adjust.method="BH")
gene_list_Control_IL1                    <- topTable(fit2, coef = "Control - IL1", p.value=0.01, lfc=log2(2), number = nrow(anno_df), 
                                                        adjust.method="BH")

#LEPTIN_IL-1 VS LEPTIN
topTable(fit2, coef="Leptin_IL1 - Leptin", p.value=0.01, lfc=log2(2), adjust.method="BH")
gene_list_LeptinIL1_Leptin_unfiltered         <- topTable(fit2, coef = "Leptin_IL1 - Leptin", number = nrow(anno_df), adjust.method="BH")
gene_list_LeptinIL1_Leptin                    <- topTable(fit2, coef = "Leptin_IL1 - Leptin", p.value=0.01, lfc=log2(2), number = nrow(anno_df), 
                                                        adjust.method="BH")

#LEPTIN_IL-1 VS IL-1
topTable(fit2, coef="Leptin_IL1 - IL1", p.value=0.01, lfc=log2(2), adjust.method="BH")
gene_list_LeptinIL1_IL1_unfiltered         <- topTable(fit2, coef = "Leptin_IL1 - IL1", number = nrow(anno_df), adjust.method="BH")
gene_list_LeptinIL1_IL1                    <- topTable(fit2, coef = "Leptin_IL1 - IL1", p.value=0.01, lfc=log2(2), number = nrow(anno_df), 
                                                        adjust.method="BH")

#LEPTIN_IL-1 VS CONTROL
topTable(fit2, coef="Leptin_IL1 - Control", p.value=0.01, lfc=log2(2), adjust.method="BH")
gene_list_LeptinIL1_Control_unfiltered         <- topTable(fit2, coef = "Leptin_IL1 - Control", number = nrow(anno_df), adjust.method="BH")
gene_list_LeptinIL1_Control                    <- topTable(fit2, coef = "Leptin_IL1 - Control", p.value=0.01, lfc=log2(2), number = nrow(anno_df), 
                                                        adjust.method="BH")

##'WRITE TEXT FILES
write.table(gene_list_Control_Leptin,    "Control_Leptin.txt",    sep="\t")
write.table(gene_list_Control_IL1,       "Control_IL1.txt",       sep="\t")
write.table(gene_list_LeptinIL1_Leptin,  "LeptinIL1_Leptin.txt",  sep="\t")
write.table(gene_list_LeptinIL1_IL1,     "LeptinIL1_IL1.txt",     sep="\t")
write.table(gene_list_LeptinIL1_Control, "LeptinIL1_Control.txt", sep="\t")
##

##HIERARCHICAL CLUSTERING
fa_data          <- combat_edata
dist_measure     <- dist(t(combat_edata), method="euclidian")
clust_by_samples <- hclust(dist_measure)
clust_by_samples$labels <- array_names
pdf("Dendogram.pdf")
plot(clust_by_samples, cex=.5)
dev.off()
##

##'VOLCANO PLOTS
library(ggplot2)
gene_lists <- c("gene_list_Control_Leptin_unfiltered", "gene_list_Control_IL1_unfiltered", "gene_list_LeptinIL1_Leptin_unfiltered", 
                "gene_list_LeptinIL1_IL1_unfiltered", "gene_list_LeptinIL1_Control_unfiltered")

pdf("Volcano_Plots.pdf")
for(i in 1:length(gene_lists))
{
  print(i)
  gene_list <- eval(parse(text=gene_lists[i]))  
  
  gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.01/nrow(gene_list))#/nrow(anno_df))
  g = ggplot(data=gene_list, aes(x=gene_list$logFC, y=-log10(gene_list$adj.P.Val), colour=gene_list$threshold)) +
    geom_point(alpha=0.4, size=1.75) +
    labs(title = gene_lists[i]) + 
    theme(legend.position = "none") +
    xlim(c(-8, 8)) + ylim(c(0, 11)) +
    xlab("log2 fold change") + ylab("-log10 p-value")
  print(g)
}
dev.off()
##

##'VENN
pdf("Venn.pdf")
results = classifyTestsP(fit2, p.value = 0.01, method = "fdr")
coeff_comparison <- c("Control - Leptin", "Control - IL1", "Leptin_IL1 - Leptin", 
                      "Leptin_IL1 - IL1", "Leptin_IL1 - Control")
vennDiagram(results[, coeff_comparison])
dev.off()
##

##'HEATMAPS OF DIFF EXP GENES
library(gplots)
gene_lists_1 <- c("gene_list_Control_Leptin", "gene_list_Control_IL1", "gene_list_LeptinIL1_Leptin", 
                  "gene_list_LeptinIL1_IL1", "gene_list_LeptinIL1_Control")

pdf("Heatmaps.pdf")
for(i in 1:length(gene_lists_1))
{
  print(i)
  gene_list1 <- eval(parse(text=gene_lists_1[5]))  
  
  diff_exp <- fa_data[rownames(gene_list1),c(10,11,12,1,2,3)]
  heatmap.2(diff_exp)
  heatmap.2(diff_exp, col="redblue", density.info="none",
            trace="none", cexRow=0.2, main=gene_lists_1[i])
}
dev.off()
##

###RAMI_GO
library(lumiHumanAll.db)
library(annotate)
library("GOstats")
library("RamiGO")

for(i in 1:length(gene_lists_1))
  #for(i in 1:5)
  {
    print(i)
    gene_list1      <- eval(parse(text=gene_lists_1[i]))  
    sig_probes      <- as.character(gene_list1$probe_list)
    entrez          <- unique(unlist(lookUp(nuIDs[sig_probes],
                                            "lumiHumanAll.db", "ENTREZID")))
    entrez          <- as.character(entrez[!is.na(entrez)])
    entrez_universe <- unique(unlist(lookUp(nuIDs,"lumiHumanAll.db", "ENTREZID")))
    entrez_universe <- as.character(entrez_universe[!is.na(entrez_universe)])
    
    params <- new("GOHyperGParams", geneIds=entrez,universeGeneIds=entrez_universe,annotation="lumiHumanAll.db",ontology="BP",
                  pvalueCutoff= 0.05, conditional=FALSE, testDirection="over")
    hyperg_result <- hyperGTest(params)
    print(hyperg_result)
    str(summary(hyperg_result))
    
    pval_go     <- pvalues(hyperg_result)
    go_fdr      <- p.adjust(pvalues(hyperg_result), method="fdr")
    sig_go_id   <- names(go_fdr[go_fdr < 0.05])
    sig_go_term <- getGOTerm(sig_go_id)[["BP"]]
    
    if(length(sig_go_id) > 0)
    {
      print(length(sig_go_id))
      top_5_go = sig_go_id
      top_5_p = go_fdr[go_fdr < 0.05]
      amigo_tree = getAmigoTree(top_5_go, pvalues=top_5_p, filename=paste(gene_lists_1[i], "__RamiGo.png", sep=" "))
    }
  }
###END RAMI_GO

###KEGG
library('pathview')
parentDir <- getwd()
for(i in 1:length(gene_lists_1))
{
  print(i)
  gene_list1      <- eval(parse(text=gene_lists_1[i]))
  
  subDir <- paste(gene_lists_1[i], "_pathview", sep="")
  dir.create(file.path(parentDir, subDir))
  setwd(file.path(parentDir, subDir))
  
  sig_probes      <- as.character(gene_list1$probe_list)
  entrez          <- unique(unlist(lookUp(nuIDs[sig_probes],
                                          "lumiHumanAll.db", "ENTREZID")))
  entrez          <- as.character(entrez[!is.na(entrez)])
  entrez_universe <- unique(unlist(lookUp(nuIDs,"lumiHumanAll.db", "ENTREZID")))
  entrez_universe <- as.character(entrez_universe[!is.na(entrez_universe)])
  
  kegg_params <- new("KEGGHyperGParams",
                     geneIds=entrez,
                     universeGeneIds=entrez_universe,
                     annotation="lumiHumanAll.db",
                     pvalueCutoff= 0.01,
                     testDirection="over")
  kegg_hyperg_result <- hyperGTest(kegg_params)
  print(kegg_hyperg_result)
  
  sig_logFCs <- gene_list1$logFC
  fc_entrez <- unlist(lookUp(nuIDs[sig_probes],
                             "lumiHumanAll.db", "ENTREZID"))
  names(sig_logFCs) <- fc_entrez
  sig_logFCs <- sig_logFCs[!is.na(entrez)]
  
  if(length(sigCategories(kegg_hyperg_result)) < 6)
  {
    end_cond <- length(sigCategories(kegg_hyperg_result))
  } else
  {
    end_cond <- 6
  }
  
  for(j in 1:end_cond)
  {
    pv <- pathview(gene.data=sig_logFCs,
                   pathway.id=sigCategories(kegg_hyperg_result)[j],
                   species="hsa",
                   limit=list(gene=c(round(min(sig_logFCs)),
                                     round(max(sig_logFCs))), cpd=1))
  }
  setwd(parentDir)
}
###END KEGG
