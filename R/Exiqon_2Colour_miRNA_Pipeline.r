##'LOAD IN GAL FILE (ANNOTATION) AND RAW DATA
biocLite("ExiMiR")
library(ExiMiR)
install.packages("miRada")
make.gal.env(galname="gal_208500,208501,208502,208510_lot35003-35003_hsa-and-related-vira_from_mb180,miRPlus.gal" , gal.path="./" )
ebatch <- ReadExi(galname="gal_208500,208501,208502,208510_lot35003-35003_hsa-and-related-vira_from_mb180,miRPlus.gal" , txtfile.path="./" )
##

##'LIMMA RAW OBJECT CREATION
library(limma)
targets        <- readTargets(path = "./" )
RGList         <- read.maimages(targets[,c("Cy3", "Cy5")], source = "imagene", path = "./" )
RGList$genes   <- readGAL(path = "./")
RGList$printer <- getLayout(RGList$genes)

###AFFY NORMALISATION
obatch         <- createAB(RGList)
spikein.all    <- grep( '^spike', featureNames(obatch), value=TRUE)
spikein.subset <- setdiff(spikein.all, 'spike_control_f')
spikein.params <- list(probeset.list=spikein.subset,
                       loess.span=0.6,
                       force.zero=TRUE,
                       figures.show=TRUE, figures.output=file)
eset.spike <- NormiR(ebatch, bgcorrect.method= 'normexp',
                     bgcorrect.param=list(offset=50),
                     normalize.method= 'spikein',
                     normalize.param= list(spikein.params),
                     summary.method= 'medianpolish')
###
##


##'BOXPLOT OF INTENSITIES FOR EACH PAIRED ARRAY
##'FROM PROVIDED TARGETS TEXT FILE
pdf("Boxplot_raw_intensities.pdf")
labs <- colnames(data.frame(cbind(log2(RGList$Gb),log2(RGList$Rb))))
lab_test <- c("IL-17", "IL-17/IFNy", "IL-17/IFNy", "IL-17", "IL-17",
              "IFNy", "IFNy", "IL-17/IFNy", "IFNy", "IL-17", "IFNy",	"IL-17/IFNy")
lab_test_Cy3 <- paste(lab_test, "_Cy3", sep="")
lab_test_Cy5 <- paste(lab_test, "_Cy5", sep="")
lab_test_final <- c(lab_test_Cy3, lab_test_Cy5)
boxplot(data.frame(cbind(log2(RGList$Gb),log2(RGList$Rb))),
              main="Boxplot of Channel Intensity Values", xaxt="n",
              col=rep(c("green","red"),each=12))
axis(1, at=seq(1, length(lab_test_final), by=1), labels = FALSE)
text(seq(1+.75, length(lab_test_final)+.75, by=1), par("usr")[3]-.1, 
     labels = lab_test_final, srt = -45, pos = 1, xpd = TRUE, cex=0.7)
dev.off()
##

##'R' CAN BE CHANGED FOR OTHER ATTRIBUTES
##'R  - RED CHANNEL
##'G  - GREEN CHANNEL
##'Rb - RED BACKGROUND CHANNEL
##'Gb - GREEN BACKGROUND CHANNEL
pdf("Pseudo_Array_Images.pdf")
for(i in 1:length(colnames(RGList$R)))
{
  imageplot(log2(RGList$R[,i]),RGList$printer, main=paste("Array",i ,"Pseudo-Image - Red Channel", sep=" "))
}
dev.off()
##

##'SCATTER PLOT OF SIGNAL VS BACKGROUND
pdf("Scatter_plots_primary_channel_vs_background_channel.pdf")
for(i in 1:length(colnames(RGList$G)))
{
  plot(log2(RGList$Gb[,i]),log2(RGList$G[,i]),
       main=paste("Scatter Plot of Array", i, "- Green Channel Signal Vs Background"))
  lines(c(-9,99),c(-9,99),col=2)
}
for(i in 1:length(colnames(RGList$R)))
{
  plot(log2(RGList$Rb[,i]),log2(RGList$R[,i]),
       main=paste("Scatter Plot of Array", i, "- Red Channel Signal Vs Background"))
  lines(c(-9,99),c(-9,99),col=2)
}
dev.off()

##'BACKGROUND CORRECTION
RGList <- backgroundCorrect(RGList, method = "normexp")

labs <- colnames(data.frame(cbind(log2(RGList$Gb),log2(RGList$Rb))))
lab_test <- c("IL-17", "IL-17/IFNy", "IL-17/IFNy", "IL-17", "IL-17",
              "IFNy", "IFNy", "IL-17/IFNy", "IFNy", "IL-17", "IFNy",  "IL-17/IFNy")
lab_test_Cy3 <- paste(lab_test, "_Cy3", sep="")
lab_test_Cy5 <- paste(lab_test, "_Cy5", sep="")
for(i in 1:12)
{
  lab_test_final_1 <- c(lab_test_final_1, lab_test_Cy5[i], lab_test_Cy3[i])
}
pdf("Boxplot_Background_Corrected.pdf")
boxplot(log2(data.frame(
  RGList$R,RGList$G))[,as.vector(rbind(1:12,13:24))],
        main="Boxplot of Background Corrected RG Channel Intensity Values", xaxt="n",
        col=rep(c("red","green"),24))
axis(1, at=seq(1, length(lab_test_final_1), by=1), labels = FALSE)
text(seq(1+.75, length(lab_test_final_1)+.75, by=1), par("usr")[3]-1.5, 
     labels = lab_test_final_1, srt = -45, pos = 1, xpd = TRUE, cex=0.7)
dev.off()
##

##'BACKGROUND CORRECTED DENSITY
pdf("Background_Corrected_Densities.pdf")
plotDensities(RGList)
dev.off()
##

##'MA PLOTS OF THE DIFFERENCE VS AVERAGE LOG BACKGROUND
##'CORRECTED SIGNAL
pdf("MA_Plots_Background_Corrected.pdf")
for(i in 1:length(lab_test))
{
  plotMA(RGList[,i], main=paste(lab_test[i], "- MA Plot"))
}
dev.off()
##

##'LOESS AND MEDIAN NORMALISATION
##'SPIKE NOT AVAILABLE AS METHOD
MAList <- normalizeWithinArrays(RGList, method = "loess")
MAList <- normalizeWithinArrays(MAList, method = "median")
##

##'POST NORMALISATION PLOTS
pdf("Normalised_Plots_Comparison.pdf")
for(i in 1:12)
{
  plotMA(RGList[,i])
  plotMA(MAList[,i])
}
plotDensities(RGList, main="Unnormalised Densities")
plotDensities(MAList, main="Normalised Densities")

boxplot(data.frame(MAList$M),
  main="Boxplot of Normalised Intensity Values", xaxt="n")
axis(1, at=seq(1, length(lab_test), by=1), labels = FALSE)
text(seq(1+.25, length(lab_test)+.25, by=1), par("usr")[3]-1.2, 
     labels = lab_test, srt = -45, pos = 1, xpd = TRUE, cex=0.7)
dev.off()
##

##'SCALE NORMALISATION
MAList <- normalizeBetweenArrays(MAList, method = "scale")
pdf("Boxplot_post_scale_normalisation.pdf")
boxplot(data.frame(MAList$M),
  main="Boxplot of Normalised (Post Scale Normalisation) Intensity Values", xaxt="n")
axis(1, at=seq(1, length(lab_test), by=1), labels = FALSE)
text(seq(1+.25, length(lab_test)+.25, by=1), par("usr")[3]-1.2, 
     labels = lab_test, srt = -45, pos = 1, xpd = TRUE, cex=0.7)
dev.off()        

pdf("Fully_Normalised_MA_Plots_Density.pdf")
for(i in 1:12)
{
  plotMA(MAList[,i])
}
plotDensities(RGList, main="Unnormalised Densities")
plotDensities(MAList, main="Normalised Densities")
dev.off()
##

##'DIFFERENTIAL EXPRESSION
treatments       <- c("IL17", "IL17_IFNy", "IFNy")
array_names      <- c("IL17", "IL17_IFNy", "IL17_IFNy", "IL17", "IL17", "IFNy",
                      "IFNy", "IL17_IFNy", "IFNy", "IL17", "IFNy", "IL17_IFNy")
design           <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design) <- treatments
num_parameters   <- ncol(design)
fit              <- lmFit(MAList, design)

cont_mat <- makeContrasts(IFNy-IL17, IL17_IFNy-IL17, IL17_IFNy-IFNy,  levels=treatments)
fit2 <- contrasts.fit(fit, contrasts=cont_mat)
fit2 <- eBayes(fit2)


gene_list_IFNy_IL17     <- topTable(fit2, coef="IFNy - IL17", p.value=0.05, lfc=log2(1.5), adjust.method="BH")
gene_list_IL17IFNy_IL17 <- topTable(fit2, coef="IL17_IFNy - IL17", p.value=0.05, lfc=log2(1.5), adjust.method="BH")
gene_list_IL17IFNy_IFNy <- topTable(fit2, coef="IL17_IFNy - IFNy", p.value=0.05, lfc=log2(1.5), adjust.method="BH")