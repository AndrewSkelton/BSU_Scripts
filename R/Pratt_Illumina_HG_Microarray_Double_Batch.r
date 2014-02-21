setwd("/data/customers/Musculoskeletal/Pratt/Rheumatoid_Arthritis_CD4_T-Cell_Signature//Raw_Data")
list.files()
raw_files <- c("Raw Data Phase I.txt", "Raw Data Phase II.txt")

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("sva")
library(limma)
library(sva)

###RAW DATASET
raw_data            <- lumiR.batch(raw_files)
removal_elements     <- c(87)
removal              <- sampleNames(raw_data)[removal_elements]
raw_data             <- raw_data[, !(sampleNames(raw_data) %in% removal)]
vst_data            <- lumiT(raw_data, method="vst")
rsn_data            <- lumiN(vst_data, method = "rsn")
analysis_ready_data <- lumiQ(rsn_data)
#arrayQualityMetrics(expressionset=analysis_ready_data, outdir="QC_INITIAL_DATASET")

new_Array_Names <- c("177_D_I_i", "304_C_I_i", "182_B_I_i", "174_A_I_i", "239_A_I_i", "264_B_I_i", "151_C_I_i", "175_D_I_i", "259_D_I_i", 
                     "314_C_I_i", "214_B_I_i", "322_A_I_i", "273_A_I_i", "315_B_I_i", "307_C_I_i", "277_D_I_i", "192_D_I_i", "193_C_I_i", 
                     "236_B_I_i", "184_B_I_i", "241_A_I_i", "221_B_I_i", "263_C_I_i", "270_D_I_ii", "210_C_I_ii", "153_C_I_ii", "195_B_I_ii", 
                     "253_A_I_ii", "254_A_I_ii", "294_B_I_ii", "231_C_I_ii", "230_D_I_ii", "248_D_I_ii", "191_C_I_ii", "218_B_I_ii", 
                     "196_A_I_ii", "246_B_I_ii", "149_C_I_ii", "237_C_I_ii", "156_D_I_ii", "260_D_I_ii", "181_C_I_ii", "256_B_I_ii", "220_A_I_ii", 
                     "305_A_I_ii", "133_B_I_ii", "272_C_I_ii", "157_D_I_iii", "291_D_I_iii", "240_C_I_iii", "268_B_I_iii", "302_A_I_iii", 
                     "243_C_I_iii", "125_C_I_iii", "188_C_I_iii", "247_B_I_iii", "269_B_I_iii", "209_A_I_iii", "142_A_I_iii", "261_B_I_iii", 
                     "303_C_I_iii", "211_B_I_iii", "176_A_I_iii", "251_B_I_iii", "250_C_I_iii", "306_B_I_iii", "292_D_I_iii", "172_B_I_iii", 
                     "134_C_I_iii", "137_C_I_iii", "204_B_I_iii", "266_C_I_iii", "183_B_I_iv", "120_A_I_iv", "144_C_I_iv", "145_C_I_iv", 
                     "207_D_I_iv", "141_D_I_iv", "118_C_I_iv", "249_B_I_iv", "203_B_I_iv", "131_A_I_iv", "299_A_I_iv", "163_B_I_iv", "170_C_I_iv", 
                     "159_D_I_iv", "161_D_I_iv", "136_C_I_iv", "300_D_I_iv", "232_A_I_iv", "288_B_I_iv", "296_A_I_iv", "212_D_I_iv", "132_C_I_iv", 
                     "135_C_I_iv", "388_A_II_v", "297_C_II_v", "154_C_II_v", "333_B_II_v", "364_A_II_v", "316_C_II_v", "336_D_II_v", "227_E_II_v", 
                     "442_C_II_v", "358_B_II_v", "155_E_II_v", "382_E_II_v", "284_E_II_v", "355_E_II_v", "242_E_II_v", "337_E_II_v", "205_E_II_v", 
                     "286_E_II_v", "165_D_II_v", "301_B_II_v", "258_D_II_v", "430_E_II_v", "323_D_II_v", "401_A_II_v", "349_B_II_v", "441_D_II_v", 
                     "391_C_II_v", "431_C_II_v", "354_C_II_v", "126_D_II_v", "213_C_II_v", "262_E_II_v", "238_C_II_v", "197_B_II_v", "367_E_II_v", 
                     "310_E_II_v", "255_E_II_vi", "311_C_II_vi", "397_A_II_vi", "283_D_II_vi", "392_D_II_vi", "366_A_II_vi", "353_C_II_vi", 
                     "346_D_II_vi", "309_E_II_vi", "287_B_II_vi", "282_C_II_vi", "334_E_II_vi", "344_A_II_v", "313_C_II_v", "340_D_II_v", 
                     "200_E_II_v", "387_C_II_v", "406_C_II_v", "433_D_II_v", "276_A_II_v", "389_A_II_v", "369_A_II_v", "359_B_II_v", "285_E_II_v",
                     "293_C_II_v", "244_A_II_v", "393_B_II_v", "439_A_II_v", "326_E_II_v", "412_C_II_v", "350_A_II_v", "405_C_II_v", "198_D_II_v", 
                     "351_C_II_v", "320_C_II_v", "378_B_II_v", "223_C_II_v", "374_B_II_v", "173_E_II_v", "328_D_II_v", "216_D_II_v", "222_B_II_v")
sampleNames(analysis_ready_data)      <- new_Array_Names

###DETECTION CALL THRESHOLD REMOVAL
exprs_data                   <- exprs(analysis_ready_data)
present_count                <- detectionCall(analysis_ready_data) 
filtered_analysis_ready_data <- exprs_data[present_count > 0, ]

##'BOTH ARRAYS NORMALISED RELATIVE TO EACH OTHER
phase_1 <- filtered_analysis_ready_data[,1:95]
phase_2 <- filtered_analysis_ready_data[,96:173]
batches <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
batches_phase_1 <- batches[1 :95 ]
batches_phase_2 <- batches[96:173]
batches_phase_2 <- gsub(5, 1, batches_phase_2); batches_phase_2 <- gsub(6, 2, batches_phase_2)

##'CHECK FOR BATCH EFFECT NOW (AMPLIFICATION BATCHES)
pdf("Post_Amplification_Correction_Plots.pdf")
PCA_data <- phase_1
colnames(PCA_data) <- batches_phase_1
#Pre Amplification Correction
#Post Phase Correction
plotSampleRelation(PCA_data, method= 'cluster' , cv.Th=0, main="Pre Amplification Correction - Phase 1", cex=.4)
plotSampleRelation(PCA_data, method= 'mds' , cv.Th=0, main="Pre Amplification Correction - Phase 1")
dev.off()
####

treatments  <- c("A", "B", "C", "D", "E")
array_names <- c("D", "C", "B", "A", "A", "B", "C", "D", "D", "C", "B", "A", "A", "B", "C", "D", "D", "C", "B", "B", "A", "B", "C", "D", "C", "C", "B", "A", "A", "B", "C", "D", "D", "C", "B", "A", "B", "C", "C", "D", "D", "C", "B", "A", "A", "B", "C", "D", "D", "C", "B", "A", "C", "C", "C", "B", "B", "A", "A", "B", "C", "B", "A", "B", "C", "B", "D", "B", "C", "C", "B", "C", "B", "A", "C", "C", "D", "D", "C", "B", "B", "A", "A", "B", "C", "D", "D", "C", "D", "A", "B", "A", "D", "C", "C", "A", "C", "C", "B", "A", "C", "D", "E", "C", "B", "E", "E", "E", "E", "E", "E", "E", "E", "D", "B", "D", "E", "D", "A", "B", "D", "C", "C", "C", "D", "C", "E", "C", "B", "E", "E", "E", "C", "A", "D", "D", "A", "C", "D", "E", "B", "C", "E", "A", "C", "D", "E", "C", "C", "D", "A", "A", "A", "B", "E", "C", "A", "B", "A", "E", "C", "A", "C", "D", "C", "C", "B", "C", "B", "E", "D", "D", "B")

##'BATCH CORRECTION (AMPLIFICATION)
pheno            <- data.frame(sample = c(1:95), outcome = array_names[1:95], batch = batches_phase_1)
rownames(pheno)  <- colnames(phase_1)
batch = pheno$batch
mod = model.matrix(~as.factor(outcome), data=pheno)
phase_1_combat_edata = ComBat(dat=phase_1, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)

pheno            <- data.frame(sample = c(1:78), outcome = array_names[1:78], batch = batches_phase_2)
rownames(pheno)  <- colnames(phase_2)
batch = pheno$batch
mod = model.matrix(~as.factor(outcome), data=pheno)
phase_2_combat_edata = ComBat(dat=phase_2, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
##'END BATCH CORRECTION (AMPLIFICATION)

merged_data_post_amp <- cbind(phase_1_combat_edata, phase_2_combat_edata)











###PERCENTAGE OF PROBES REMOVED [40.519%]
(nrow(filtered_analysis_ready_data)/nrow(analysis_ready_data))*100
###END OF DETECTION CALL THRESHOLD REMOVAL

##'BATCH EFFECT - HEATMAP
library(gplots)
batch_map_names <- c("I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", 
                     "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", 
                     "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", 
                     "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", 
                     "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", 
                     "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", 
                     "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", 
                     "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II", "II")

batch_map_names <- c("i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", "i", 
                     "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii", 
                     "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii", "ii", 
                     "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii", 
                     "iii", "iii", "iii", "iii", "iii", "iii", "iii", "iii", 
                     "iv", "iv", "iv", "iv", "iv", "iv", "iv", "iv", "iv", "iv", "iv", "iv", "iv", "iv", 
                     "iv", "iv", "iv", "iv", "iv", "iv", "iv", "iv", "iv", 
                     "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", 
                     "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", 
                     "vi", "vi", "vi", "vi", "vi", "vi", "vi", "vi", "vi", "vi", "vi", "vi", 
                     "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", 
                     "v", "v", "v", "v", "v", "v", "v", "v")

batch_map <- filtered_analysis_ready_data[1:1000,]
colnames(batch_map) <- new_Array_Names#batch_map_names
heatmap.2(batch_map, col="redblue", density.info="none",
          trace="none", cexRow=0.2, main="Batch Effect")

PCA_data <- as.matrix(combat_edata)#combat_edata

pdf("Post_Amplification_Correction_Plots.pdf")
colnames(PCA_data) <- batch_map_names
#Pre Amplification Correction
#Post Phase Correction
plotSampleRelation(PCA_data, method= 'cluster' , cv.Th=0, main="Post Amplification Correction", cex=.4)
plotSampleRelation(PCA_data, method= 'mds' , cv.Th=0, main="Post Amplification Correction")
dev.off()
####

###ARRAY ANNOTATION
probe_list <- rownames(filtered_analysis_ready_data)
probe_list <- probe_list[-match(3450064, probe_list)]
nuIDs      <- probeID2nuID(probe_list)[, "nuID"]
symbol     <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name       <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df    <- data.frame(ID = nuIDs, probe_list, symbol, name)
###END OF ARRAY ANNOTATION



treatments  <- c("A", "B", "C", "D", "E")
array_names <- c("D", "C", "B", "A", "A", "B", "C", "D", "D", "C", "B", "A", "A", "B", "C", "D", "D", "C", "B", "B", "A", "B", "C", "D", "C", "C", "B", "A", "A", "B", "C", "D", "D", "C", "B", "A", "B", "C", "C", "D", "D", "C", "B", "A", "A", "B", "C", "D", "D", "C", "B", "A", "C", "C", "C", "B", "B", "A", "A", "B", "C", "B", "A", "B", "C", "B", "D", "B", "C", "C", "B", "C", "B", "A", "C", "C", "D", "D", "C", "B", "B", "A", "A", "B", "C", "D", "D", "C", "D", "A", "B", "A", "D", "C", "C", "A", "C", "C", "B", "A", "C", "D", "E", "C", "B", "E", "E", "E", "E", "E", "E", "E", "E", "D", "B", "D", "E", "D", "A", "B", "D", "C", "C", "C", "D", "C", "E", "C", "B", "E", "E", "E", "C", "A", "D", "D", "A", "C", "D", "E", "B", "C", "E", "A", "C", "D", "E", "C", "C", "D", "A", "A", "A", "B", "E", "C", "A", "B", "A", "E", "C", "A", "C", "D", "C", "C", "B", "C", "B", "E", "D", "D", "B")
design <- model.matrix(~0 + factor(array_names, levels = treatments))

batches <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
batches <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)

combat_expr_data = ComBat(dat=analysis_ready_data, batch=batches, mod=design, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)