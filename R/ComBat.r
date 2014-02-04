pheno            <- data.frame(sample = c(1:173), outcome = array_names, batch = batches)
rownames(pheno)  <- colnames(combat_edata)
batch = pheno$batch
mod = model.matrix(~as.factor(outcome), data=pheno)
combat_edata = ComBat(dat=combat_edata, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)

plotSampleRelation(combat_edata, method= 'mds' , cv.Th=0, main="ComBat")
