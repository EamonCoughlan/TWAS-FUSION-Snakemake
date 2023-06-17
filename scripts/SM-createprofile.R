#script was not working, I think cv.performance requires a library. Trying the libraries from compute_weights script
#no, cv.performance is in the Rdat file. Maybe it's an issue with everything being set to look for blup. Try adding blup to models.

#take rows from input file as a list:
models <- readLines(snakemake@input[[1]])
profile<-matrix(0,length(models), 11)
colnames(profile)<-c("id", "nsnps", "hsq", "hsq.se", "hsq.pv", "top1.r2", "enet.r2", "lasso.r2", "top1.pv", "enet.pv", "lasso.pv")
profile[,1]<-gsub("models/","", models)
profile[,1]<-gsub(snakemake@wildcards['condition'],"", profile[,1], fixed=T)
profile[,1]<-gsub("/","", profile[,1], fixed=T)
profile[,1]<-gsub(".wgt.RDat","", profile[,1], fixed=T)
rownames(profile)<-models

for(i in models){
	#try to load the model, but if it fails, skip it
	tryCatch({
		load(i)
		profile[i,"nsnps"]<-nrow(snps)
		profile[i,"hsq"]<-hsq[1]
		profile[i,"hsq.se"]<-hsq[2]
		profile[i,"hsq.pv"]<-hsq.pv
		profile[i, "top1.r2"]<-cv.performance["rsq", "top1"]
		profile[i, "enet.r2"]<-cv.performance["rsq", "enet"]
		profile[i, "lasso.r2"]<-cv.performance["rsq", "lasso"]
		profile[i, "lasso.pv"]<-cv.performance["pval", "lasso"]
		profile[i, "enet.pv"]<-cv.performance["pval", "enet"]
		profile[i, "top1.pv"]<-cv.performance["pval", "top1"]
	}, error=function(e){})
}

#write.table(profile, file=snakemake@output[1], quote=FALSE, col.names=FALSE)
write.table(profile, file=snakemake@output[[1]], quote=FALSE, col.names=TRUE, row.names=FALSE)