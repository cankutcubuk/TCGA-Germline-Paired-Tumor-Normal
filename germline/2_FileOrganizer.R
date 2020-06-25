rm(list=ls())

# This script organizes the .bam files (downladed from api.gdc.cancer.gov) and distrubutes them into the folders of different cancer types.
# Usage: Rscript 2_FileOrganizer.R --workdir "/mnt/scratch/DGE/MOPOPGEN/ccubuk/gdc/bam_splice/allbams/" --annot /mnt/scratch/DGE/MOPOPGEN/ccubuk/gdc/files/SampleAnnots.tsv

library(R.utils,quietly=T,verbose=F,warn.conflicts=F)
args <- commandArgs(trailingOnly = F, asValues = T,excludeEnvVars = F)
workdir <- args[["wdir"]] # the folder where the all .bam files downloaded. 
annot <- args[["annot"]] # sample annotation file
setwd(workdir)

downs <- list.files(workdir, pattern = "bam")
downs <- gsub("[.]bam","",downs)
metadata <- read.delim(annot, row.names = NULL, stringsAsFactors = F, quote = "")
metadata <- metadata[which(metadata$id %in% downs),]

cat("Discarded from the analysis. Unclear eviddence for being tumor or normal. \nSample type with sample numbers:")
print(table(metadata$sample_type[which(is.na(metadata$simpletype))]))
metadata <- metadata[which(!is.na(metadata$simpletype)),]

for(p in unique(metadata$project_id)){

  tfile <- metadata$id[metadata$project_id==p & metadata$simpletype=="Normal"]
  nfile <- metadata$id[metadata$project_id==p & metadata$simpletype=="Tumor"]
  
  nt <- length(tfile)
  nn <- length(nfile)
  
  if(nt==0 | nn==0){ warning(paste0(p, " is discarded from the analysis. It has 0 ", ifelse(nn==0, "normal samples", "tumor samples"))); next}
  if(nt/(nt+nn)<0.1 | nt/(nt+nn)>0.9){warning(paste0(p, " is discarded from the analysis. It has unbalanced tumor/normal sample size"))  }
  print(paste0(length(tfile), "-",length(nfile)))
  
   pcreate <- gsub("-","_",gsub("[.]","",p))
   nfolder <- normalizePath(paste0(workdir,"/../tidy/",pcreate,"/Normals/"))
   tfolder <- normalizePath(paste0(workdir,"/../tidy/",pcreate,"/Tumors/"))
   
   dir.create(nfolder,recursive = T,showWarnings = F)
   dir.create(tfolder,recursive = T,showWarnings = F)
   
   system(paste0("cp ", paste0(paste0(workdir,"/",nfile,".bam"), collapse = " "), " ", nfolder)) 
   system(paste0("cp ", paste0(paste0(workdir,"/",tfile,".bam"), collapse = " "), " ", tfolder)) 

   write.table(paste0(nfile,".bam"), file = paste0(normalizePath(paste0(workdir,"/../tidy/",pcreate)),"/allNormals.txt"), row.names = F,quote = F, col.names = F)
   write.table(paste0(tfile,".bam"), file = paste0(normalizePath(paste0(workdir,"/../tidy/",pcreate)),"/allTumors.txt"), row.names = F,quote = F, col.names = F)
   
   print(paste0(p, "...DONE"))
}
message(paste0("Your files ordered under ", normalizePath(paste0(workdir,"/../")), "/tidy/{project_id}/"))
