library(dplyr)
setwd("/home/centos/USS/erik/PPM/env_corrs/bayenv2/")

# ALLELES ----
mafs = read.table("/home/centos/USS/erik/PPM/env_corrs/maf_changes/revised_mafs/newSnpCalls/AllPopsDF_joinedPPM_realigned_newSnpCalls_MAFs0.2.txt",header=T)
names(mafs) = gsub("_MAF","",names(mafs))
popmap = read.table("/home/centos/USS/erik/PPM/joined_historical_and_contemporary/realigned/AllPPM_realigned_samplemap.txt",header = F,stringsAsFactors = T)
popmap$uniqueID = paste0(popmap$V1,popmap$V2)
#pop_order = gsub("-","",as.character(levels(as.factor(popmap$V1))))
pop_order = unique(popmap$uniqueID)
pop_size = as.numeric(dplyr::count(popmap, uniqueID)[,2])

# Optionally prune the SNP set based on an arbitrary distance to make the panel more wieldy 
chrs = unique(mafs$chromo)
fout = "mafs_1kb_thinned.txt"
write.table(paste(names(mafs),collapse='\t'),fout,quote=F,row.names = F,col.names=F,append = T)
for (i in chrs){
  df = mafs[which(mafs$chromo == i),]
  pos = 0
  for (j in 1:nrow(df)){
    if (df[j,"position"] > pos){
      write.table(df[j,],fout,quote=F,row.names = F,col.names=F,append = T)
      pos = df[j,"position"] + 1000
    }
  }
}

mafs_coords = mafs[,c(1,2)]
mafs_counts = round(mapply('*',mafs[,pop_order],pop_size))
per_site_counts = apply(mafs_counts,1,"sum")
nrow(mafs_counts)==length(per_site_counts)
mafs_counts = mafs_counts[per_site_counts > 2 & per_site_counts < (sum(pop_size)-2),] # Require at least three gene copies of each allele
mafs_coords = mafs_coords[per_site_counts > 2 & per_site_counts < (sum(pop_size)-2),]
major_counts = abs(sweep(mafs_counts,2,pop_size))
mafs_counts = cbind(mafs_coords,mafs_counts)
major_counts = cbind(mafs_coords,major_counts)
row.names(mafs_counts) = seq(1:nrow(mafs_counts)) * 2 - 1
row.names(major_counts) = seq(1:nrow(major_counts)) * 2

all_counts = as.data.frame(rbind(mafs_counts,major_counts))
all_counts = all_counts[order(as.numeric(row.names(all_counts))), ]
write.table(all_counts, "ppm_counts_with_coords.txt",quote=F,row.names = F,col.names = F,sep='\t') # If subsequently splitting this file, need to keep coords with counts by writing all columns, not just 3-
#write.table(mafs_coords,"maf_coords",row.names=F,col.names=F,quote=F)

#Subsample for covariance matrix estimation
nsnps = round(runif(10000,1,nrow(mafs_coords)))
coord_sub = mafs_coords[nsnps,]
min_sub = mafs_counts[nsnps,]
maj_sub = major_counts[nsnps,]
all_counts_sub = as.data.frame(rbind(min_sub,maj_sub))
all_counts_sub = all_counts_sub[order(as.numeric(row.names(all_counts_sub))), ]
write.table(coord_sub,"Subset_maf_coordinates_for_covmat.txt",row.names = F,quote=F,col.names=F,sep="\t")
write.table(all_counts_sub[,3:ncol(all_counts_sub)],"Subset_maf_counts_for_covmat.txt",quote = F, row.names=F,col.names=F,sep = "\t")

# Or use the coords to break up the file in to Chrs
for(i in 1:30){ # Set equal to the number of chrs wanting to include
  out_dir = paste0("CHR",i)
  scaff = paste0("HiC_scaffold_",i)
  write.table(all_counts[all_counts==scaff,3:ncol(all_counts)],file=paste0(out_dir,"/",out_dir,"allele_counts"),quote=F,row.names=F,col.names = F,sep="\t")
  write.table(mafs_coords[which(mafs_coords$chromo==scaff),],file = paste0(out_dir,"/snp_coords"),quote=F,row.names=F,col.names=F)
}
# Or evenly
nfiles = 10



# CLIMATE ----
clim_means = read.table("/home/centos/USS/erik/PPM/env_corrs/summerClimateMeans.txt",header=F,stringsAsFactors = F)
row.names(clim_means) = clim_means[,1]
clim_trimmed = clim_means[pop_order,c(4:ncol(clim_means))]
clim_trimmed = t(clim_trimmed)
for (i in 1:nrow(clim_trimmed)){ # Standardize each climate variable (x - mean / SD)
  clim_trimmed[i,] = (clim_trimmed[i,] - mean(clim_trimmed[i,])) / sd(clim_trimmed[i,])
}
write.table(clim_trimmed,"Climate_standardized_means.txt",row.names = F,col.names=F,quote=F,sep='\t')
