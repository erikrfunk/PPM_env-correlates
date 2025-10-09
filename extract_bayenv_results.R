library(dplyr)
results = read.table("AllResults_withCoords.txt",header=F)
idxs = c(1,2,3,4,5)
vars = c("tmin","tmax","tdmean","vpdmin","vpdmax")

for(i in 1:length(vars)){ #Extract the columns corresponding to the desired variables
  cuts = c(idxs[i]*3+1,idxs[i]*3+2)
  bfs = results[,c(1,2,cuts)]
  bfs$clim = vars[i]
  names(bfs) = c("Chr","Pos","BF","Rho","Clim")
  write.table(bfs,paste0(vars[i],"results.txt"),row.names = F,quote=F)
}

f_results = list.files(".","*results.txt$",full.names = F)
BFpercent = 0.01
Rhopercent = 0.05
for(i in f_results){
  results = read.table(i,header=T)
  variable = gsub("results.txt","",i)
  topBF = results[order(results[,3],decreasing = T),][1:round(nrow(results)*BFpercent),c(1,2,3,5)]
  topRho = results[order(abs(results[,4]),decreasing = T),][1:round(nrow(results)*Rhopercent),c(1,2,4,5)]
  topsBoth = inner_join(topBF,topRho,by=c("Chr","Pos","Clim"))
  write.table(topsBoth[,c("Chr","Pos","Clim","BF","Rho")],file = paste0(variable,"_top1p5p.txt"),quote=F,row.names=F)
}

