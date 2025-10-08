# Function definitions at end of script

clim = read.table("summerClimateMeans.txt",header=T)
row.names(clim) = clim$PopFull
names(clim)[4:8] = c("tmin","tmax","tdmean","vpdmin","vpdmax")
pop_ord = c("DPcontemporary","DPhistorical","LAhistorical","SMcontemporary","SMhistorical","SSMcontemporary","TJhistorical")

#---- Empirical changes
quasi_with_popAFs = read.table("MAF_master_primaryCandidates_wLCWP2.5k.txt",header=T,stringsAsFactors = F)
quasi_with_popAFs$lcwp_predicted_glm = glm_predict_climate(quasi_with_popAFs,pop_ord,clim,clim["LCWPcontemporary",])
quasi_with_popAFs$Rs = get_Pearson_corr(clim,pop_ord,quasi_with_popAFs,which(names(quasi_with_popAFs)=="Var"))
# Compare to a predicted values
quasi_with_popAFs$concordant2017 = with(quasi_with_popAFs,abs(start_2017-lcwp_predicted_glm)>abs(end_2017-lcwp_predicted_glm) | p2017>0.05)
quasi_with_popAFs$concordant2018 = with(quasi_with_popAFs,abs(start_2018-lcwp_predicted_glm)>abs(end_2018-lcwp_predicted_glm) | p2018>0.05)
quasi_with_popAFs$concordant2020 = with(quasi_with_popAFs,abs(start_2020-lcwp_predicted_glm)>abs(end_2020-lcwp_predicted_glm) | p2020>0.05)
quasi_with_popAFs$sigP = apply(quasi_with_popAFs[,c("p2017","p2018","p2020")], 1, min) < 0.05
quasi_with_popAFs$concordantAll = with(quasi_with_popAFs,concordant2017==TRUE & concordant2018==TRUE & concordant2020==TRUE)
quasi_with_popAFs$concordantAllwithP = with(quasi_with_popAFs,(concordant2017==TRUE | p2017>0.05) & (concordant2018==TRUE | p2018>0.05) & (concordant2020==TRUE | p2020>0.05) & sigP==TRUE)
# Calculate the gene-wise proportion of concordant changes
gene_concordants = aggregate(concordantAllwithP~Var,data=quasi_with_popAFs,FUN=sum)
gene_all = aggregate(concordantAllwithP~Var,data=quasi_with_popAFs,FUN=length)
names(gene_concordants)[2] = "concordant"
names(gene_all)[2] = "all"
gene_counts = merge(gene_concordants,gene_all,by=c("Var"))
gene_counts$prop = gene_counts$concordant/gene_counts$all
gene_counts$Binom = NA
for(i in 1:nrow(gene_counts)){
  res = binom.test(gene_counts[i,"concordant"],gene_counts[i,"all"],p=0.12,alternative = "greater")
  gene_counts[i,"Binom"] = res$p.value
}
print(gene_counts)

#---- Null
lcwp_non_assocs = read.table("LCWPnonAssociatedSitesListV6.txt",header=T)
null_dist = generate_null(lcwp_non_assocs,100,FALSE)
get_p(na.omit(null_dist$concordantWithP),0.30)



#---- Function definitions
glm_predict_climate = function(majors, pop_ord, climate, second_climate){
  predictions = NULL
  for (i in 1:nrow(majors)){
    clim_var = as.character(majors$Var[i])
    clim_values = climate[pop_ord,clim_var]
    afs = majors[i,pop_ord]
    mod = glm(as.numeric(afs)~clim_values,family="binomial")
    predictions[i] = predict(mod,newdata=data.frame(clim_values = second_climate[1,clim_var]))
  }
  return(predictions)
}

generate_null = function(lcwp_non_assocs,nloci,win=FALSE) {
  concordant_counts = as.data.frame(matrix(NA,nrow=1000,ncol=5))
  for(j in 1:1000){
    print(paste0(j))
    all_assoc_windows = lcwp_non_assocs[sample(nrow(lcwp_non_assocs),nloci),]
    if(win==TRUE){
      all_assoc_windows$Start = all_assoc_windows$pos - 2500
      all_assoc_windows$End = all_assoc_windows$pos + 2500
      merged_mafs_wVar = c()
      for(i in 1:nrow(all_assoc_windows)){
       sub_df = lcwp_non_assocs[lcwp_non_assocs$marker == all_assoc_windows[i,1] & lcwp_non_assocs$pos > all_assoc_windows[i,"Start"] & lcwp_non_assocs$pos < all_assoc_windows[i,"End"],]
       if(dim(sub_df)[1]>0){
         sub_df$Var = all_assoc_windows[i,"Var"]
         #sub_df$Gene = all_assoc_windows[i,"Gene"]
         sub_df$Distance = abs(sub_df$pos - all_assoc_windows[i,"pos"])
         merged_mafs_wVar = rbind(merged_mafs_wVar,sub_df)
       }
      }
      df_transformed = unique(merged_mafs_wVar)
    } else {
      df_transformed = unique(all_assoc_windows)
    }
    df_transformed$Rs = runif(nrow(df_transformed),-1,1)
    df_transformed$Predicted = runif(1,0,1)
    df_transformed$SigP = apply(df_transformed[,c("p2017","p2018","p2020")], 1, min) < 0.05
    df_transformed$concordant2017 = with(df_transformed,abs(start_2017-Predicted)>abs(end_2017-Predicted))
    df_transformed$concordant2018 = with(df_transformed,abs(start_2018-Predicted)>abs(end_2018-Predicted))
    df_transformed$concordant2020 = with(df_transformed,abs(start_2020-Predicted)>abs(end_2020-Predicted))
    df_transformed$concordantAll = with(df_transformed,concordant2017==TRUE & concordant2018==TRUE & concordant2020==TRUE)
    df_transformed$concordantAllwithP = with(df_transformed,(concordant2017==TRUE | p2017>0.05) & (concordant2018==TRUE | p2018>0.05) & (concordant2020==TRUE | p2020>0.05) & SigP==TRUE)
    df_transformed$Rconcordant2017 = with(df_transformed,(slope2017<0)==(Rs<0))
    df_transformed$Rconcordant2018 = with(df_transformed,(slope2018<0)==(Rs<0))
    df_transformed$Rconcordant2020 = with(df_transformed,(slope2020<0)==(Rs<0))
    df_transformed$RconcordantAll = with(df_transformed,Rconcordant2017==TRUE & Rconcordant2018==TRUE & Rconcordant2020==TRUE)
    df_transformed$RconcordantAllwithP = with(df_transformed,(Rconcordant2017==TRUE | p2017>0.05) & (Rconcordant2018==TRUE | p2018>0.05) & (Rconcordant2020==TRUE | p2020>0.05) & SigP==TRUE)
    concordant_counts[j,1] = sum(df_transformed$concordantAll)/nrow(df_transformed)
    concordant_counts[j,2] = sum(df_transformed$concordantAllwithP)/nrow(df_transformed)
    concordant_counts[j,3] = sum(df_transformed$RconcordantAll)/nrow(df_transformed)
    concordant_counts[j,4] = sum(df_transformed$RconcordantAllwithP)/nrow(df_transformed)
    concordant_counts[j,5] = nrow(df_transformed)
  }
  names(concordant_counts) = c("concordantAll","concordantWithP","RconcordantAll","RconcordantWithP","N")
  return(concordant_counts)
}

get_p = function(dist,x){
  p=1-(sum(dist<x)/length(dist))
  return(p)
}
