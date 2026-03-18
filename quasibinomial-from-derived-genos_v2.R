# Get derived gprobs from extract-individual-deriveds.sh
start_genos = read.table("start_2017_derived.txt",header=T)
end_genos = read.table("end_2017_derived.txt",header=T)
nrow(start_genos)==nrow(end_genos)

reg_results = data.frame("slope"=rep(NA,nrow(start_genos)),
                         "p"=rep(NA,nrow(start_genos)))
for(i in 1:nrow(start_genos)){
    cat(paste0(i,"\n"))
    locus = data.frame("Time"=c(rep("start",ncol(start_genos[5:ncol(start_genos)])),
                                rep("end",ncol(end_genos[5:ncol(end_genos)]))),
                       "Derived"=c(as.numeric(start_genos[i,5:ncol(start_genos)]),
                                   as.numeric(end_genos[i,5:ncol(end_genos)])))
    locus$Time = factor(locus$Time,levels=c("start","end"))
    glm_results = glm(Derived~Time,data=locus,family=quasibinomial)
    reg_results[i,1] = glm_results$coefficients[2]
    reg_results[i,2] = summary(glm_results)$coefficients[2,4]
}

reg_results = cbind(start_genos[,1:4],reg_results)
write.table(reg_results,file="2017.txt",row.names = F,quote = F)