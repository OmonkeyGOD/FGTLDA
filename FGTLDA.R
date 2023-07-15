#####Pairwise FGT & LET 
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Must provide 2 arguments in the following order: input file and output file", call.=FALSE)
}
Pvalue <- c()
k <- 1
x_n <- read.table(args[1], sep = "\t",
                  header=TRUE, row.names = 1) #SNP results in binary format
SNP_count=ncol(x_n)
x_n_t <- transpose(x_n)
rownames(x_n_t) <- colnames(x_n)
colnames(x_n_t) <- rownames(x_n)
x_n = x_n_t
x_n <- x_n[order(rownames(x_n)),]
nucle_eff_frame <- data.frame(matrix("", ncol = 0, nrow = SNP_count)) #row number set as the sample number
for (i in 1 : length(colnames(x_n))){
  POS=colnames(x_n[i])
  nucle_eff_frame[,POS]<-x_n[,i]
}
sink(args[2]) #output file
for(p in 1:(length(colnames(nucle_eff_frame))-1)){  # pairwise comparison of SNP loci
  for (q in (p+1):length(colnames(nucle_eff_frame))){
    p_na=nucle_eff_frame[,p]
    q_na=nucle_eff_frame[,q]
    p_x=which(p_na=="X") 
    p_n=which(p_na=="N")
    q_x=which(q_na=="X")
    q_n=which(q_na=="N")
    pq_na=unique(c(p_x,p_n,q_x,q_n)) # samples with missing data at any of the two SNP loci
    if(length(pq_na)>SNP_count-2){ #skip SNP loci missing data in over 535 samples
      next
    }
    if(length(pq_na)!=0){
      p_na_r=as.numeric(as.character(unlist(p_na[-pq_na])))
      q_na_r=as.numeric(as.character(unlist(q_na[-pq_na])))
    }else{
      p_na_r=as.numeric(as.character(unlist(p_na)))
      q_na_r=as.numeric(as.character(unlist(q_na)))
    }
    len=length(p_na_r) # sample number after removing samples with missing data at any of the two loci
    p_na_r_sum=sum(p_na_r) #number of alternative allele at site 1
    q_na_r_sum=sum(q_na_r) #number of alternative allele at site 2
    v1=c(p_na_r_sum,len-p_na_r_sum)#v1=site1(#1,#0)
    p1=p_na_r_sum/len # frequency of alternative allele at site 1
    p2=(len-p_na_r_sum)/len #frequency of reference allele at site 1
    v2=c(q_na_r_sum,len-q_na_r_sum)#v1=site2(#1,#0)
    q1=q_na_r_sum/len # frequency of alternative allele at site 2
    q2=(len-q_na_r_sum)/len #frequency of reference allele at site 2
    m=rbind(v1,v2) 
    flag0=0 #four combinations of two sites
    flag2=0
    flag10=0
    flag01=0
    for(i in 1:len){
      s=p_na_r[i]+q_na_r[i]
      if(s==0){
        flag0=flag0+1
      }else if (s==2){
        flag2=flag2+1
      }else if (p_na_r[i]>q_na_r[i]){
        flag10=flag10+1
      }else if (p_na_r[i]<q_na_r[i]){
        flag01=flag01+1
      }
    }
    if(flag0>0 & flag2>0 & flag01>0 & flag10>0 ){ # SNP pairs past the FGT 
      cat("OO:")
      cat(flag0)
      p22=flag0/len
      cat("\t")
      cat("II:")
      cat(flag2)
      p11=flag2/len
      cat("\t")
      cat("IO:")
      cat(flag10)
      p12=flag10/len
      cat("\t")
      cat("OI:")
      cat(flag01)
      p21=flag01/len
      cat("\t")
      cat(colnames(nucle_eff_frame[p]))
      cat("\t")
      cat(colnames(nucle_eff_frame[q]))
      cat("\t")
      D = p11*p22-p12*p21
      if(D>0){
        Dmax=min(p1*q2,p2*q1)
      }else{
        Dmax=max(-p1*q1,-p2*q2)
      }
      D1=D/Dmax #normalized D
      r=D/((p1*p2*q1*q2)^(1/2)) # r
      pvalue=pchisq(r*r*len, df=1, lower.tail=FALSE) 
      Pvalue[k] <- pvalue
      k <- k + 1
      cat("D:")
      cat("\t")
      cat(D1)
      cat("\t")
      cat("r2:")
      cat("\t")
      cat(r*r)
      cat("\t")
      cat("p.value:")
      cat("\t")
      cat(pvalue)
      cat("\n")
    }
  }
}
FDR=p.adjust(Pvalue,"fdr")
sink()
lapply(FDR, write, paste(args[2], ".FDR"), append=TRUE, ncolumns=1000)




