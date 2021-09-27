library(magrittr)
#step1
data = read.delim("example_paml_output.txt", stringsAsFactors = F, header = F)
seq = grep("\\.{3,3}",data$V1,value = T)
seq1 = strsplit(seq,split = "\\s+") %>% sapply(function(x) x[2]) %>% sapply(function(x) substr(x, 2, (nchar(x)-1)))
seq2 = strsplit(seq,split = "\\s+") %>% sapply(function(x) x[5]) %>% sapply(function(x) substr(x, 2, (nchar(x)-1)))
lnL = grep("lnL",data$V1,value = T)
lnL= strsplit(lnL,split = "\\s+") %>% sapply(function(x) x[3])
readns = function(x){
  t = as.numeric(x[2])
  S = as.numeric(x[4])
  N = as.numeric(x[6])
  dNdS = as.numeric(x[8])
  dN = as.numeric(x[11])
  dS = as.numeric(x[14])
  return(c(t,S,N,dNdS,dN,dS))
}
ns = grep("dN/dS=",data$V1,value = T)
ns = t(strsplit(ns,split = "\\s+") %>% sapply(readns))
all = as.data.frame(cbind(seq1, seq2, lnL, ns))
colnames(all) = c("seq1","seq2","lnL","t","S","N","dNdS","dN","dS")
rownames(all) = c(1:nrow(all))

#step2
library(data.table)
groupfile = read.csv("example_groupfile.csv")
grouplist = as.character(as.data.frame(table(groupfile$group))$Var1)
group = data.table(groupfile[, c("genebank_id", "group")], key = "genebank_id")
all$group1 = group[all$seq1]$group
all$group2 = group[all$seq2]$group

all$dNdS[(all$dN == 0 & all$dS != 0)] = 0
pairs = as.data.frame(NULL)
z = length(grouplist)
for(x in 1:z){
  group_num1 = grouplist[x]
  for(y in 1:z){
    group_num2 = grouplist[y]
    a_pair = all[which((all$group1 == group_num1 & all$group2 == group_num2) | (all$group1 == group_num2 & all$group2 == group_num1)),]
    a_pair$N = as.numeric(a_pair$N)
    a_pair$dN = as.numeric(a_pair$dN)
    a_pair$S = as.numeric(a_pair$S)
    a_pair$dS = as.numeric(a_pair$dS)
    N = mean((a_pair$N*a_pair$dN), na.rm = T)
    S = mean((a_pair$S*a_pair$dS), na.rm = T)
    a_pair = cbind(group_num1,group_num2,N,S)
    pairs = rbind(pairs,a_pair)
  }
}

pairs = pairs[-which(pairs$group_num1 == pairs$group_num2),]
groupfile$group = as.character(groupfile$group)
cluster = unique(data.table(groupfile[, c("group", "cluster")], key = "group"))
pairs$cluster1 = cluster[pairs$group_num1]$cluster
pairs$cluster2 = cluster[pairs$group_num2]$cluster
between = pairs[which(pairs$cluster1 != pairs$cluster2),]
within = pairs[which(pairs$cluster1 == pairs$cluster2),]

between_N = sum(as.numeric(between$N),na.rm = T)/2
between_S = sum(as.numeric(between$S),na.rm = T)/2
within_N = sum(as.numeric(within$N),na.rm = T)/2
within_S = sum(as.numeric(within$S),na.rm = T)/2
pMK = matrix(c(round(between_N),round(between_S),round(within_N),round(within_S)),ncol = 2,nrow = 2)
colnames(pMK) = c("beween","within")
rownames(pMK) = c("N","S")
write.csv(pMK,"pMK_output.csv")
