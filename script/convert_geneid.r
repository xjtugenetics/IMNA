library(org.Hs.eg.db)
library(clusterProfiler)

argv <- commandArgs(TRUE)
data <- read.table(argv[1],header=T)  ##input
oup <- argv[2]  ##output
corvert_col <- as.numeric(argv[3]) # target col
geneid <- argv[4]   ##SYMBOL or ENTREZID


dim(data)
colnames(data)[corvert_col] <- geneid

if(geneid=="SYMBOL"){

data2 = bitr(data[,corvert_col], fromType="SYMBOL", toType="ENTREZID",OrgDb="org.Hs.eg.db")
data_cov <- merge(data2,data,by='SYMBOL')
data_cov <- data_cov[,-1]
head(data_cov)

}else if(geneid=="ENTREZID"){


data2 = bitr(data[,corvert_col], fromType="ENTREZID", toType="SYMBOL",OrgDb="org.Hs.eg.db")
data_cov <- merge(data2,data,by='ENTREZID')
data_cov <- data_cov[,-1]
head(data_cov)

}

write.table(data_cov,file=oup,quote=F,col=T,row=F,sep="\t")

