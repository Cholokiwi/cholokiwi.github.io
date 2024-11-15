# NOTE: this data is only for illustration purposes
# It does not reflect real data and genotypes have been manually edited (i.e. 'invented') to help tease out concepts we wanted to show in the chapter 
# Do not use this to make any inferences about 'reality'   

# read in data
samples=read.table("3breeds.txt",header=F,sep=" ")
dim(samples)

# calculate frequencies
M=matrix(NA,nrow(samples),3)
colnames(M)=c("Hanwoo","Angus","Brahman")

# hanwoo
M[,1]=apply (samples[,1:25],1,function(x) sum(x)/(length(x)*2))

# angus
M[,2]=apply (samples[,26:50],1,function(x) sum(x)/(length(x)*2))

# brahman
M[,3]=apply (samples[,51:75],1,function(x) sum(x)/(length(x)*2))

# FST
meansB=rowMeans(M) # average allele frequency across populations
alleleVar=meansB*(1-meansB) # p*q variance
meanDevB=M-meansB # deviation of each population from mean
FST=meanDevB^2/alleleVar # deviation squared divided by var

library(lokern)
smoothHanwoo= lokerns(FST[,1], n.out=77)
smoothAngus = lokerns(FST[,2], n.out=77)
smoothBrahman = lokerns(FST[,3], n.out=77)

plot(FST[,1],type="l",xlab="SNP",ylab="Fst",col="gray")
lines(smoothHanwoo$est,type="l",col="red",lwd=6)

write.table(FST,"FST3breeds.txt",quote=F,row.names=F)
smoothed=data.frame(SNPindex=smoothHanwoo$x.out,Hanwoo=smoothHanwoo$est,Angus=smoothAngus$est,Brahman=smoothBrahman$est)
write.table(smoothed,"FSTsmooth3breeds.txt",quote=F,row.names=F)

# GRM
M=as.matrix(samples)
p=apply (M,1,function(x) sum(x)/(length(x)*2)) # frequency of second (minor) allele = num alleles / total num alleles

M=M-1
P=2*(p-0.5) # deviation from 0.5 - P should use base population frequencies!
Z=M-P
ZtZ = t(Z) %*% Z
d=2*sum(p*(1-p))
G=ZtZ/d

cols=c(rep("red",25),rep("blue",25),rep("green",25))
# heatmap
heatmap(G,symm=T,col=gray.colors(16,start=0,end=1),RowSideColors=cols,ColSideColors=cols)
legend("topright",c("Hanwoo","Angus","Brahman"),fil=c("red","blue","green"),cex=1)

# singular value decomposition of GRM
SVD=svd(G)
plot(SVD$v[,2],-1*SVD$v[,1],cex.main=0.9,main="singular value decomposition",xlab="PC1",ylab="PC2",col=cols,pch=20)
legend("topright",c("Hanwoo","Angus","Brahman"),fil=c("red","blue","green"),cex=1)

