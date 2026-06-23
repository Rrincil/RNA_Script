library(gmodels)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ggrepel)
library(data.table)

args <- commandArgs()

# expr <- read.table(args[6],header=T,row.names=1,sep="\t",check.names=F)
expr <- fread(args[6])
expr <- expr[!duplicated(expr[,1]),]

# expr = expr[,c(-2,-13,-7,-11)]
head(expr)
expr <- as.data.frame(expr)

if(length(args)==8){
colnum1 <- grep(paste("^",args[7],"-[0-9]+",sep=""),colnames(expr))
colnum2 <- grep(paste("^",args[8],"-[0-9]+",sep=""),colnames(expr))
expr <- expr[,c(colnum1,colnum2)]
}


group <- rep(0,ncol(expr))
sorted <- rep(0,(length(args)-6))
for(m in 7:length(args)) {
M <- m-6
sorted[M] <- args[m]
colnum <- grep(args[m],colnames(expr))
for(n in 1:length(colnum)) {
group[colnum[n]] <- args[m]
}
}

keep <- apply(expr,1,mean)>0.5
expr <- expr[keep,]
# data <- t(as.matrix(log2(expr+1)))

#expr <- t(scale(t(expr)))
data <- t(as.matrix(expr))

###########################################################

data.pca <- fast.prcomp(data,retx=T,scale=F,center=T)

a <- summary(data.pca)
tmp <- a[4]$importance
pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100
pro2 <- as.numeric(sprintf("%.3f",tmp[2,2]))*100

pc <- as.data.frame(a$x)
pc$group <- group
pc$names <- rownames(pc)

xlab <- paste("PC1(",pro1,"%)",sep="")
ylab <- paste("PC2(",pro2,"%)",sep="")
pca <- ggplot(pc,aes(PC1,PC2)) +
geom_hline(yintercept=0,linetype="dashed",color="grey") +
geom_vline(xintercept=0,linetype="dashed",color="grey") +
stat_ellipse(geom="polygon",alpha=1/10,
aes(fill=group,colour=group)) +
geom_point(size=2,pch=21,colour="black",aes(fill=group)) +
labs(x=xlab,y=ylab) +
theme_bw() +
theme(axis.text=element_text(size=6,color="black"),
axis.title=element_text(size=6,color="black"),
plot.title=element_text(size=6,color="black"),
legend.title=element_text(size=6,color="black"),
legend.text=element_text(size=6,color="black"),
legend.background=element_rect(fill="white",colour="black")) +
labs(color="Sample") +
geom_rug(aes(colour=group)) +
scale_color_aaas(limits=sorted) +
scale_fill_aaas(limits=sorted) +
#scale_x_continuous(expand=c(0,0),limits=c(-40,40)) +
#scale_y_continuous(expand=c(0,0),limits=c(-40,40))
geom_text_repel(aes(label=names),size=6/3,force=10)
out1 <- paste(sorted,collapse="-vs-")
out2 <- paste("PCA.",out1,sep="")
ggsave(paste(out2,".pdf",sep=""),pca,width=6,height=5)

###########################################################

library(Rtsne)

data_unique <- unique(data)
set.seed(12345)

tsne_out <- Rtsne(as.matrix(data_unique),perplexity=1)
#tsne_out <- Rtsne(as.matrix(data_unique),perplexity=length(sorted))

dm <- as.data.frame(tsne_out$Y)
dm$group <- group
dm$names <- rownames(data_unique)
colnames(dm)[1:2] <- c("DM1","DM2")

tsne <- ggplot(dm,aes(DM1,DM2)) +
geom_hline(yintercept=0,linetype="dashed",color="grey") +
geom_vline(xintercept=0,linetype="dashed",color="grey") +
geom_point(size=2,aes(colour=group)) +
labs(x=c("Dimension 1"),y=c("Dimension 2")) +
theme_bw() +
theme(axis.text=element_text(size=6,color="black"),
axis.title=element_text(size=6,color="black"),
plot.title=element_text(size=6,color="black"),
legend.title=element_text(size=6,color="black"),
legend.text=element_text(size=6,color="black"),
legend.background=element_rect(fill="white",colour="black")) +
labs(color="Sample") +
stat_ellipse(geom="polygon",alpha=1/10,linetype="dashed",aes(fill=group,colour=group)) +
geom_rug(aes(color=group)) +
scale_color_nejm(limits=sorted,guide=FALSE) +
scale_fill_nejm(limits=sorted,guide=FALSE) +
geom_text_repel(aes(label=names),size=6/3,force=10)
out3 <- paste("t-sne.",out1,sep="")
ggsave(paste(out3,".pdf",sep=""),tsne,width=2.5,height=2.5)

###########################################################

library(dendextend) ## using chaining operator
library("ape")

method <- c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
for(i in 1:length(method)){
hc <- hclust(dist(data),method=method[i])

labels <- names(cutree(hc,1))
color <- rep(0,length(labels))
for(m in 7:length(args)) {
M <- m-6
colnum <- grep(args[m],labels)
for(n in 1:length(colnum)) {
color[colnum[n]] <- rainbow(length(args)-6+1)[M]
}
}

out4 <- paste("hclust_",method[i],".",out1,sep="")
#pdf(file=paste(out4,".circle.pdf",sep=""),width=2.5,height=2.5)
#par(mar=c(0,0,0,0),oma=c(1,1,1,1))
#plot(as.phylo(hc),type="fan",tip.color=color,cex.axis=0.5,cex.lab=0.5,cex=0.5)
#dev.off()

color <- rep(0,length(labels(hc)))
for(m in 7:length(args)) {
M <- m-6
colnum <- grep(args[m],labels(hc))
for(n in 1:length(colnum)) {
color[colnum[n]] <- pal_aaas("default")(length(args)-6+1)[M]
#color[colnum[n]] <- rainbow(length(args)-6+1)[M]
}
}

dend <- as.dendrogram(hc) %>% 
set("labels_col",color) %>% 
# set("leaves_pch",17) %>% 
set("leaves_pch",17) %>% 
set("labels_cex",0.5) %>% 
set("leaves_cex",0.5) %>%  #下面三角形大小
set("leaves_col",color)
out4 <- paste("hclust_",method[i],".",out1,sep="")
rw <- length(labels(hc))*30
pdf(file=paste(out4,".pdf",sep=""),width=rw/100,height=2.5)
# par(mar=c(3,3.5,0,0),oma=c(0,1,0,0))
par(mar=c(2.5,2.5,0.5,0),oma=c(0,1,0,0)) #下，左、上、右

library(ggdendro)

# ggdendrogram(hc)
plot(dend,ylab="Height",cex.axis=0.5,cex.lab=0.5)
dev.off()
}

#library(mixOmics)
#plsda <- plsda(data,group,ncomp=length(group))
##set ncomp for performance assessment later
#out5 <- paste("pls-da.",out1,sep="")
#pdf(file=paste(out5,".pdf",sep=""),
#width=5.5,height=5)
#plotIndiv(plsda,comp=1:2,group=group,
#ind.names=FALSE,ellipse=TRUE,legend=TRUE,title='PLSDA')
#dev.off()
