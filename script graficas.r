WD <- "/home/user/Escritorio/Biologia/TESIS/Datos_Maiz/DATA/"
setwd("/home/user/Escritorio/MAIZ/Articulo_tesis/Final/Figures")

## Load required libraries
library (adegenet)
library(poppr)
library(ape)
library(NMF)
library(snpStats)
library(ggplot2)
library(gdsfmt)
library(SNPRelate)
library(grid)
library(gridExtra)
require(maptools)
require(ggmap)

##################################################################################################
## Load SNP data
bed<-paste(WD,"/maiz_listo.bed",sep="")
bim<-paste(WD,"/maiz_listo.bim",sep="")
fam<-paste(WD,"/maiz_listo.fam",sep="")
ped<-paste(WD,"/maiz_listo.ped",sep="")
raw<-paste(WD,"/maiz_listo.raw",sep="")
map <-paste(WD, "/maiz_listo.map",sep="")

## Load sample data
subject.support<-read.csv(paste(WD,"/subject.support.csv",sep=""),header=TRUE) 
# Load information about domestication and improvement SNPs
snpsextract<-read.csv(paste(WD,"/snp.support.csv",sep=""),header=TRUE)

## Integrate SNP data for SNPStats package handling
maiz<-read.plink(bed,bim,fam)
snps<-maiz$genotype
snp.support<-maiz$map
## Integrate SNP data for Adegenet package handling
maizgl <- read.PLINK(file=raw,mapfile=map)
## Remove duplicated SNPs
glmatrix <- as.matrix(maizgl)
glmatrix <- glmatrix[,-grep("HET",colnames(glmatrix))]
maizgl<-as.genlight(glmatrix)

## Integrate data for SNPRelate package handling
genofile <- snpgdsOpen(paste(WD,"maiz_listo.gds",sep=""))

snpset <- snpgdsLDpruning(genofile,ld.threshold = 1, autosome.only=F) ## Remove monomorphic SNPs
snpset.id <- unlist(snpset) #get all selected snp id

## Subset SNP matrix keeping only Domestication and improvement SNPs
glmatrixdom <- glmatrix[, pmatch(snpsextract[,2],colnames(glmatrix))]
glmatrixdom <- glmatrixdom[, colSums(is.na(glmatrixdom)) != nrow(glmatrixdom)]

### Fst measurement for every SNP among races, keeping SNP name and removing NAs
fraza <- Fst (snps,subject.support[,'raza'])
names(fraza$Fst)<-colnames(snps)
fraza$Fst<-fraza$Fst[!is.na(fraza$Fst)]

### Function to identify SNPs with an Fst value higher than certain percent treshold. Takes as input a named vector with the Fst value for SNPs and a givena treshold. 
fst.alta <- function (snpsfst, percent){ ## snpsfst=resultado de fst por snp para cierto agrupamiento, percent=.99 i.e 99%
  alta <- quantile(snpsfst,probs=percent) ## obtener valor de Fst para el percent deseado
  fst.snp <- data.frame (SNPs=c(1:length(snpsfst[snpsfst>=0])), Fst = snpsfst[snpsfst>=0], Dif = snpsfst[snpsfst>=0]>alta) # armar matriz con snps(Fst>=0), sus valores de Fst y su condicion(>/< percent)
  grafica <- ggplot(fst.snp, aes(x=SNPs, y=Fst, colour=Dif)) + geom_point(shape=19)+scale_colour_brewer(palette="Set1")+guides(colour=FALSE)+ylab(expression(paste("F"[ST],sep=""))) # scatterplot Fst < percent
  snps.alta <- snpsfst[snpsfst>alta] #lista de snps Fst > treshold
  return(list(Fst.value=alta, snp.info=snps.alta,fst.snp=fst.snp, grafica=grafica))
}
################################################################################################################################
## Identify high Fst SNPs among races
raza.99<-fst.alta(fraza$Fst,.99)
snps_Fstalta<-data.frame(names(raza.99$snp.info))
###########################################################################
### Subset SNP matrix keeping only landrace high Fst  SNPs
snpsfst99 <- names(raza.99$snp.info)
glmatrixfst <- glmatrix[, pmatch(snpsfst99, colnames(glmatrix))]
glmatrixfst <- glmatrixfst[, colSums(is.na(glmatrixfst)) != nrow(glmatrixfst)]
################################################################################################################################
### PCA of the three SNP sets
#### Function to create PCA using a given set SNPs and grouping of samples.
PCA <- function(genofile, grouping, snps=snpset.id){ 
  used.snps <- snps[match(snpset.id,snps)]
  used.snps <- used.snps[complete.cases(used.snps)]
  pca <- snpgdsPCA (genofile, snp.id=used.snps)
  
  pc.percent<- pca$varprop*100
  
  col<-c("red","gold","green3","dodgerblue","magenta4")
  sample.id <- read.gdsn (index.gdsn(genofile, "sample.id"))
  pop_code <- as.vector (subject.support[,grouping])        
  tab <- data.frame (sample.id=pca$sample.id,
                     Landrace=factor(pop_code)[match(pca$sample.id, sample.id)],
                     PC1 = pca$eigenvect[,1], PC2 = pca$eigenvect[,2],
                     stringsAsFactors = FALSE)
  plot <- ggplot(tab, aes(x=PC2, y=PC1, colour=Landrace)) + geom_point(shape=19)+ scale_colour_manual(values = col)+theme_bw()
  return(list(PCA=pca,PCpercent=pc.percent,table=tab,grafica=plot,snps=used.snps))
}
##############################################################
snpsdom <- as.character(snpsextract[,2])
snpsdom <- snpsdom[match(snpset.id,snpsdom)] ## select the domestication snps present in the snpset
snpsdom <- snpsdom[complete.cases(snpsdom)]

snpsfst99 <- names(raza.99$snp.info) ## names of high Fst SNPs

pca.all.raza <- PCA(genofile=genofile, snps = snpset.id, grouping = 'raza') # PCA with all SNPs
pca.dom.raza <- PCA(genofile=genofile, snps = snpsdom, grouping = 'raza')# PCA with domestication and improvement SNPs
pca.fst99.raza <- PCA(genofile=genofile, snps = snpsfst99, grouping = 'raza') # PCA with high Fst SNPs

################################################################################################################################
## Clustering analysis
gldom <- as.genlight(glmatrixdom)
glfst <- as.genlight(glmatrixfst,centers=3)

set.seed(12)
allclusters <- find.clusters.genlight(maizgl,n.pca = 48)
domclusters <- find.clusters.genlight(gldom,n.pca = 48)
fstclusters <- find.clusters.genlight(glfst,n.pca = 48)


fstclusters <- find.clusters.genlight(glfst,n.pca = 48, n.clust = 3)
tab<-data.frame(Sample=rownames(glmatrix),Cluster=fstclusters$grp, Landrace=subject.support[,'raza'],row.names = NULL)
tab<-tab[order(tab[,'Cluster'], tab[,'Landrace']),]

################################################################################################################################
## DAPC with high Fst SNPs
######### DAPC
maizclusters <- as.numeric(as.factor(subject.support[,'raza'])) ## define groups by landrace analysis
names(maizclusters)<-rownames(glmatrix)
maizclusters<-as.factor(maizclusters)

dapc.fst <- dapc.genlight(glfst,maizclusters,n.pca=4,n.da=2)
temp <- optim.a.score(dapc.fst) # analysis to select the correct number of PCA to use in the DAPC 

################################################################################################################################
### Altitude analysis.
## Classify samples according to altitude
subject.support<-cbind(subject.support,Altoobajo=cut(subject.support[,"Altitud"],breaks<-c(0,750,2241),labels=c("bajo","alto")))
### Fst measurement for every SNP among altitudinal categories, keeping SNP name and removing NAs
faltura<-Fst(snps,subject.support[,"Altoobajo"])
names(faltura$Fst)<-colnames(snps)
faltura$Fst<-faltura$Fst[!is.na(faltura$Fst)]
##### Subset SNP matrix keeping only altitude high Fst  SNPs 
altura.99<-fst.alta(faltura$Fst,.99)
snpsaltura99 <- names(altura.99$snp.info)
glmatrixalt <- glmatrix[, pmatch(snpsaltura99, colnames(glmatrix))]
glmatrixalt <- glmatrixalt[, colSums(is.na(glmatrixalt)) != nrow(glmatrixalt)]
############################################################################
## DAPC
glalt<- as.genlight(glmatrixalt)
set.seed(2)
altclusters<-find.clusters.genlight(glalt,n.clust=2,n.pca = 48)

dapc.alt <- dapc.genlight(glalt,altclusters$grp,n.pca=1, n.da=1)

################################################################################################################
## Figure_1 Maps

map <- get_googlemap( center = c(-97, 19), zoom = 6, maptype = "terrain", style = 'element:labels|visibility:off', color="bw")
p = ggmap(map)
p <- p + scale_y_continuous(limits = c(13.5,22.5))+ xlab("") + ylab("")

comiteco<- subset(subject.support,raza=="Comiteco")
conejo<- subset(subject.support,raza=="Conejo")
tehua<- subset(subject.support,raza=="Tehua")
zchico<- subset(subject.support,raza=="Zapalote Chico")
zgrande<- subset(subject.support,raza=="Zapalote Grande")

distcomi<-fortify(sh <- readShapePoly("/home/user/Escritorio/MAIZ/Articulo_tesis/mapas/Shapefiles_Ana/COMIT_7C.SHP"))
distcone<-fortify(sh <- readShapePoly("/home/user/Escritorio/MAIZ/Articulo_tesis/mapas/Shapefiles_Ana/CONEJ_6C.SHP"))
distte<-fortify(sh <- readShapePoly("/home/user/Escritorio/MAIZ/Articulo_tesis/mapas/Shapefiles_Ana/TEHUA_7C.SHP"))
distzc<-fortify(sh <- readShapePoly("/home/user/Escritorio/MAIZ/Articulo_tesis/mapas/Shapefiles_Ana/ZAPAC_6C.SHP"))
distzg<-fortify(sh <- readShapePoly("/home/user/Escritorio/MAIZ/Articulo_tesis/mapas/Shapefiles_Ana/ZAPAG_3C.SHP"))

col<-c("red","gold","green3","dodgerblue","magenta4")

com <- p +geom_polygon(aes(x=long,y=lat, group=group, alpha=0.05), data=distcomi, fill=col[1])+
  geom_point(data=comiteco, aes(x=Longitud, y=Latitud),size=1,color="black")+
  geom_point(data=comiteco, aes(x=Longitud, y=Latitud),size=.5,color=col[1])+
  theme(legend.position="none",plot.title = element_text(size=10),plot.margin=unit(c(0,0.3,-.3,0), "cm"),text = element_text(size=8))+ ggtitle("Comiteco")

con <- p + geom_polygon(aes(x=long,y=lat, group=group, alpha=0.5), data=distcone, fill=col[2])+
  geom_point(data=conejo, aes(x=Longitud, y=Latitud),size=1,color="black") +
  geom_point(data=conejo, aes(x=Longitud, y=Latitud),size=.5,color=col[2])+ 
  theme(legend.position="none",plot.title = element_text(size=10),plot.margin=unit(c(0,0.3,-.3,0), "cm"),text = element_text(size=8))+ ggtitle("Conejo")

teh <- p +  geom_polygon(aes(x=long,y=lat, group=group, alpha=0.5), data=distte, fill=col[3])+
  geom_point(data=tehua, aes(x=Longitud, y=Latitud),size=1) + 
  geom_point(data=tehua, aes(x=Longitud, y=Latitud),size=.5,color=col[3])+ 
  theme(legend.position="none",plot.title = element_text(size=10),plot.margin=unit(c(0,0.3,-.3,0), "cm"),text = element_text(size=8))+ ggtitle("Tehua")

zchi <- p + geom_polygon(aes(x=long,y=lat, group=group, alpha=0.5), data=distzc, fill=col[4]) +
  geom_point(data=zchico, aes(x=Longitud, y=Latitud),size=1) + 
  geom_point(data=zchico, aes(x=Longitud, y=Latitud),size=.5,color=col[4])+ 
  theme(legend.position="none",plot.title = element_text(size=10),plot.margin=unit(c(0,0.3,-.3,0), "cm"),text = element_text(size=8))+ ggtitle("Zapalote Chico")

zgra <- p +geom_polygon(aes(x=long,y=lat, group=group, alpha=0.5), data=distzg, fill=col[5]) +
  geom_point(data=zgrande, aes(x=Longitud, y=Latitud),size=1) + geom_point(data=zgrande, aes(x=Longitud, y=Latitud),size=.5,color=col[5])+
  theme(legend.position="none",plot.title = element_text(size=10),plot.margin=unit(c(0,0.3,-.3,0), "cm"),text = element_text(size=8))+ ggtitle("Zapalote Grande")

setEPS()
cairo_ps("Figure_1.eps",width=8,height = 10)
grid.arrange(com,con,teh,zchi,zgra,ncol=2,nrow=3)
dev.off()

############################################
## Figur2 2  SNP set
Fuente<- rownames(raza.99$fst.snp)%in%snpsextract[,2]

raza.99$grafica<-  raza.99$grafica + 
  geom_point(data=raza.99$fst.snp[Fuente,], aes(x=SNPs, y=Fst),size=2.5,shape=19,colour="black")+
  geom_point(data=raza.99$fst.snp[Fuente,], aes(x=SNPs, y=Fst,colour=Dif),size=1,shape=19)
postscript("Figure_2.eps",height=5,width=5)
  raza.99$grafica  
dev.off()
############################################################################################################3
### Figure_3
# PCA graphs
all<-pca.all.raza$grafica+ theme(legend.position="", plot.title = element_text(family="Arial",size=11))+ggtitle("All SNPs")+ xlab(paste("PC2 (",round(pca.all.raza$PCpercent,2)[2],"%)",sep="")) + ylab(paste("PC1 (",round(pca.all.raza$PCpercent,2)[1],"%)",sep=""))
dom<-pca.dom.raza$grafica+ theme(legend.position="", plot.title = element_text(family="Arial",size=11))+ggtitle(" Domestication and improvement SNPs")+ xlab(paste("PC2 (",round(pca.dom.raza$PCpercent,2)[2],"%)",sep="")) + ylab(paste("PC1 (",round(pca.dom.raza$PCpercent,2)[1],"%)",sep=""))
fst<-pca.fst99.raza$grafica +ggtitle(expression(paste("High F"[ST]," SNPs",sep="")))+theme(legend.position="", plot.title = element_text(family="Arial", face="bold",size=11))+ xlab(paste("PC2 (",round(pca.fst99.raza$PCpercent,2)[2],"%)",sep="")) + ylab(paste("PC1 (",round(pca.fst99.raza$PCpercent,2)[1],"%)",sep=""))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend <- g_legend(pca.fst99.raza$grafica)

postscript("Figure_3.eps",height=6,width=7)
grid.arrange(all,dom,fst,legend,ncol=2,nrow=2)
dev.off()
############################################################################################################3
### Figure_4
glmatrixdom <- glmatrix[, pmatch(snpsextract[,2],colnames(glmatrix))]
glmatrixdom <- glmatrixdom[, colSums(is.na(glmatrixdom)) != nrow(glmatrixdom)]

gldom <- as.genlight(glmatrixdom)
glfst <- as.genlight(glmatrixfst,centers=3)

fstclusters <- find.clusters.genlight(glfst,n.pca = 48)

cairo_ps("Figure_4.eps",height=8,width=3)
par(mfrow=c(3,1))
plot(allclusters$Kstat,xlab="Number of clusters", ylab="BIC",type="b",col="blue",font.main=1,main="All SNPs")
plot(domclusters$Kstat,xlab="Number of clusters", ylab="BIC",type="b",col="blue",font.main=1,main="Domestication and improvement SNPs")
plot(fstclusters$Kstat,xlab="Number of clusters", ylab="BIC",type="b",col="blue",font.main=1,main=expression(paste("High F"[ST]," SNPs",sep="")))
dev.off()
############################################################################################################3
### Figure_5
col<-c("red","gold","green3","dodgerblue","magenta4")

cairo_ps("Figure_5a.eps",height=3,width=7)
par(mfrow=c(1,1),cex=0.8)
scatter(dapc.fst,col=col ,bg="white",cstar=1,scree.da=FALSE, solid=.7, cex=1.5, clab=0, leg=TRUE, txt.leg=c("Comiteco","Conejo","Tehua", "Zapalote Chico", "Zapalote Grande"),posi.leg = "topleft",cleg=.8,cellipse=1,xlab="DF1")
dev.off()
## Ordenar las probabilidades de la función por raza para las gráficas
dapc.fst$posterior<-cbind(dapc.fst$posterior, subject.support[,"raza"])
dapc.fst$posterior<-dapc.fst$posterior[order(dapc.fst$posterior[,6]),]
colnames(dapc.fst$posterior)<-levels(subject.support[,"raza"])
dapc.fst$posterior<-dapc.fst$posterior[,-6]

cairo_ps("Figure_5b.eps",height=4,width=7)
par(mar=c(3.1, 4.1, 2.1, 2.1),mfrow=c(1,1),xpd=F,cex=0.8)
assignplot(dapc.fst,cex.lab=.6,pch="")
points(x=sort(as.integer(subject.support[,"raza"])),y=seq(49.5,.5,-1),pch=20) 
dev.off()

cairo_ps("Figure_5c.eps",height=4,width=7)
par(mar=c(4.1,4.1,5.1,.5),cex=0.8)
compoplot(dapc.fst,lab=rownames(dapc.fst$posterior), txt.leg=levels(subject.support[,"raza"]),cex.names=0.6,posi=1.22)
dev.off()

#####
cairo_ps("Figure_5.eps",height=10,width=6)
par(mar=c(5.1, 4.1, 4.1, 2.1),mfrow=c(3,1),cex=0.8,oma=c(0.2,0,1,0))
scatter(dapc.fst,col=col ,bg="white",cstar=1,scree.da=FALSE, solid=.7, cex=1.5, clab=0, leg=TRUE, txt.leg=c("Comiteco","Conejo","Tehua", "Zapalote Chico", "Zapalote Grande"),posi.leg = "topleft",cleg=.8,cellipse=1)
text("DF1",pos=4,x=-9.2,y=0.2,cex=0.8)
text("DF2",pos=4,x=-0.6,y=-4,cex=0.8,srt=90)
mtext(text="a)", side = 3, line = 0, outer=T,adj=0,cex=0.8)

par(xpd=F,mar=c(1.1, 4.1, 2.1, 2.1),cex=0.7)
assignplot(dapc.fst,cex.lab=.5,pch="")
points(x=sort(as.integer(subject.support[,"raza"])),y=seq(49.5,.5,-1),pch=20) 
par(xpd=T)
text("b)",pos=4,x=-.1,y=51,cex=1.1)


par(mar=c(4.1,4.1,5.1,.5),cex=0.8)
compoplot(dapc.fst,lab=rownames(dapc.fst$posterior), txt.leg=levels(subject.support[,"raza"]),cex.names=0.6,posi=1.22)
text("c)",pos=4,x=-11,y=1.3,cex=1.1)

dev.off()
################################################################################################################
## Supplementary figures. Landrace and altitude informative SNPs

################################################################################################################################
## Identify landrace informative SNPs

rownames(dapc.fst$var.contr)<-colnames(glmatrixfst)

## Supplementary figure 1: Landrace Informative SNPs
cairo_ps("Sup_Fig1.eps",height=9,width=6)
par(mfrow=c(3,2),mar=c(5.1,4.1,2.1,2.1),xpd=F)

plot(sort(dapc.fst$var.contr[,"LD1"]),xlab="SNPs",ylab="Loading", main="DF1")
abline(h=0.0095, col="red")

loci_x<-loadingplot (dapc.fst$var.contr,threshold = 0.0095, lab.jitter = 1, axis=1,xlab="SNPs",ylab="Loading",main="DF1")

plot(sort(dapc.fst$var.contr[,"LD2"]),xlab="SNPs",ylab="Loading", main="DF2")
abline(h=0.01464034 , col="red")
loci_y<-loadingplot (dapc.fst$var.contr,threshold = 0.015,cex.lab = 0.6, lab.jitter = 1, axis=2,xlab="SNPs",ylab="Loading",main="DF2")

boxplot(data.frame(DF1=dapc.fst$var.contr[,"LD1"],DF2=dapc.fst$var.contr[,"LD2"]),ylab="Loading",main="SNPs loading distribution")
lines(x=.5:1.5,y=rep(0.0095,2), col="red")
lines(x=1.5:2.5,y=rep(0.01173871 ,2), col="red")
dev.off()

### Altitude graphs
rownames(dapc.alt$posterior)<-make.names(subject.support[,'raza'],unique=TRUE)
dapc.alt$posterior<- dapc.alt$posterior[order(rownames(dapc.alt$posterior)),]
colnames(dapc.alt$posterior)<-c("Low","High")

## Altitude DAPC graphs
cairo_ps("Sup_Fig2.eps",height=7,width=8)
par(mar=c(5.1, 4.1, 4.1, 2.1),mfrow=c(1,2),xpd=T)
scatter(dapc.alt, legend=T, txt.leg=c("Low","High"),posi.leg="topleft")
text("a)",pos=4,x=-9,y=.38,cex=1)

par(xpd=F)
assignplot(dapc.alt,cex.lab = 0.6,pch=NA)
points(x=as.integer(subject.support[order(subject.support[,'raza']),"Altoobajo"]),y=seq(49.5,.5,-1),pch=20)
par(xpd=T)
text("b)",pos=4,x=-.1,y=53,cex=1)
dev.off()

rownames(dapc.alt$var.contr) <- colnames(glmatrixalt)

cairo_ps("Sup_Fig3.eps",height=9,width=8)
par(mfrow=c(2,2))
plot(sort(dapc.alt$var.contr[,"LD1"]),xlab="SNPs",ylab="Loading")
abline(h=0.012, col="red")
loadingplot(dapc.alt$var.contr,threshold = 0.012,lab.jitter=1,xlab="SNPs",ylab="Loading",main="")

boxplot(data.frame(dapc.alt$var.contr),ylab="Loading",main="SNPs loading distribution",xlab="SNPs")
dev.off()
### Check for SNPs shared between the landrace and altitude informative SNPs
sum(snpsalt$var.names%in%loci_x$var.names)
sum(snpsalt$var.names%in%loci_y$var.names)

snpsalt$var.names[(snpsalt$var.names%in%loci_x$var.names)]
