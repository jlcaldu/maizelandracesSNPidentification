require(maptools)
require(ggmap)
require(gridExtra)

setwd("/home/user/Escritorio/MAIZ/Articulo_tesis/Datos_Maiz/DATA")

map <- get_googlemap( center = c(-97, 19), zoom = 6, maptype = "terrain", style = 'element:labels|visibility:off', color="bw")
p = ggmap(map)
p <- p + scale_y_continuous(limits = c(13.5,22.5))+ xlab("") + ylab("")
p

subject.support=read.csv("subject.support.csv")
subject.support= subject.support[order(subject.support[,"raza"]),c("raza","Longitud","Latitud")]

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

setwd("~/Escritorio/MAIZ/Articulo_tesis/mapas")
setEPS()
cairo_ps("Figure1.eps",width=8,height = 10)
grid.arrange(com,con,teh,zchi,zgra,ncol=2,nrow=3)
dev.off()
