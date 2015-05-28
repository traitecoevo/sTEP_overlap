library(data.table)

get_cont_raster<-function(){
  library(sp) #Load your libraries
  library(maptools)
  library(raster)
  #Download the continents shapefile
  download.file("http://baruch.cuny.edu/geoportal/data/esri/world/continent.zip",
                "cont.zip")
  #Unzip it
  unzip("cont.zip")
  #Load it
  cont <- readShapeSpatial("continent.shp")
  worldclim<-getData('worldclim', var='tmean', res=2.5)
  cont.raster<-rasterize(cont,worldclim,field="CONTINENT")
  writeRaster(cont.raster,"cont2.5.grd")
}
library(dplyr)


a<-fread("saving_to_be_safe.csv")
cont<-raster("cont2.5.grd")
a<-filter(a,!is.na(as.numeric(decimalLongitude)))

sp<-SpatialPoints(cbind(as.numeric(a$decimalLongitude),as.numeric(a$decimalLatitude)))

cont.x<-extract(cont,sp)
a$cont<-cont.x

endemic<-summarize(group_by(a,species),num.cont=length(unique(cont)))
table(endemic$num.cont)
one.cont.endemics<-filter(endemic,num.cont==1)

filter(a,species%in%one.cont.endemics$species)%>%
  dplyr::select(species,cont)->one.cont.list

one.cont.list<-one.cont.list[!duplicated(one.cont.list$species),]

filter(one.cont.list,cont==4)%>%
  dplyr::select(species)->known.aussie.endemics

filter(one.cont.list,cont==7)%>%
  dplyr::select(species)->known.oceania
