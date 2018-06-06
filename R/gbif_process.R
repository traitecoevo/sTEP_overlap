

get_cont_raster<-function(){
  library(sp) #Load your libraries
  library(maptools)
  library(raster)
  #Download the continents shapefile
  download.file("http://faculty.baruch.cuny.edu/geoportal/data/esri/world/continent.zip",
                "cont.zip")
  
  #Unzip it
  unzip("cont.zip")
  #Load it
  cont.shp <- readShapeSpatial("continent.shp")
  worldclim<-getData('worldclim', var='tmean', res=2.5)
  worldclim10<-getData('worldclim', var='tmean', res=10)
  
  cont.raster<-rasterize(cont.shp,worldclim,field="CONTINENT")
  writeRaster(cont.raster,"cont2.5.grd")
}

gbif_tpl<-function(gbif){
  syns<-fread("../../../srv/scratch/z3484779/taxonomicResources/plantList11syns.csv")
  syns$correct.names<-sub("_"," ",syns$correct.names)
  syns$all.names<-sub("_"," ",syns$all.names)
  gbif$corrected.name<-syns$correct.names[match(gbif$species,syns$all.names)]
  gbif<-filter(gbif,!is.na(corrected.name))
  gbif<-select(gbif,species=corrected.name,lat=decimalLatitude,long=decimalLongitude)
  return(gbif)
}

get_gbif<-function(){
  if (Sys.info()[[1]]=="Linux") a<-fread("../../../srv/scratch/z3484779/overlap-data/raw_data/cooked.csv")
  else a<-fread("occur.csv")
  names(a)<-c("species","lat","long")
  read_csv("../../../srv/scratch/z3484779/taxonomicResources//plantList11syns.csv")%>%
    dplyr::select(correct.names)%>%
    mutate(correct.names=tolower(gsub("_", " ", correct.names)))->goodNames
    accNames<-unique(goodNames$correct.names)
  b<-filter(a,species%in%accNames) 

  return(b)
}

get_genbank<-function(){
  genbank<-read.csv("genBankList.txt",header=FALSE,as.is=TRUE)
  genbank.scrubbed<-tolower(scrub(genbank$V1))
  out<-genbank.scrubbed[!is.na(genbank.scrubbed)]
  return(out)
}

make_sampling_map<-function(a){
  genbank.scrubbed<-get_genbank()
  a$genbank.yes.no<-a$species%in%genbank.scrubbed
  worldclim10<-getData('worldclim', var='tmean', res=10)
  sp<-SpatialPoints(cbind(as.numeric(a$decimalLongitude),as.numeric(a$decimalLatitude)))
  sp<-SpatialPointsDataFrame(coords=sp,data=data.frame(ing=a$genbank.yes.no))
  genbank.sampling.map<-rasterize(sp,worldclim10,field="ing",fun=mean)
  pdf("figures/genbank.pdf")
  plot(genbank.sampling.map,col=brewer.pal(9,"Blues"),main="Proportion of GBIF observations that are in GenBank")
  #plot(cont.shp,add=TRUE,lwd = 0.3)
  dev.off()
}

zae.diaz.analysis<-function(b,type){
  zae<-read.tree("zanne_tpl_1_1.tre")
  z<-scrub(zae$tip.label)
  b$genbank.yes.no<-b$species%in%z
  ou<-makeCluster(15,type="SOCK")
  gam.genbank<-bam(genbank.yes.no~s(lat),family=binomial(),data=b,cluster=ou,gc.level=2)
  #try_sp<-read.delim("TryAccSpecies.txt",as.is=TRUE)
  #try_super<-try_sp$AccSpeciesName[which(try_sp$TraitNum>5)]
  #ts<-scrub(try_super)
  diaz<-read_csv("diaz_etal_names.csv")$name_TLP_TRY30_resolved
  diaz<-use.synonym.lookup(scrub(diaz))
  b$try.yes.no<-b$species%in%diaz
  out_try<-bam(try.yes.no~s(lat),family=binomial(),data=b,cluster=ou,gc.level=2)
  
  gam.df<-data.frame(lat=c(b$lat,b$lat),fit=c(fitted(out_try),fitted(gam.genbank)),dataset=c(rep("Well sampled TRY",length(b$lat)),rep("Zanne",length(b$lat))),type=type)
  stopCluster(ou)
  return(gam.df)
}

do.gam.analysis<-function(b,type){
  genbank.scrubbed<-get_genbank()
  b$genbank.yes.no<-b$species%in%genbank.scrubbed
  ou<-makeCluster(15,type="SOCK")
  gam.genbank<-bam(genbank.yes.no~s(lat),family=binomial(),data=b,cluster=ou,gc.level=2)
  try_sp<-read.delim("TryAccSpecies.txt",as.is=TRUE)
  #try_sp$sp_scrubb<-scrub(try_sp$AccSpeciesName)
  b$try.yes.no<-b$species%in%tolower(try_sp$AccSpeciesName)
  out_try<-bam(try.yes.no~s(lat),family=binomial(),data=b,cluster=ou,gc.level=2)
  gam.df<-data.frame(lat=c(b$lat,b$lat),fit=c(fitted(out_try),fitted(gam.genbank)),dataset=c(rep("TRY",length(b$lat)),rep("genbank",length(b$lat))),type=type)
  stopCluster(ou)
  return(gam.df)
}


plot_gbif_bins<-function(){
  a<-get_gbif()
  
  #random subsample
  random.obs<-a[sample(1:dim(a)[1],2*10^6,replace=F),]
  
  #species median dataset
  by.species<-summarize(group_by(a,species),lat=median(lat))
  
  #memmory clean up
  rm(a)
  gc()
  
  #run the gams
  gam.df.obs<-do.gam.analysis(random.obs,type="by gbif observation")
  gam.df.sp<-do.gam.analysis(by.species,type="by species")
 
  #stick data frame together
  out<-rbind(gam.df.obs,gam.df.sp)
  
  #plot
  png("figures/multi_gam.png")
  print(ggplot(out,aes(x=lat,y=fit))+
          ylab("Proportion in database")+
          geom_line(aes(col=dataset,linetype=type))+
        theme_classic())
  dev.off()
}



mean_gbif<-function(a){
  
  png("figures/multi_gam_species_based.png")
  print(ggplot(gam.df,aes(x=lat,y=fit))+geom_point(aes(col=dataset)))
  dev.off()
}



add_continent<-function(){
  #NOT WORKING YET
  cont<-raster("cont2.5.grd")
  if (Sys.info()[[1]]=="Linux") a<-fread("../../../srv/scratch/z3484779/gbif/gbif_cleaner.csv")
  else a<-fread("occur.csv")
  sp<-SpatialPoints(cbind(as.numeric(a$long),as.numeric(a$lat)))
  cont.x<-extract(cont,sp)
  a$cont<-cont.x
  write_csv(a,"../../../srv/scratch/z3484779/gbif/gbif_cleaner_geo_data.csv")
  return(table(a$cont))
}

other_stuff<-function(){
  #get_cont_raster()
  cont<-raster("cont2.5.grd")
  #a<-filter(a,!is.na(as.numeric(decimalLongitude)))
  a$in.genbank<-a$species%in%genbank
  a$in.try<-a$species%in%try.sp
  
  
  sp<-SpatialPoints(cbind(as.numeric(a$decimalLongitude),as.numeric(a$decimalLatitude)))
  sp<-SpatialPointsDataFrame(coords=sp,data=data.frame(ing=a$in.genbank,int=a$in.try))
  genbank.sampling.map<-rasterize(sp,worldclim10,field="ing",fun=mean)
  try.sampling.map<-rasterize(sp,worldclim10,field="int",fun=mean)
  
  pdf("genbank.pdf")
  plot(genbank.sampling.map,col=brewer.pal(9,"Blues"),main="Proportion of GBIF observations that are in GenBank")
  #plot(cont.shp,add=TRUE,lwd = 0.3)
  dev.off()
  
  pdf("try-genbank.pdf")
  plot(try.sampling.map,col=brewer.pal(9,"Blues"),main="Proportion of GBIF observations that are in TRY")
  #plot(cont.shp,add=TRUE,lwd = 0.3)
  dev.off()
  
  sp<-SpatialPoints(cbind(as.numeric(a$decimalLongitude),as.numeric(a$decimalLatitude)))
  
  cont.x<-extract(cont,sp)
  a$cont<-cont.x
}
  
get_endemics<-function(){
  a<-read_csv("../../../srv/scratch/z3484779/gbif/gbif_cleaner_geo_data.csv")
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
  write_csv(one.cont.list,"one_cont_list.csv")
}

check_endemics<-function(){
  one.cont.list<-read_csv("one_cont_list.csv")
  genbank.scrubbed<-get_genbank()
  one.cont.list$one.cont.list.genbank.yes.no<-one.cont.list$species%in%genbank.scrubbed
  tt<-table(one.cont.list$cont,one.cont.list$one.cont.list.yes.no)
  tt[,2]/(tt[,1]+tt[,2])
  try_sp<-read.delim("TryAccSpecies.txt",as.is=TRUE)
  try_sp$sp_scrubb<-scrub(try_sp$AccSpeciesName)
  one.cont.list$try.yes.no<-one.cont.list$species%in%try_sp$sp_scrubb
  tt<-table(one.cont.list$cont,one.cont.list$try.yes.no)
  tt[,2]/(tt[,1]+tt[,2])
  write_csv(one.cont.list,"one_cont_list.csv")
}

