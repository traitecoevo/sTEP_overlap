# Libraries
library(sp)
library(maptools)
library(raster)

get_cont_raster<-function(){
  cont.shp <- readShapeSpatial("raw_data/continent.shp")
  worldclim<-raster::getData('worldclim', var='tmin', res=2.5)
  worldclim10<-raster::getData('worldclim', var='tmin', res=10)
  
  cont.raster<-rasterize(cont.shp,worldclim,field="CONTINENT")
  writeRaster(cont.raster,"clean_data/cont2.5.grd")
}

gbif_tpl<-function(gbif){
  syns<-fread("raw_data/tpl_names.txt")
  syns$correct.names<-sub("_"," ",syns$correct.names)
  syns$all.names<-sub("_"," ",syns$all.names)
  gbif$corrected.name<-syns$correct.names[match(gbif$species,syns$all.names)]
  gbif<-filter(gbif,!is.na(corrected.name))
  gbif<-select(gbif,species=corrected.name,lat=decimalLatitude,long=decimalLongitude)
  return(gbif)
}


get_gbif_names<-function(){
  a <- fread("clean_data/gbif_spp.txt")
  names(a)<-c("species","lat","long")
  read_csv("raw_data/tpl_names.txt")%>%
    dplyr::select(correct.names)%>%
    mutate(correct.names=tolower(gsub("_", " ", correct.names)))->goodNames
    accNames<-unique(goodNames$correct.names)
  b<-filter(a,species%in%accNames) 

  return(b)
}

get_gbif<-function(){
  a <- fread("clean_data/gbif_tpl_locations.csv")
  names(a)<-c("species","lat","long")
  read_csv("raw_data/tpl_names.txt")%>%
    dplyr::select(correct.names)%>%
    mutate(correct.names=tolower(gsub("_", " ", correct.names)))->goodNames
    accNames<-unique(goodNames$correct.names)
  b<-filter(a,species%in%accNames) 

  return(b)
}

get_genbank<-function(){
  genbank<-read.csv("clean_data/genbank_spp_clean.txt",header=FALSE,as.is=TRUE)
  genbank.scrubbed<-tolower(genbank$V1)
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


do.gam.analysis<-function(b,type){
  genbank.scrubbed<-get_genbank()
  b$genbank.yes.no<-b$species%in%genbank.scrubbed
  ou<-makeCluster(15,type="SOCK")
  gam.genbank<-bam(genbank.yes.no~s(lat),family=binomial(),data=b,cluster=ou,gc.level=2)
  try_sp<-read.csv("clean_data/try_spp_clean.txt", header=FALSE, as.is=TRUE)[,1]
  #try_sp$sp_scrubb<-scrub(try_sp)
  b$try.yes.no<-b$species%in%tolower(try_sp)
  out_try<-bam(try.yes.no~s(lat),family=binomial(),data=b,cluster=ou,gc.level=2)
 
  b$try.genbank.yes.no<-b$species%in%tolower(try_sp)&
                        b$species%in%genbank.scrubbed
 
   gam_try_plus_genbank<-bam(try.genbank.yes.no~s(lat),family=binomial(),data=b,cluster=ou,gc.level=2)
 
  gam.df<-data.frame(lat=c(b$lat,b$lat,b$lat),
                    fit=c(fitted(out_try),
                    fitted(gam.genbank),
                    fitted(gam_try_plus_genbank)),
                    dataset=c(rep("TRY",length(b$lat)),
                    rep("genbank",length(b$lat)),
                    rep("genbank_plus_try",length(b$lat))),
                    type=type)
  stopCluster(ou)
  return(gam.df)
}


plot_gbif_sampling<-function(){
   b<-get_gbif()
   genbank.scrubbed<-get_genbank()
  b$genbank.yes.no<-b$species%in%genbank.scrubbed
  ou<-makeCluster(15,type="SOCK")
  try_sp<-read.csv("clean_data/try_spp_clean.txt", header=FALSE, as.is=TRUE)[,1]
  #try_sp$sp_scrubb<-scrub(try_sp)
  b$try.yes.no<-b$species%in%tolower(try_sp)
 
  b$try.genbank.yes.no<-b$species%in%tolower(try_sp)&
                        b$species%in%genbank.scrubbed
  group_by(b,species)%>%
  summarize(number_of_gbif_obs=n(),try_presence=mean(try.yes.no),genbank_presence=mean(genbank.yes.no))->z

  z$database_status <- ifelse(z$try_presence==1&z$genbank_presence==0, "Only Try",
                              ifelse(z$try_presence==0&z$genbank_presence==1, "Only Genbank",
                                            ifelse(z$try_presence==1&z$genbank_presence==1, "Both Try and Genbank",
                                                   ifelse(z$try_presence==0&z$genbank_presence==0, "Neither Try nor Genbank","error"
                                                   ))))

  z$database_status <- factor(z$database_status,levels = c("Both Try and Genbank", "Only Try", "Only Genbank", "Neither Try nor Genbank"))
     pdf("figures/hist_gbif.pdf",width=8.5,height=5)
  print(ggplot(z,aes(x=number_of_gbif_obs,fill=database_status))+
          geom_histogram(position = "stack", binwidth=0.4)+
          scale_x_log10()+
        theme_classic()+
        scale_fill_brewer(palette="Set2",type="qual")
        )
  dev.off()                    

}


run_gam_df<-function(){
  a<-get_gbif()
  
  #random subsample
  random.obs<-a[sample(1:dim(a)[1],2*10^7,replace=F),]
  
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
  ss<-out[sample(1:dim(out)[1],1*10^6,replace=F),]
  ggplot(ss,aes(x=lat,y=fit))+
          ylab("Proportion in database")+
          geom_line(aes(col=dataset,linetype=type))+
        theme_classic()
  ggsave("figures/multi_gam.png",width=8.5,height=5)
}

plot_gbif_bins<-function(out){
  #plot
  ss<-out[sample(1:dim(out)[1],1*10^6,replace=F),]
  png("figures/multi_gam.png",width=8.5,height=5)
  print(ggplot(ss,aes(x=lat,y=fit))+
          ylab("Proportion in database")+
          geom_line(aes(col=dataset,linetype=type))+
        theme_classic())
  dev.off()
}

plot_gam_2<-function(gam_df){
gam_df<-filter(gam_df,abs(lat)<67)
ss<-gam_df[sample(1:dim(gam_df)[1],1*10^6,replace=F),]

ggplot(ss,aes(x=lat,y=fit))+
          ylab("Proportion in database")+
          geom_line(aes(col=dataset,linetype=type))+
        theme_classic()
ggsave("figures/multi_gam_2.png",width=8.5,height=5)
}

mean_gbif<-function(a){
  
  png("figures/multi_gam_species_based.png")
  print(ggplot(gam.df,aes(x=lat,y=fit))+geom_point(aes(col=dataset)))
  dev.off()
}



add_continent<-function(){
  #NOT WORKING YET
  cont<-raster("clean_data/cont2.5.grd")
  a <- fread("clean_data/gbif_tpl_locations.txt")
  sp<-SpatialPoints(cbind(as.numeric(a$long),as.numeric(a$lat)))
  cont.x<-extract(cont,sp)
  a$cont<-cont.x
  write_csv(a,"clean_data/gbif_cleaner_geo_data.csv")
  return(table(a$cont))
}

other_stuff<-function(){
  #get_cont_raster()
  cont<-raster("clean_data/cont2.5.grd")
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
  a<-read_csv("clean_data/gbif_tpl_locations.txt")
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

#check_endemics<-function(){
#  one.cont.list<-read_csv("one_cont_list.csv")
#  genbank.scrubbed<-get_genbank()
#  one.cont.list$one.cont.list.genbank.yes.no<-one.cont.list$species%in%genbank.scrubbed
  #tt<-table(one.cont.list$cont,one.cont.list$one.cont.list.yes.no)
  #tt[,2]/(tt[,1]+tt[,2])
#  try_sp<-read_csv("TryAccSpecies.txt",col_names="AccSpeciesName")
#  try_sp$sp_scrubb<-try_sp
#  one.cont.list$try.yes.no<-one.cont.list$species%in%try_sp$sp_scrubb
  #tt<-table(one.cont.list$cont,one.cont.list$try.yes.no)
  #tt[,2]/(tt[,1]+tt[,2])
#  write_csv(one.cont.list,"one_cont_list.csv")
#}

