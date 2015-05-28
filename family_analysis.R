library(data.table)
library(TaxonLookup)
library(Deducer)
library(dplyr)
library(taxize)
library(parallel)
library(stargazer)
source("data_manipulation_functions.R")


sync.species.lists<-function(sampled.list){
  fread("taxonomicResources/plantList11syns.csv")%>%
    dplyr::select(correct.names)%>%
    mutate(correct.names=gsub("_", " ", correct.names))->goodNames
  goodNames<-data.frame(gs=unique(goodNames$correct.names),stringsAsFactors=FALSE)
  goodNames$genera<-sapply(as.character(goodNames$gs),FUN=function(x) strsplit(x," ")[[1]][1],USE.NAMES=F)
  pl<-plant_lookup()
  goodNames$family<-pl$family[match(goodNames$genera,pl$genus)]
  goodNames$family[is.na(goodNames$family)]<-"bryo"
  goodNames$in.list<-goodNames$gs%in%sampled.list
  return(goodNames)
}

prepare.sampling.df<-function(sampled.list,ref.list){
  goodNames<-data.frame(gs=unique(ref.list),stringsAsFactors=FALSE)
  goodNames$genera<-sapply(as.character(goodNames$gs),FUN=function(x) strsplit(x," ")[[1]][1],USE.NAMES=F)
  pl<-plant_lookup()
  goodNames$family<-pl$family[match(goodNames$genera,pl$genus)]
  goodNames$in.list<-goodNames$gs%in%sampled.list
  goodNames<-filter(goodNames,!is.na(family))
  return(goodNames)
}




test.family<-function(family.in,goodNames=goodNames){
  #print(family.in)
  goodNames$family.of.interest<-goodNames$family==family.in
  if(sum(goodNames$family.of.interest,na.rm=T)>0){
    return(likelihood.test(goodNames$in.list,goodNames$family.of.interest))
  }
  return(NA)
}

calc.proportion<-function(family.in,goodNames=goodNames){
  #print(family.in)
  fam.only<-filter(goodNames,family==family.in)
  return(sum(fam.only$in.list)/length(fam.only$in.list))
}

read.in.try<-function(){
  require(dplyr)
  fread("TryAccSpecies.txt")%>%
    dplyr::select(AccSpeciesName)%>%
    mutate(sp=scrub(AccSpeciesName))%>%
    mutate(sp.fix=use.synonym.lookup(sp))->try.all
  try.sp<-unique(try.all$sp.fix)
  return(try.sp)
}

process.endemic.list<-function(sp.names){
  sp.names%>%
    mutate(species=scrub(species))%>%
    mutate(sp.fix=use.synonym.lookup(species))->endem.fixed
  endemic.out<-unique(endem.fixed$sp.fix)
  return(endemic.out)
}

read.genBank<-function(){
# this is the genbank species list from the NCBI Browser website
read.delim("genBankList.txt",header=FALSE,as.is=T)%>%
  mutate(V2=scrub(V1))%>%
  mutate(V3=use.synonym.lookup(V2))->genbank.all
genbank<-unique(genbank.all$V3)
return(genbank)
}

read.csv("")

# fread("species_centers.txt")%>%
#   dplyr::select(sp=V2)%>%
#   mutate(sp=scrub(sp))%>%
#   mutate(sp=use.synonym.lookup(sp))->gbif
#goodNames<-sync.species.lists(gbif$sp)

genbank<-read.genBank()
oceania<-process.endemic.list(known.oceania)
oceania.try<-prepare.sampling.df(sampled.list = try.sp,ref.list = oceania)
aussie.try<-prepare.sampling.df(sampled.list = try.sp,ref.list = aussie)
ocenaia.genbank<-prepare.sampling.df(sampled.list=genbank,ref.list=oceania)
#goodNames<-sync.species.lists(try.sp)

oceania.try<-filter(oceania.try,!is.na(family))
family.list<-as.list(unique(oceania.try$family))
test<-lapply(family.list,FUN=test.family,goodNames=ocenaia.genbank)

prop<-unlist(lapply(family.list,FUN=calc.proportion,goodNames=aussie.try))
#sr<-summarize(group_by(goodNames,family),length(gs))
g<-unlist(lapply(test,function(x)x$statistic))
p<-unlist(lapply(test,function(x)x$p.value))
sr<-unlist(lapply(test,function(x)sum(x$observed[,2])))

a<-test.family("Euphorbiaceae",ocenaia.genbank)
calc.proportion("Euphorbiaceae",ocenaia.genbank)

#NOT WORKING NOW
data.table(family=unlist(family.list),prop.sampled=prop,sr=sr,g=g,p=p)%>%
  arrange(g)->ranking

mean(ranking$prop.sampled,na.rm=T)

filter(ranking,prop.sampled<mean(prop.sampled,na.rm=T))%>%
  arrange(desc(g))->oceania.try.fam

stargazer(aussie.try.fam[1:5,], summary=FALSE)
