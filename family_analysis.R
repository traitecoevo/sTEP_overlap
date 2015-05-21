library(data.table)
library(TaxonLookup)
library(Deducer)
library(dplyr)
library(taxize)
library(parallel)
source("data_manipulation_functions.R")


sync.species.lists<-function(sampled.list){
  fread("taxonomicResources/plantList11syns.csv")%>%
    select(correct.names)%>%
    mutate(correct.names=gsub("_", " ", correct.names))->goodNames
  goodNames<-data.frame(gs=unique(goodNames$correct.names),stringsAsFactors=FALSE)
  goodNames$genera<-sapply(as.character(goodNames$gs),FUN=function(x) strsplit(x," ")[[1]][1],USE.NAMES=F)
  pl<-plant_lookup()
  goodNames$family<-pl$family[match(goodNames$genera,pl$genus)]
  goodNames$family[is.na(goodNames$family)]<-"bryo"
  goodNames$in.list<-goodNames$gs%in%sampled.list
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

read.delim("TryAccSpecies.txt",as.is=TRUE)%>%
  select(AccSpeciesName)%>%
  mutate(sp=scrub(AccSpeciesName))%>%
  mutate(sp.fix=use.synonym.lookup(sp))->try.all
try.sp<-unique(try.all$sp.fix)
rm(try.all)

goodNames<-sync.species.lists(try.sp)

fread("species_centers.txt")%>%
  select(sp=V2)%>%
  mutate(sp=scrub(sp))%>%
  mutate(sp=use.synonym.lookup(sp))->gbif
goodNames<-sync.species.lists(gbif$sp)

family.list<-as.list(unique(goodNames$family))
test<-mclapply(family.list,FUN=test.family,goodNames=goodNames)



prop<-unlist(mclapply(family.list,FUN=calc.proportion,goodNames=goodNames))
#sr<-summarize(group_by(goodNames,family),length(gs))
g<-unlist(lapply(test,function(x)x$statistic))
p<-unlist(lapply(test,function(x)x$p.value))
sr<-unlist(lapply(test,function(x)sum(x$observed[,2])))

ranking<-data.frame(family=unique(goodNames$family),prop.sampled=prop,sr=sr,g=g,p=p,stringsAsFactors=FALSE)
write.csv(ranking,"temp.csv")


