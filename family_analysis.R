library(TaxonLookup)
library(Deducer)
library(dplyr)
library(taxize)
source("data_manipulation_functions.R")

read.csv("taxonomicResources/plantList11syns.csv")%>%
  select(correct.names)%>%
  mutate(correct.names=gsub("_", " ", correct.names))->goodNames
goodNames<-data.frame(gs=unique(goodNames$correct.names),stringsAsFactors=FALSE)

goodNames$genera<-sapply(as.character(goodNames$gs),FUN=function(x) strsplit(x," ")[[1]][1],USE.NAMES=F)
pl<-plant_lookup()
goodNames$family<-pl$family[match(goodNames$genera,pl$genus)]
goodNames$family[is.na(goodNames$family)]<-"bryo"


read.delim("TryAccSpecies.txt",as.is=TRUE)%>%
  select(AccSpeciesName)%>%
  mutate(sp=scrub(AccSpeciesName))%>%
  mutate(sp.fix=use.synonym.lookup(sp))->try.all
try.sp<-unique(try.all$sp.fix)
rm(try.all)

goodNames$in.try<-goodNames$gs%in%try.sp

test.family<-function(family.in,goodNames=goodNames){
  print(family.in)
  goodNames$family.of.interest<-goodNames$family==family.in
  if(sum(goodNames$family.of.interest,na.rm=T)>0){
  return(likelihood.test(goodNames$in.try,goodNames$family.of.interest))
  }
  return(NA)
}

calc.proportion<-function(family.in,goodNames=goodNames){
  #print(family.in)
  fam.only<-filter(goodNames,family==family.in)
  return(sum(fam.only$in.try)/length(fam.only$in.try))
}

family.list<-as.list(unique(goodNames$family))
test<-lapply(family.list,FUN=test.family,goodNames=goodNames)
prop<-unlist(lapply(family.list,FUN=calc.proportion,goodNames=goodNames))
#sr<-summarize(group_by(goodNames,family),length(gs))
g<-unlist(lapply(test,function(x)x$statistic))
p<-unlist(lapply(test,function(x)x$p.value))
sr<-unlist(lapply(test,function(x)sum(x$observed[,2])))

ranking<-data.frame(family=unique(goodNames$family),prop.sampled=prop,sr=sr,g=g,p=p,stringsAsFactors=FALSE)
write.csv(ranking,"temp.csv")

pchisq(8,1)

for(i in 1:length(test)){
  print(i)
  
}
