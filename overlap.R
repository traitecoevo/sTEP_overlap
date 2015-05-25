library(dplyr)
library(Hmisc)
library(ape)
source("data_manipulation_functions.R")

# this is the genbank species list from the NCBI Browser website
read.delim("genBankList.txt",header=FALSE,as.is=T)%>%
  mutate(V2=scrub(V1))%>%
  mutate(V3=use.synonym.lookup(V2))->genbank.all
genbank<-unique(genbank.all$V3)

# this is the genbank species list from Will's parsing of the file from the ftp server
read.delim("names_subset.txt",as.is=T,header=F) %>%
  mutate(V2=scrub(V1))%>%
  mutate(V3=use.synonym.lookup(V2))->genbank.2
genbank2<-unique(genbank.2$V3)

# this is the species list from try
read.delim("TryAccSpecies.txt",as.is=TRUE)%>%
  select(AccSpeciesName)%>%
  mutate(sp=scrub(AccSpeciesName))%>%
  mutate(sp.fix=use.synonym.lookup(sp))->try.all
try.sp<-unique(try.all$sp.fix)

#gets the accepted names from plant list 1.1
read.csv("taxonomicResources/plantList11syns.csv")%>%
  select(correct.names)%>%
  mutate(correct.names=gsub("_", " ", correct.names))->pl
pl<-unique(pl$correct.names)

#gets the names from the zanne tree
tree <- read.tree("Vascular_Plants_rooted.dated.tre")
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree.names<-scrub(tree$tip.label)
tree.names<-use.synonym.lookup(tree.names)

#clean up a bit
rm(list=c("try.all","tree","genbank.2","genbank.all"))
gc()   

#if i were a better programer these two would be one function
calculate.overlap<-function(y){
  pw<-combn(y,2,FUN=function(x)sum(x[[1]]%in%x[[2]]))
  names(pw)<-combn(y,2,FUN=function(x)paste(names(x)[[1]],names(x)[[2]],sep="-"))
  return(pw)
}
calculate.overlap.3<-function(y){
  pw<-combn(y,3,FUN=function(x)sum(x[[1]]%in%x[[2]]&x[[1]]%in%x[[3]]))
  names(pw)<-combn(y,3,FUN=function(x)paste(names(x)[[1]],names(x)[[2]],names(x)[[3]],sep="-"))
  return(pw)
}


y<-list(zanne.tree=tree.names,try.all.names=try.sp,genbank.parse1=genbank,genbank.parse2=genbank2,plant.list.1.1=pl)
out<-c(sapply(y,length),calculate.overlap(y),calculate.overlap.3(y))
write.csv(out,"output/two_and_three_way_comparisons.csv")

#FOR FUTURE CHECKS:
#tree.names[!tree.names%in%genbank]
#good.names<-sapply(as.character(syn$V1),FUN=function(x) strsplit(x,",")[[1]][1],USE.NAMES=F)

