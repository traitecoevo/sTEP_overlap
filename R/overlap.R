# library(dplyr)
# library(Hmisc)
# library(ape)
# library(data.table)
# source("data_manipulation_functions.R")

# this is the genbank species list from the NCBI Browser website
#read.delim("genBankList.txt",header=FALSE,as.is=T)%>%
#  mutate(V2=scrub(V1))%>%
#  mutate(V3=use.synonym.lookup(V2))->genbank.all
#genbank<-unique(genbank.all$V3)

# this is the genbank species list from Will's parsing of the file from the ftp server
#read.delim("names_subset.txt",as.is=T,header=F) %>%
#  mutate(V2=scrub(V1))%>%
#  mutate(V3=use.synonym.lookup(V2))->genbank.2
#genbank2<-unique(genbank.2$V3)

# this is the species list from try
#read.delim("TryAccSpecies.txt",as.is=TRUE)%>%
#  dplyr::select(AccSpeciesName)%>%
#  mutate(sp=scrub(AccSpeciesName))%>%
#  mutate(sp.fix=use.synonym.lookup(sp))->try.all
#try.sp<-unique(try.all$sp.fix)

#gets the accepted names from plant list 1.1
#read.csv("taxonomicResources/plantList11syns.csv")%>%
#  select(correct.names)%>%
#  mutate(correct.names=gsub("_", " ", correct.names))->pl
#pl<-unique(pl$correct.names)

#gets the names from the zanne tree
#tree <- read.tree("Vascular_Plants_rooted.dated.tre")
#tree$tip.label <- gsub("_", " ", tree$tip.label)
#tree.names<-scrub(tree$tip.label)
#tree.names<-use.synonym.lookup(tree.names)

get_gbif_names<-function(){
  a<-get_gbif()
  return(unique(a$species))
}


#thanks to peeps on stackoverflow
calculate.overlap <- function(y, i){
  pw <- combn(seq_along(y), i, FUN= function(x) {
    res <- length(Reduce(intersect, y[x]))
    names(res) <- paste(names(y[x]), collapse = "-")
    res
  }, simplify = FALSE)
  return(do.call(c, pw))
}


do_overlap_analysis<-function(){
y<-list(genbank=get_genbank(),
        try.all.names=read.in.try(),
        gbif=get_gbif_names())
out<-c(sapply(y,length),calculate.overlap(y,2),calculate.overlap(y,3))
write.csv(out,"tables/two_and_three_way_comparisons.csv")
}

#FOR FUTURE CHECKS:
#tree.names[!tree.names%in%genbank]
#good.names<-sapply(as.character(syn$V1),FUN=function(x) strsplit(x,",")[[1]][1],USE.NAMES=F)

