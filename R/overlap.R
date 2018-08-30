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



#thanks to peeps on stackoverflow
calculate.overlap <- function(y, i){
  pw <- combn(seq_along(y), i, FUN= function(x) {
    res <- length(Reduce(intersect, y[x]))
    names(res) <- paste(names(y[x]), collapse = "-")
    res
  }, simplify = FALSE)
  return(do.call(c, pw))
}

read.in.try<-function(){
   try_sp<-read.csv("clean_data/try_spp_clean.txt", header=FALSE, as.is=TRUE)[,1]
    read_csv("raw_data/tpl_names.txt")%>%
    filter(status=="Accepted")->acc_names
    correct.names <-tolower(gsub("_", " ", unique(acc_names$gs)))
   b<-try_sp[try_sp%in% correct.names]
   return(b)
   }

do_overlap_analysis<-function(){
  read_csv("raw_data/tpl_names.txt")%>%
  filter(status=="Accepted") ->lookup
    tpl_names <- tolower(unique(sub("_"," ",lookup$gs)))
  gbif<-get_gbif_names()
  try_names<-read.in.try()
  y<-list(
    genbank=get_genbank()[get_genbank()%in%tpl_names],
    try.all.names=try_names[try_names%in%tpl_names],
    gbif=gbif[gbif%in%tpl_names]
  )
  return(y)
}

write_overlap_table<-function(y){
  Number_of_overlapping_species<-c(sapply(y,length),calculate.overlap(y,2),calculate.overlap(y,3))
  write.csv(Number_of_overlapping_species,"tables/two_and_three_way_comparisons.csv")
  print(xtable(data.frame(Number_of_overlapping_species),caption="this is a caption"),file="tables/two_and_three_way_comparisons.tex",booktabs=TRUE,floating=FALSE,caption.placement="top", label = "tab:two_and_three_way_comparisons")
}




#FOR FUTURE CHECKS:
#tree.names[!tree.names%in%genbank]
#good.names<-sapply(as.character(syn$V1),FUN=function(x) strsplit(x,",")[[1]][1],USE.NAMES=F)

