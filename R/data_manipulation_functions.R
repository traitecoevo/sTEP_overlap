pacman::p_load(dplyr, raster,data.table,maptools,ape,readr,ggplot2,taxonlookup, snow, parallel,mgcv)

get_good_names<-function(){
  require(dplyr)
  syn<-read_csv("raw_data/tpl_names.txt")
  good.names<-syn$gs
  good.names<-sub("_"," ",good.names)
  return(unique(good.names))
}

get.traitfill.format<-function(plant.list.names,growthFormCode){
  gf<-cbind(split_genus(plant.list.names$sp),plant.list.names)
  gf.growth.form.count<-count.growth.form(gf,growthFormCode)
  merged_pl_gf<-merge(get_pl_species_by_genus_counts(),gf.growth.form.count,by="genus",all.x=T)
  merged_pl_gf$n1[is.na(merged_pl_gf$n1)]<-0
  merged_pl_gf$n0[is.na(merged_pl_gf$n0)]<-0
  merged_pl_gf$n0[merged_pl_gf$genus=="Angiopteris"]<-1
  gf_with_fam_order<-add.family.order(merged_pl_gf)
  gf_without_mosses<-filter(gf_with_fam_order,!is.na(family))
  gf_without_mosses$order[gf_without_mosses$order==""]<-paste0("Unplaced",gf_without_mosses$family[gf_without_mosses$order==""])
  return(gf_without_mosses)
}

run.traitfill<-function(plant.list.names){
  gf_without_mosses<-get.traitfill.format(plant.list.names,"C")  
  results<-traitfill(gf_without_mosses,500,with_replacement = T)
  results.F<-traitfill(gf_without_mosses,500,with_replacement = F)
  climbers<-rbind(results$overall,results.F$overall)
  climbers$gf<-"C"
  climbers$replace<-c("T","F")

  gf_without_mosses<-get.traitfill.format(plant.list.names,"E")  
  results<-traitfill(gf_without_mosses,500,with_replacement = T)
  results.F<-traitfill(gf_without_mosses,500,with_replacement = F)
  epiphytes<-rbind(results$overall,results.F$overall)
  epiphytes$gf<-"E"
  epiphytes$replace<-c("T","F")

  gf_without_mosses<-get.traitfill.format(plant.list.names,"A")  
  results<-traitfill(gf_without_mosses,500,with_replacement = T)
  results.F<-traitfill(gf_without_mosses,500,with_replacement = F)
  aquatic<-rbind(results$overall,results.F$overall)
  aquatic$gf<-"A"
  aquatic$replace<-c("T","F")
  
  write.csv(rbind(climbers,epiphytes,aquatic),"output/traitfillResults.csv")
}


get_pl_species_by_genus_counts<-function(){
  pl_names<-get_good_names()
  out<-split_genus(pl_names)
  num.of.species.by.genus<-summarize(group_by(out,genus),N=length(species))
  return(num.of.species.by.genus)
}

split_genus <- function(str) {
  x <- strsplit(str," +")
  n <- sapply(x, length)
  x[n > 2] <- lapply(x[n > 2], function(el) el[1:2])
  g <- do.call("rbind", x, quote=TRUE)
  g <- matrix(unlist(x), ncol=2, byrow=TRUE)
  colnames(g) <- c("genus", "species")
  as.data.frame(g)
}

add.family.order<-function(gf){
  require(TaxonLookup)
  lookUp<-plant_lookup()
  gf$family<-lookUp$family[match(gf$genus,lookUp$genus)]
  gf$order<-lookUp$order[match(gf$genus,lookUp$genus)]
  return(gf)
}



plant.list.subset<-function(df){
  require(dplyr)
  syn<-read.delim("raw_data/tpl_names.txt",sep = "\t",header = FALSE,as.is=T,nrows=350699)
  good.names<-sapply(as.character(syn$V1),FUN=function(x) strsplit(x,",")[[1]][1],USE.NAMES=F)
  good.names<-sub("_"," ",good.names)
  out<-filter(df,sp%in%good.names)
  return(out)
}

use.synonym.lookup<-function(orig.names){
  orig.names<-sub(" ","_",orig.names)
  #names(lookup)<-c("correct.names","all.names")
  #write.csv(lookup,"taxonomicResources//plantList11syns.csv",row.names=FALSE,quote=F)
  lookup<-read_csv("raw_data/tpl_names.txt")
  out.sp<-orig.names
  matched.spp<-lookup$correct.names[match(orig.names,lookup$all.names)]
  
  #change all the species in all.names
  out.sp[!is.na(matched.spp)]<-matched.spp[!is.na(matched.spp)]
  
  #change back all the species that are also accepted names (this is for a case where some subspecies have listed synonyms)
  out.sp[orig.names%in%lookup$correct.names]<-orig.names[orig.names%in%lookup$correct.names]
  
  out.sp<-sub("_"," ",out.sp)

  return(out.sp)
}



use.spelling.lookup<-function(orig.names){
  orig.names<-sub(" ","_",orig.names)
  lookup<-read.csv("taxonomicResources//spelling_fixes.csv",as.is=T)
  out.sp<-orig.names
  matched.spp<-lookup$good[match(orig.names,lookup$bad)]
  matched.spp<-sub("_"," ",matched.spp)
  orig.names<-sub("_"," ",orig.names)
  out.sp<-sub("_"," ",out.sp)
  out.sp[!is.na(matched.spp)]<-matched.spp[!is.na(matched.spp)]
  df<-data.frame(new.sp=matched.spp[!is.na(matched.spp)],old.name=orig.names[!is.na(matched.spp)])
  df<-subset(df,as.character(df$new.sp)!=as.character(df$old.name))
  write.csv(df,"output/spelling_fixes_done.csv",row.names=F)
  return(out.sp)
}



scrub <- function(sp.names) {
  
  ## ###################################################
  ## FILTER 1: Identify common mistakes and            #
  ## fix names in sp.names                             #
  ## ###################################################
  
  ##---------------------------------------------------#
  ## Turn all underscores and >1 space between names
  ## into single space
  ##---------------------------------------------------#
  re.underscore <- "[_ \n]+"
  fix.underscore <- grepl(re.underscore, sp.names)
  sp.names[fix.underscore] <-
    gsub(re.underscore," ", sp.names[fix.underscore])
  
  ##---------------------------------------------------#
  ## Delete cf collapsing to binomial
  ## Delete aff. or sp. aff., collapsing to binomial
  ##---------------------------------------------------#
  re.aff <- "(.*) (cf\\.?|aff\\.?|sp\\.? aff\\.?) (.*)"
  fix.aff <- grepl(re.aff, sp.names)
  sp.names[fix.aff] <- gsub(re.aff, "\\1 \\3", sp.names[fix.aff])
  
  re.cf2 <- "(.*)(\\.cf.*\\.)(.*)"
  fix.cf2 <- grepl(re.cf2, sp.names)
  sp.names[fix.cf2] <- gsub(re.cf2, "\\1 \\3", sp.names[fix.cf2])
  
  ##---------------------------------------------------#
  ## Remove everything following a var. subsp. ssp. subvar.
  ##---------------------------------------------------#
  rm.str1 <- paste("\\.*",c("var", "subsp", "ssp", "sp", "subvar",
                            "agg", "cv"), " ", sep="")
  rm.str2 <- paste("\\.*",c("var", "subsp", "ssp", "sp", "subvar",
                            "agg", "cv"), "\\.", sep="")
  rm.str <- c(rm.str1, rm.str2)
  re.ssp <- paste("^(.*) (", paste(rm.str, collapse="|"), ").*$", sep="")
  fix.ssp <- grepl(re.ssp, sp.names)
  sp.names[fix.ssp] <- sub(re.ssp, "\\1", sp.names[fix.ssp])
  
  ##---------------------------------------------------#
  ## Bionomials with alternate genus in parentheses
  ## Dump things in parentheses, collapsing to binomial
  ##---------------------------------------------------#
  
  ## ending in parentheses
  re.end.in.par <- " ?\\(.*\\)$"
  fix.end.in.par <- grepl(re.end.in.par, sp.names)
  sp.names[fix.end.in.par] <-
    sub(re.end.in.par, "", sp.names[grepl(re.end.in.par, sp.names)])
  
  ## parentheses embedded in name
  re.syn <- "(.*)(\\(.+\\))(.+)$"
  fix.syn <- grepl(re.syn, sp.names)
  sp.names[fix.syn] <- gsub(re.syn, "\\1 \\3", sp.names[fix.syn])
  
  ## Turn all spaces between names into single space
  sp.names <- gsub("( )+"," ", sp.names)
  
  ##---------------------------------------------------#
  ## Extract potential binomials, trinomials, etc.
  ## and collapse to binomails
  ##---------------------------------------------------#
  
  ## All binomials or trinomials
  re.bi.tri <- "(^[A-Z][a-z]+) ([-a-z]+)( .*)+$"
  fix.bi.tri <- grepl(re.bi.tri, sp.names)
  sp.names[fix.bi.tri] <-
    gsub(re.bi.tri, "\\1 \\2", sp.names[fix.bi.tri])
  
  ## Nuke trailing whitespace:
  sp.names <- sub(" +$", "", sp.names)
  ## Nuke trailing periods from possible species names (not
  ## abbreviations)
  sp.names <- sub("([a-z]{4,})\\.+$", "\\1", sp.names)
  
  ## ###################################################
  ## FILTER 2: Identify names to exclude               #
  ## ###################################################
  ## Logical vector for tracking which names are good binomials
  good.names <- rep(FALSE, length(sp.names))
  bad.names <- rep(FALSE, length(sp.names))
  
  re.potentials <- "^[A-Z][a-z]+(-[A-Z][a-z]+)? [-a-z]+( .*)*$"
  good.names[grep(re.potentials, sp.names)] <- TRUE
  
  ## Exclude species name "indet"
  ex.indet <- grepl("[A-Za-z]+ indet", sp.names)
  bad.names[ex.indet] <- TRUE
  
  ## Exclude species name "nondet"
  ex.nondet <- grepl("[A-Za-z]+ nondet", sp.names)
  bad.names[ex.nondet] <- TRUE
  
  ## Exclude species name "sect"
  ex.sect <- grepl("[A-Za-z]+ sect", sp.names)
  bad.names[ex.sect] <- TRUE
  
  ## Exclude species name "x"
  ex.hybrid <- grepl("[A-Za-z]+ x$", sp.names)
  bad.names[ex.hybrid] <- TRUE
  
  ## Exclude species that start with parentheses
  re.start.in.par <- "^\\((.*)$"
  fix.start.in.par <- grepl(re.start.in.par, sp.names)
  bad.names[fix.start.in.par] <- TRUE
  
  ## Exclude rows with species name that end in sp. or sp
  ## and followed by numbers
  re.sp <- "(.+) (sp.|sp|spp.|sp[0-9]+|sp.[0-9]+)( |$)"
  fix.sp <- grepl(re.sp, sp.names)
  bad.names[fix.sp] <- TRUE
  
  ## Exclude species names that include question marks
  re.quest <- "\\?"
  fix.quest <- grepl(re.quest, sp.names)
  bad.names[fix.quest] <- TRUE
  
  ## Exclude hybrids
  re.x <- " (x|X) "
  fix.x <- grepl(re.x, sp.names)
  bad.names[fix.x] <- TRUE
  
  ## Exclude mono names:
  bad.names[!grepl(" ", sp.names)] <- TRUE
  ## Exclude anything with a number
  bad.names[grepl("[0-9]", sp.names)] <- TRUE
  
  
  
  sort(unique(sp.names[!good.names & !bad.names]))
  
  sp.names[!(good.names & !bad.names)] <- NA
  
  
  return(sp.names)
}


