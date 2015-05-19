
get_good_names<-function(){
  require(dplyr)
  syn<-read.delim("taxonomicResources//TPL1.1_synonymy_list",sep = "\t",header = FALSE,as.is=T,nrows=350699)
  good.names<-sapply(as.character(syn$V1),FUN=function(x) strsplit(x,",")[[1]][1],USE.NAMES=F)
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
  syn<-read.delim("taxonomicResources//TPL1.1_synonymy_list",sep = "\t",header = FALSE,as.is=T,nrows=350699)
  good.names<-sapply(as.character(syn$V1),FUN=function(x) strsplit(x,",")[[1]][1],USE.NAMES=F)
  good.names<-sub("_"," ",good.names)
  out<-filter(df,sp%in%good.names)
  return(out)
}

use.synonym.lookup<-function(orig.names){
  orig.names<-sub(" ","_",orig.names)
  #names(lookup)<-c("correct.names","all.names")
  #write.csv(lookup,"taxonomicResources//plantList11syns.csv",row.names=FALSE,quote=F)
  lookup<-read.csv("taxonomicResources//plantList11syns.csv",as.is=T)
  out.sp<-orig.names
  matched.spp<-lookup$correct.names[match(orig.names,lookup$all.names)]
  matched.spp<-sub("_"," ",matched.spp)
  orig.names<-sub("_"," ",orig.names)
  out.sp<-sub("_"," ",out.sp)
  out.sp[!is.na(matched.spp)]<-matched.spp[!is.na(matched.spp)]
  df<-data.frame(new.sp=matched.spp[!is.na(matched.spp)],old.name=orig.names[!is.na(matched.spp)])
  df<-subset(df,as.character(df$new.sp)!=as.character(df$old.name))
  write.csv(df,"output/syn_fixes_done.csv",row.names=F)
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

decode.nouranges<-function(N){
  support<-N$support
  support[support=="Aquatic herb"]<-"A"
  support[support=="Emergent tree"]<-"F"
  support[support=="Epilithic"]<-"F"
  support[support=="Epilithic herb"]<-"F"
  support[support=="Epilithic treelet"]<-"F"
  support[support=="Epiphyte"]<-"E"
  support[support=="Epiphyte & epilithic"]<-"H"
  support[support=="Epiphyte & terrestrial"]<-"H"
  support[support=="Epiphyte & tree"]<-"H"
  support[support=="Epiphyte herb"]<-"E"
  support[support=="Epiphytic"]<-"E"
  support[support=="Epiphytic & epilithic"]<-"E"
  support[support=="Epiphytic herb"]<-"E"
  support[support=="Hemipiphytic tree"]<-"H"
  support[support=="Herb"]<-"F"
  support[support=="Liana"]<-"C"
  support[support=="Medium tree"]<-"F"
  support[support=="Parasite"]<-"P"
  support[support=="Shrub"]<-"F"
  support[support=="Tree"]<-"F"
  support[support=="Vine"]<-"C"
  support[support=="tree"]<-"F"
  out<-data.frame(sp=N$sp,support=support,source="Bongers2004")
  out<-filter(out,support%in%c("C","E","H","P","A","F"))
  return(out)
}

decode.GF_Cornwell<-function(W){
  support<-W$support
  support[support=="tree"]<-"F"
  support[support=="subshrub"]<-"F"
  support[support=="shrub/tree"]<-"F"
  support[support=="shrub"]<-"F"
  support[support=="palm"]<-"F"
  support[support=="mistletoe"]<-"P"
  support[support=="liana/tree"]<-"C"
  support[support=="liana/shrub"]<-"C"
  support[support=="liana"]<-"C"
  support[support=="herb/subshrub"]<-"F"
  support[support=="herb"]<-"F"
  support[support=="graminoid"]<-"F"
  support[support=="fern"]<-"F"
  support[support=="cycad"]<-"F"
  support[support=="climber/shrub"]<-"C"
  support[support=="climber/liana/shrub"]<-"C"
  support[support=="climber/herb"]<-"C"
  support[support=="climber"]<-"C"
  out<-data.frame(sp=W$sp,support=support,source="Gallagher2015")
  out<-filter(out,support%in%c("C","P","F"))
  return(out)
}


loadData<-function(){
  # read in the raw data files
  apg<-read.csv("rawData/apgGF_data.csv",as.is=T)
  kew<-read.csv("rawData/kewLifeform_data.csv",as.is=T)
  decoder<-read.csv("rawData/gf_decoder.csv",as.is=T)
  TRY<-read.csv("rawData/TRY_Categorical_Traits_Lookup_Table_2012_03_17_TestRelease.csv",as.is=T)
  Couvreur<-read.csv("rawData/Couvreur2015Data.csv",as.is=T)
  Parthasararthy<-read.csv("rawData/Parthasararthy2004Data.csv",as.is=T)
  Mascaro<-read.csv("rawData/Mascaro2004Data.csv",as.is=T)
  Senbeta<-read.csv("rawData/Senbeta2005Data.csv",as.is=T)
  Dewalt<-read.csv("rawData/Dewalt2000Data.csv",as.is=T)
  Killeen<-read.csv("rawData/Killeen1998Data.csv",as.is=T)
  NabeNielsen<-read.csv("rawData/NabeNielsen2001Data.csv",as.is=T)
  Santos<-read.csv("rawData/Santos2009Data.csv",as.is=T)
  Uddin<-read.csv("rawData/Uddin2010Data.csv",as.is=T)
  Muthumperumal<-read.csv("rawData/Muthumperumal2009Data.csv",as.is=T)
  leda<-read.csv("rawData/leda.csv",as.is=T)
  strangler_fig<-read.csv("rawData/stangler_fig_list.csv")
  wiki.trees.shrubs<-load.wikipedia.trees("rawData//wikipedia_tree_shrub_list.txt")
  eMonocot<-read.csv("rawData/eMonocotData.csv",as.is=T)
  Cucurbitaceae<-read.csv("rawData/CucurbitaceaeStevens2001.csv",as.is=T)
  soares<-read.csv("rawData/soares2013.csv",as.is=T)
  stri.raw<-read.csv("rawData/stri_plant_species.txt",as.is=T)
  Ponnuchamy<-read.csv("rawData/Ponnuchamy2013Data.csv",as.is=T)
  Muoghalu<-read.csv("rawData/Muoghalu2005Data.csv",as.is=T)
  AddoFordjour<-read.csv("rawData/AddoFordjour2008Data.csv",as.is=T)
  Kelly<-read.csv("rawData/Kelly1985Data.csv",as.is=T)
  AfricanEpiphytes<-read.csv("rawData/AfricanepiphytesData.csv",as.is=T)
  bolivia<-read.csv("rawData//boliviaFlora.csv")
  aquaria<-read.csv("rawData//aquariumPlants.csv")
  cornwell<-read.csv("rawData//cornwellJR_data.csv")
  AddoFordjourRahmad<-read.csv("rawData/AddoFordjourRahmad2014Data.csv",as.is=T)
  Durgion<-read.csv("rawData/Durgion2011Data.csv",as.is=T)
  ElGhani<-read.csv("rawData/ElGhani2011Data.csv",as.is=T)
  FloraIsrael<-read.csv("rawData/FloraofIsraelOnlineData.csv",as.is=T)
  grassBase<-read.csv("rawData/grassBase.csv",as.is=T)
  Chettri<-read.csv("rawData/Chettri2010.csv",as.is=T)
  ClimbersEXtreeNA<-read.csv("rawData/ClimbersEXtreeNA.csv", as.is=T)
  Ding<-read.csv("rawData/Ding2014.csv", as.is=T)
  Githae<-read.csv("rawData/Githae2007.csv", as.is=T)
  Han<-read.csv("rawData/Han2010.csv", as.is=T)
  Hearn<-read.csv("rawData/Hearn2009.csv", as.is=T)
  InsideWood<-read.csv("rawData/InsideWoodData.csv", as.is=T)
  Krings<-read.csv("rawData/Krings2000.csv", as.is=T)
  Lu<-read.csv("rawData/Lu2008.csv", as.is=T)
  MolinaFreaner<-read.csv("rawData/MolinaFreaner1997.csv", as.is=T)
  Reddy<-read.csv("rawData/Reddy2003.csv", as.is=T)
  Solorzano<-read.csv("rawData/Solorzano2002.csv", as.is=T)
  Yuan<-read.csv("rawData/Yuan2009.csv", as.is=T)
  POSA<-read.csv("rawData/POSA.csv", as.is=T)
  nzflora<-read.csv("rawData/nzflora.csv", as.is=T)
  FloraMyanmar<-read.csv("rawData/FloraMyanmar.csv", as.is=T)
  AVH<-read.csv("rawData/AVHData.csv", as.is=T)
  ALA<-read.csv("rawData/AOLAData.csv", as.is=T)
  FloraMadagascar<-read.csv("rawData/MadagascarData.csv", as.is=T)
  SthSthAmericanEpi<-read.csv("rawData/Southsouthamericanepiphytes.csv", as.is=T)
  WstSthAmericanEpi<-read.csv("rawData/Westsouthamerica.csv", as.is=T)
  NthAmericanEpi<-read.csv("rawData/Northamericanepiphytes.csv", as.is=T)
  CrbEpi<-read.csv("rawData/Carribeanepiphytes.csv", as.is=T)
  AfrEpi<-read.csv("rawData/Africanepiphytes.csv", as.is=T)
  AsianTropiEpi<-read.csv("rawData/Asiatropical.csv", as.is=T)
  CtrlAmericanEpi<-read.csv("rawData/Centralamerica.csv", as.is=T)
  OzEpi<-read.csv("rawData/Australasianepiphytes.csv", as.is=T)
  BrazilEpi<-read.csv("rawData/Brazilepiphytes.csv", as.is=T)
  TempAsiaEpi<-read.csv("rawData/Temperateasia.csv", as.is=T)
  EuroEpi<-read.csv("rawData/Europeanepiphytes.csv", as.is=T)
  PacEpi<-read.csv("rawData/Pacificepiphytes.csv", as.is=T)
  FloraColumbia<-read.csv("rawData/FloraColumbia.csv", as.is=T)
  FloraParaguay<-read.csv("rawData/FloraParaguay.csv", as.is=T)
  FloraPeru<-read.csv("rawData/FloraPeru.csv", as.is=T)
  FloraEcuador<-read.csv("rawData/FloraEcuador.csv", as.is=T)
  FloraPalestina<-read.csv("rawData/FloraPalestina.csv", as.is=T)
  FloraPakistan<-read.csv("rawData/FloraPakistan.csv", as.is=T)
  eFloraChina<-read.csv("rawData/eFloraChina.csv", as.is=T)
  FloraPuertoRico<-read.csv("rawData/FloraPuertoRico.csv", as.is=T)
  NthSthAmericanEpi<-read.csv("rawData/Northsouthamerica.csv", as.is=T)
  AfrOrchids<-read.csv("rawData/Africanorchids.csv", as.is=T)
  SthAfrEpi<-read.csv("rawData/SouthernAfricanEpiphytes.csv", as.is=T)
  PeruEpi<-read.csv("rawData/Peruvianepiphytes.csv", as.is=T)
  EcuadorEpi<-read.csv("rawData/EcuadorEpiphytes.csv", as.is=T)
  EcuadorHemiEpi<-read.csv("rawData/EcuadorHemiepiphytes.csv", as.is=T)
  FloraWestAfrica<-read.csv("rawData/FloraWestAfrica.csv", as.is=T)
  Gentry<-read.csv("rawData/Gentry.csv", as.is=T)
  Geo<-read.csv("rawData/Geo.csv", as.is=T, colClasses=c(support="character"))
  Proteaceae<-read.csv("rawData/Proteaceae.csv", as.is=T, colClasses=c(support="character"))
  Lauraceae<-read.csv("rawData/Lauraceae.csv", as.is=T)
  Casuarinaceae<-read.csv("rawData/Casuarinaceae.csv", as.is=T, colClasses=c(support="character"))
  MadEpi<-read.csv("rawData/MadagascarEpi.csv", as.is=T)
  ParasiticPlantList<-read.csv("rawData/ParasiticPlantList.csv", as.is=T)
  AraceaeEpi<-read.csv("rawData/Araceae Epiphytes.csv", as.is=T)
  HemiEpiMisc<-read.csv("rawData/Hemiepiphytesmisc.csv", as.is=T)
  WstAfricanEpi<-read.csv("rawData/WestAfricanEpiphytes.csv", as.is=T)
  Nouranges<-read.csv("rawData/nouranges.csv", as.is=T)
  Goodeniaceae<-read.csv("rawData/Goodeniaceae.csv", as.is=T, colClasses=c(support="character"))
  Azioaceae<-read.csv("rawData/Azioaceae.csv", as.is=T, colClasses=c(support="character"))
  HelioEpi<-read.csv("rawData/HeliocereusEpi.csv", as.is=T)
  PeruCactiFree<-read.csv("rawData/PeruCactiFree.csv", as.is=T, colClasses=c(support="character"))
  EcuadorCactiFree<-read.csv("rawData/EcuadorCactiFree.csv", as.is=T, colClasses=c(support="character"))
  UndersampledFamilies<-read.csv("rawData/UndersampledFamilies.csv", as.is=T)
  Asteraceae<-read.csv("rawData/AsteraceaeFreeStanding.csv", as.is=T, colClasses=c(support="character"))
  Cactaceae<-read.csv("rawData/FreeStandCactaceae.csv", as.is=T, colClasses=c(support="character"))
  Mazaceae<-read.csv("rawData/Mazaceae.csv", as.is=T, colClasses=c(support="character"))
  eMonocot<-read.csv("rawData/eMonocot.csv", as.is=T)
  Plumbaginaceae<-read.csv("rawData/Plumbaginaceae.csv", as.is=T, colClasses=c(support="character"))
  UnderSampledFamilies2<-read.csv("rawData/UnderSampledFamilies2.csv", as.is=T)
  Balsaminaceae<-read.csv("rawData/Balsaminaceae.csv", as.is=T, colClasses=c(support="character"))
  refloraBrasil<-read.csv("rawData/refloraBrasil.csv", as.is=T)
  FranksAdditions<-read.csv("rawData/FranksAdditions.csv", as.is=T)
  AusHabit<-read.csv("rawData/GF_Cornwell.csv",as.is=T)
  
  # do intermediate data cleaning
  kew.decoded<-decode.kew(kew,decoder,data.source="kew")
  apg.ss<-subset.apg(apg)
  TRY.climbers<-find.try.climbers(TRY)
  TRY.epiphytes<-find.try.epiphytes(TRY)
  TRY.freestanding<-find.try.freestanding(TRY)
  leda<-decode.leda(leda,gw.decoder=decoder,data.source="leda")
  stri<-stri.cleaning(stri.raw)
  Nouranges<-decode.nouranges(Nouranges)
  AusHabit<-decode.GF_Cornwell(AusHabit)
  
  #stick everything together into one dataFrame
  
  all.data<-rbind(apg.ss,kew.decoded,TRY.climbers,TRY.epiphytes,Couvreur,Parthasararthy,Mascaro,Senbeta,Dewalt,Killeen,NabeNielsen,Uddin,Santos,Muthumperumal,leda,strangler_fig,wiki.trees.shrubs,eMonocot,TRY.freestanding,Cucurbitaceae,soares,stri,Ponnuchamy,Muoghalu,AddoFordjour,Kelly,bolivia,aquaria,AddoFordjourRahmad,Durgion,ElGhani,FloraIsrael,cornwell,AfricanEpiphytes,grassBase,Chettri,ClimbersEXtreeNA,Ding,Githae,Han,Hearn,InsideWood,Krings,Lu,MolinaFreaner,Reddy,Solorzano,Yuan,POSA,nzflora,FloraMyanmar,AVH,ALA,FloraMadagascar,SthSthAmericanEpi,WstSthAmericanEpi,NthAmericanEpi,CrbEpi,AfrEpi,AsianTropiEpi,CtrlAmericanEpi,OzEpi,BrazilEpi,TempAsiaEpi,EuroEpi,PacEpi,FloraColumbia,FloraParaguay,FloraPeru,FloraEcuador,FloraPalestina,FloraPakistan,eFloraChina,FloraPuertoRico,NthSthAmericanEpi,AfrOrchids,SthAfrEpi,PeruEpi,EcuadorEpi,EcuadorHemiEpi,FloraWestAfrica,Gentry,Geo,Proteaceae,Lauraceae,Casuarinaceae,MadEpi,ParasiticPlantList,AraceaeEpi,HemiEpiMisc,WstAfricanEpi,Nouranges,Goodeniaceae,Azioaceae,HelioEpi,PeruCactiFree,EcuadorCactiFree,UndersampledFamilies,Asteraceae,Cactaceae,Mazaceae,eMonocot,UnderSampledFamilies2,Balsaminaceae,refloraBrasil,Plumbaginaceae,FranksAdditions,AusHabit)
  return(all.data)
}

cleanData<-function(all.data=all.data){
  #fix the known errors
  all.data.fixed.support.errors<-fix.errors(all.data)
  
  #try to correct species names
  potential.fixes<-scrub(all.data.fixed.support.errors$sp)
  #try to save names
  #all.data.fixed.support.errors$sp[is.na(potential.fixes)]
  all.data.fixed.support.errors$sp<-potential.fixes
  all.data.fixed<-filter(all.data.fixed.support.errors,!is.na(sp))
  all.data.fixed<-fix.errors(all.data.fixed)
  
  #find conflicts
  all.data.fixed<-discard.internal.inconsistencies(all.data.fixed)
  conflicts<-find.conflicts(all.data.fixed)
  write.csv(conflicts,"output/conflicts.csv")
  
  
  #fix known spelling mistakes
  all.data.fixed$sp<-use.spelling.lookup(all.data.fixed$sp)
  #fix synonyms
  all.data.fixed$sp<-use.synonym.lookup(all.data.fixed$sp)
  
  #discard internal.inconsistencies
  all.data.fixed<-discard.internal.inconsistencies(all.data.fixed)
}

outputRawData<-function(data.out){
  write.csv(data.out,"output/growthFormData.csv",row.names=FALSE)
}

outputUniqueSpeciesData<-function(data.out){
  write.csv(data.out,"output/uniqueSpeciesGFData.csv",row.names=FALSE)
}
