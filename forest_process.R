library(taxonlookup)
library(sp)
library(raster)

#a<-filter(a,!is.na(as.numeric(a$decimalLongitude)))
a<-fread("occur.csv")
sp<-SpatialPoints(cbind(as.numeric(a$decimalLongitude),as.numeric(a$decimalLatitude)))
forests<-raster("tropicalForestBins.grd")
forests.x<-extract(forests,sp)

america<-filter(a,forests.x==2)
#sort(table(america$species),decreasing = T)[1:10]

mutate(america,species=scrub(species))%>%
  mutate(species.fix=use.synonym.lookup(species))->america

america$genus<-taxonlookup:::split_genus(america$species.fix)

#sapply(as.character(america$species.fix),FUN=function(x) strsplit(x," ")[[1]][1],USE.NAMES=F)



america$family<-plant_lookup()$family[match(america$genus,plant_lookup()$genus)]

america.out<-data.frame(species=names(sort(table(america$species),decreasing = T)[1:20]),
                        species.count=sort(table(america$species),decreasing = T)[1:20],           
                        genus=names(sort(table(america$genus),decreasing = T)[1:20]),
                        genus.count=sort(table(america$genus),decreasing = T)[1:20],
                        family=names(sort(table(america$family),decreasing = T)[1:20]),
                        family.count=sort(table(america$family),decreasing = T)[1:20],stringsAsFactors=F)
write.csv(america.out,"high_elevation_dry_america.csv",row.names=F)


