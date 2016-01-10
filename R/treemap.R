

make_treemaps<-function(y){
  
  lookup <- read_csv("../../../srv/scratch/z3484779/taxonomicResources//plantList11syns.csv")
  out <- sub("_"," ",unique(lookup$correct.names))
  lt<-lookup_table(out,by_species = TRUE)
  lt$species<-row.names(lt)
  
  for (i in 1:3){
    lt$prop.sampled<-lt$species%in%y[[i]]
    ranking<-summarize(group_by(lt,family),sr=length(species),prop.sampled=mean(prop.sampled),group=group[1])
    #ranking$group<-plant_lookup()$group[match(ranking$family,plant_lookup()$family)]
    #ranking$order<-plant_lookup()$order[match(ranking$family,plant_lookup()$family)]
    
    #head(ranking)
    
    pdf(sprintf("figures/treemap_%s.pdf",names(y[i])))
    treemap(ranking,
            index=c("group", "family"),
            vSize="sr",
            vColor="prop.sampled",
            type="manual",
            title=sprintf("%s",names(y[i])),
            palette=brewer.pal(11,"RdBu"),
            range=c(0,1))
    dev.off()
  }
}

make_well_known_treemap<-function(y){

  lookup <- read_csv("../../../srv/scratch/z3484779/taxonomicResources//plantList11syns.csv")
  out <- sub("_"," ",unique(lookup$correct.names))
  lt<-lookup_table(out,by_species = TRUE)
  lt$species<-row.names(lt)
  
  lt$well_known <- lt$species%in%y$gbif & lt$species%in%y$zae & lt$species%in%y$diaz
  ranking<-summarize(group_by(lt,family),sr=length(species),prop.sampled=mean(well_known),group=group[1])
  
  r2<-subset(ranking,prop.sampled<0.2)
  
  pdf("figures/treemap_well_known.pdf")
  treemap(r2,
          index=c("group", "family"),
          vSize="sr",
          vColor="prop.sampled",
          type="manual",
          title="Taxonomic distribution of well-characterized species",
          palette=brewer.pal(11,"RdBu"),
          range=c(0,0.2))
  dev.off()
}
