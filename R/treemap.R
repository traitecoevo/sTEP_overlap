

make_treemaps<-function(y){
  
  lookup <- read_csv("../../../srv/scratch/z3484779/taxonomicResources//plantList11syns.csv")
  out <- sub("_"," ",unique(lookup$correct.names))
  lt<-lookup_table(out,by_species = TRUE)
  lt$species<-row.names(lt)
  
  y.smaller[1]
  
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
            title=sprintf(%s,names(y[i]))
            palette=brewer.pal(9,"RdGy"),
            range=c(0,1))
    dev.off()
  }
}