

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
    
    names(y)[names(y)=="try.all.names"]<-"try-database"
    
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
  
  names(y)[names(y)=="try.all.names"]<-"try"
  
  lt$not_known <-  !lt$species%in%y$genbank & !lt$species%in%y$try
  ranking<-summarize(group_by(lt,family),sr=length(species),prop.sampled=mean(not_known),group=group[1])
  ranking$prop.sampled <- 1 - ranking$prop.sampled
  
  pdf("figures/treemap_not_known.pdf")
  treemap(ranking,
          index=c("group", "family"),
          vSize="sr",
          vColor="prop.sampled",
          type="manual",
          title="Taxonomic distribution of sampling for either genetic or functional data",
          palette=brewer.pal(11,"RdBu"),
          range=c(0,1))
  dev.off()
  
}
read.in.try<-function(){
t_try<-read_csv("TryAccSpecies.txt",col_names="AccSpeciesName")
return(tolower(t_try$AccSpeciesName))
}

makephylo<-function(){
  t_try<-read_csv("TryAccSpecies.txt",col_names="AccSpeciesName")
  try_out<-sync.species.lists(firstup(t_try$AccSpeciesName))
  tree<-readRDS("phy/phy.o.rds")
  group_by(try_out,order)%>%
    summarise(prop.sampled=mean(in.list))%>%
    filter(!is.na(order))->try.sampled.df
  
  t_gen<-get_genbank()
  
  gen_out<-sync.species.lists(firstup(t_gen))
  group_by(gen_out,order)%>%
    summarise(prop.sampled=mean(in.list))%>%
    filter(!is.na(order))->gen.sampled.df
  
gbif<-unique(get_gbif()))
  gbif_out<-sync.species.lists(firstup(gbif))
  group_by(gbif_out,order)%>%
    summarise(prop.sampled=mean(in.list))%>%
    filter(!is.na(order))->gbif.sampled.df


well_studied<-gbif_out[gbif%in%t_gen&gbif%in%t_try$AccSpeciesName]
well_studied_out<-sync.species.lists(firstup(well_studied))
  group_by(well_studied_out,order)%>%
    summarise(prop.sampled=mean(in.list))%>%
    filter(!is.na(order))->wellstudied.sampled.df

  phy.df<-data.frame(try=try.sampled.df$prop.sampled,
                      gen=gen.sampled.df$prop.sampled,
                      gbif=gbif.sampled.df$prop.sampled,
                      well_studied=wellstudied.sampled.df$prop.sampled)
  row.names(phy.df) <- try.sampled.df$order
  



  phy.o.dt <-
    ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% try.sampled.df$order])
  po <-
    ggtree(phy.o.dt) + geom_tiplab(size =1.5)
  #po<-ggtree(phy.o.dt)+geom_tiplab()
  cols = sample(rainbow(10), 40, replace = T)
  

  

  
  gheatmap(po, phy.df, offset = 2.5, width = 0.2)  + viridis::scale_fill_viridis(option ="D")
  ggsave("figures/order_phy.pdf")
}