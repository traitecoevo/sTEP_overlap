

make_treemap<-function(){

ranking$group<-plant_lookup()$group[match(ranking$family,plant_lookup()$family)]
ranking$order<-plant_lookup()$order[match(ranking$family,plant_lookup()$family)]

head(ranking)

pdf("sampling_in_try_by_family.pdf")
treemap(ranking,
        index=c("group", "family"),
        vSize="sr",
        vColor="prop.sampled",
        type="manual",
        palette=brewer.pal(9,"Greens"),
        range=c(0,1))
dev.off()
}
