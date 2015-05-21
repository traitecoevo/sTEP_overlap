centers<-read.table("species_centers.txt")

pdf("temp.pdf")
  plot(x=centers$long.mean,y=centers$lat.mean,pch=16,col=rgb(0, 1, 0,0.1))  
dev.off()
