###You need to decide how you would like to graph your data
##For example, do you need to square root transform it?
##Is your data rarefied?
##Make sure you do not have any taxonomy data in this file
library(vegan)

OTUMatrix_Original <- read.csv("bact_alldata_taxatable.csv", header = TRUE)
str(OTUMatrix_Original)

#Remove first column because it is OTU names
test2 <- OTUMatrix_Original[,-1]
str(test2)

##This could beging with the OTU Matrix from Vegan
#Create a squareroot transformed matrix
sqrt_OTUmatrix_bact <-  sqrt(test2)

#now transpose the data so that Samples are Rows not Columns
OTU_Matrix<-t(sqrt_OTUmatrix_bact)
#str(OTU_Matrix)
#create the NMS object
Bact_NMS_data<-metaMDS(OTU_Matrix,distance = "bray", k=2,try=100,autotransform = TRUE,maxit=1000)
#data are square root transformed, and then rescaled ("Wisoncin rescaling")

#create vectors with the NMS attributes
NMS_coordinates<-scores(Bact_NMS_data,display="sites")
Bact_NMS_axes<-as.data.frame(NMS_coordinates)
write.csv(Bact_NMS_axes,"Bact_NMS_axes.csv")
NMS_OTUscores<-scores(Bact_NMS_data,display="species")

##Load Mapping File
Mapping_File <- read.csv("bact_alldata_mapfile.csv")
str(Mapping_File)
attach(Mapping_File)
vertical = as.factor(Vertposition)

#create dataframe with NMS Coordinates and Mapping File information
#add the proportion of total sequences of each vector to the "for_plotting" object below
for_ploting<-as.data.frame(cbind(NMS_coordinates,Mapping_File))
str(for_ploting)


#now plot these data
png(file="Bact_NMS_black_small.png", width = 6000, height = 4800, res = 1200)
par(mar=c(4,4,1,1))
plot(for_ploting$NMDS2 ~ for_ploting$NMDS1,
     xlab = "NMS1",
     ylab = "NMS2",
     font=2,
     font.lab=2,
     cex.axis=1,
     pch = c(0, 1, 2, 5, 4)[as.factor(Mapping_File$Vertposition)], cex=.8,  # different 'pch' types 
     data = for_ploting)
ordiellipse(Bact_NMS_data, group=Mapping_File$Vertposition,kind = "se", 
            conf=0.95, lwd=1.9, lty= c(1,5,2,4,3))
legend(
  x ="bottomright",
  legend = c("Litter","Understory","Subcanopy","Canopy","Emergent"), # for readability of legend
  pch = c(0, 1, 2,5,4),
  cex = .60 # scale the legend to look attractively sized
)

dev.off()

#plot the stress
stressplot(Bact_NMS_data)

#now these figure colorful!!! It might be useful to have colorful figures
png(file="Bact_NMS_color.png", width = 6000, height = 4800, res = 1200)
par(mar=c(4,4,1,1))
plot(for_ploting$NMDS2 ~ for_ploting$NMDS1,
     xlab = "NMS1",
     ylab = "NMS2",
     font=2,
     font.lab=2,
     cex.axis=1,
     pch = c(0, 1, 2, 5, 4)[as.factor(Mapping_File$Vertposition)], cex=.8, 
     col =c("black","saddlebrown","tan2","green1","green4")[as.factor(Mapping_File$Vertposition)],  # different 'pch' types 
     data = for_ploting)
ordiellipse(Bact_NMS_data, group=Mapping_File$Vertposition,kind = "se", 
            conf=0.95, lwd=1.9, col =c("black","saddlebrown","tan2","green1","green4"))
legend(
  x ="bottomright",
  legend = c("Litter","Understory","Subcanopy","Canopy","Emergent"), # for readability of legend
  pch = c(0, 1, 2,5,4),
  col =c("black","saddlebrown","tan2","green1","green4"),
  cex = .60 # scale the legend to look attractively sized
)

dev.off()

