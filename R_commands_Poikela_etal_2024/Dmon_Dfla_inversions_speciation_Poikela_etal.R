#installing libraries needed
require(ggplot2)
library(ggplot2)
library(reshape2)
library(ggridges)
library(gghalves)
library(scales)
library(dplyr)
library(ggforce)
library(Publish)
library(patchwork)
library(devtools)
library(ggbiplot)
library(rgdal)
library(dismo)
library(raster)
library(PupillometryR)
library(forcats)
library(tidyverse)
library(ggrepel)
library(FactoMineR)
library(factoextra)

#### GENOMIC PCA ####

pca <-read.delim("gimble2.montana_flavomontana.nopruning.eigenvec", header=T) #read in eigenvectors
eigenval <-read.delim("gimble2.montana_flavomontana.nopruning.eigenval", header=F) #read in eignevalues

# variances
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
pve

sum_eig<-sum(eigenval$V1)
sum_eigs<-lapply(eigenval$V1,function(x){
  rt<-(x/sum_eig)*100
  rt<-round(rt)
  return(rt)
})

sum_eigs

legend_title<-"Population"

pca$POP <- factor(pca$POP, levels = c("alaska_mon", "west_mon", "west_fla", "rocky_mon", "rocky_fla"))

pca <- pca %>%
  mutate(SITE = if_else(SITE == "Honolulu_Creek", "Honolulu Creek", SITE))
pca <- pca %>%
  mutate(SITE = if_else(SITE == "Fall_Creek", "Fall Creek", SITE))

a<-ggplot(pca , aes(PC1, PC2)) +
  labs(tag="A") +
  geom_point(aes(shape=POP, color=POP, fill=POP), size=4, alpha=0.6) +
  xlab(paste0("PC1: ",sum_eigs[[1]],"% variance")) +
  ylab(paste0("PC2: ",sum_eigs[[2]],"% variance")) +
  geom_text_repel(aes(label=SITE), size=3.5, max.overlaps=21) + #The part with geom_text_repel adds easy-readable labels.
  scale_color_manual(values=c("black", "black", "black", "black", "black"), labels = c(expression('Alaska'~italic("D. mon")~''),expression('western coast'~italic("D. mon")~''),expression('western coast'~italic("D. fla")~''),expression('Rocky mountains'~italic("D. mon")~''), expression('Rocky mountains'~italic("D. fla")~''))) +
  scale_fill_manual(values=c("#1a1aff", "#1a1aff", "#FFCC00", "#1a1aff", "#FFCC00"), labels = c(expression('Alaska'~italic("D. mon")~''),expression('western coast'~italic("D. mon")~''),expression('western coast'~italic("D. fla")~''),expression('Rocky mountains'~italic("D. mon")~''), expression('Rocky mountains'~italic("D. fla")~''))) +
  scale_shape_manual(values=c(22, 21, 21, 24, 24), labels = c(expression('Alaska'~italic("D. mon")~''),expression('western coast'~italic("D. mon")~''),expression('western coast'~italic("D. fla")~''),expression('Rocky mountains'~italic("D. mon")~''), expression('Rocky mountains'~italic("D. fla")~''))) +
  theme_minimal() +
  theme(legend.position="none", panel.background = element_rect(colour='grey')) +
  theme(text=element_text(size=11))
a

ggsave(a, filename = "Genomic_PCA.png", width = 3.5, height = 3.5)


#### Bioclimatic variables ####

# http://www.worldclim.org/formats1
# https://www.rdocumentation.org/packages/dismo/versions/1.1-4/topics/biovars
# http://www.worldclim.org/formats1

# file size ~40Mb
bio1 <- raster(x="wc2.1_2.5m_bio_1.tif")
bio2 <- raster(x="wc2.1_2.5m_bio_2.tif")
bio3 <- raster(x="wc2.1_2.5m_bio_3.tif")
bio4 <- raster(x="wc2.1_2.5m_bio_4.tif")
bio5 <- raster(x="wc2.1_2.5m_bio_5.tif")
bio6 <- raster(x="wc2.1_2.5m_bio_6.tif")
bio7 <- raster(x="wc2.1_2.5m_bio_7.tif")
bio8 <- raster(x="wc2.1_2.5m_bio_8.tif")
bio9 <- raster(x="wc2.1_2.5m_bio_9.tif")
bio10 <- raster(x="wc2.1_2.5m_bio_10.tif")
bio11 <- raster(x="wc2.1_2.5m_bio_11.tif")
bio12 <- raster(x="wc2.1_2.5m_bio_12.tif")
bio13 <- raster(x="wc2.1_2.5m_bio_13.tif")
bio14 <- raster(x="wc2.1_2.5m_bio_14.tif")
bio15 <- raster(x="wc2.1_2.5m_bio_15.tif")
bio16 <- raster(x="wc2.1_2.5m_bio_16.tif")
bio17 <- raster(x="wc2.1_2.5m_bio_17.tif")
bio18 <- raster(x="wc2.1_2.5m_bio_18.tif")
bio19 <- raster(x="wc2.1_2.5m_bio_19.tif")

coordinates <- read.table("coordinates.txt", header = T)
coordinates

x <- coordinates$Latitude
y <- coordinates$Longitude

bio1_results<-extract(bio1, SpatialPoints(cbind(y, x)))
bio2_results<-extract(bio2, SpatialPoints(cbind(y, x)))
bio3_results<-extract(bio3, SpatialPoints(cbind(y, x)))
bio4_results<-extract(bio4, SpatialPoints(cbind(y, x)))
bio5_results<-extract(bio5, SpatialPoints(cbind(y, x)))
bio6_results<-extract(bio6, SpatialPoints(cbind(y, x)))
bio7_results<-extract(bio7, SpatialPoints(cbind(y, x)))
bio8_results<-extract(bio8, SpatialPoints(cbind(y, x)))
bio9_results<-extract(bio9, SpatialPoints(cbind(y, x)))
bio10_results<-extract(bio10, SpatialPoints(cbind(y, x)))
bio11_results<-extract(bio11, SpatialPoints(cbind(y, x)))
bio12_results<-extract(bio12, SpatialPoints(cbind(y, x)))
bio13_results<-extract(bio13, SpatialPoints(cbind(y, x)))
bio14_results<-extract(bio14, SpatialPoints(cbind(y, x)))
bio15_results<-extract(bio15, SpatialPoints(cbind(y, x)))
bio16_results<-extract(bio16, SpatialPoints(cbind(y, x)))
bio17_results<-extract(bio17, SpatialPoints(cbind(y, x)))
bio18_results<-extract(bio18, SpatialPoints(cbind(y, x)))
bio19_results<-extract(bio19, SpatialPoints(cbind(y, x)))


write(bio1_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio2_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio3_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio4_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio5_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio6_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio7_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio8_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio9_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio10_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio11_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio12_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio13_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio14_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio15_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio16_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio17_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio18_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")
write(bio19_results, file = "bioclimatic_data_2.5min.txt", ncolumns=16,
      append = TRUE, sep = "\t")



#### Climate PCA ####

# read df
ccrt <- read.table(("bioclimatic_variables_14populations.txt"), fill = TRUE, header=T, dec=".")

ccrt2 <- ccrt[,-1]
rownames(ccrt2) <- ccrt[,1]

rownames(ccrt2)

rownames(ccrt2) <- c("Fairbanks", "Honolulu Creek", "Seward", "Terrace", "Vancouver", "Ashford",
                     "Fall Creek", "Azalea", "McBride", "Cranbrook", "Livingston", "Jackson",
                     "Afton", "Liberty")

res.pca <- PCA(ccrt2[,4:22], graph = FALSE)

#  name               description                          
#1  "$eig"             "eigenvalues"                        
#2  "$var"             "results for the variables"          
#3  "$var$coord"       "coord. for the variables"           
#4  "$var$cor"         "correlations variables - dimensions"
#5  "$var$cos2"        "cos2 for the variables"             
#6  "$var$contrib"     "contributions of the variables"     
#7  "$ind"             "results for the individuals"        
#8  "$ind$coord"       "coord. for the individuals"         
#9  "$ind$cos2"        "cos2 for the individuals"           
#10 "$ind$contrib"     "contributions of the individuals"   
#11 "$call"            "summary statistics"                 
#12 "$call$centre"     "mean of the variables"              
#13 "$call$ecart.type" "standard error of the variables"    
#14 "$call$row.w"      "weights for the individuals"        
#15 "$call$col.w"      "weights for the variables"  

res.pca$ind
contr<-res.pca$var$contrib

eig.val <- res.pca$eig
barplot(eig.val[, 2], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Explained by PCs (%)",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.val), eig.val[, 2], 
      type = "b", pch = 19, col = "red")


plot(res.pca, choix = "ind", autoLab = "yes")

plot(res.pca, choix = "var", autoLab = "yes")

a<-fviz_pca_biplot(res.pca,
                   col.var = "#FF99FF",
                   col.ind = rownames(ccrt2),
                   repel = TRUE, title = "",
                   labelsize = 5,
                   ylab="PC2: 29% variance", xlab="PC1: 51% variance",
                   pointsize = 2) +
  labs(tag="B") +
  scale_shape_manual(values=c(17,16,16,17,15,16,15,17,17,17,17,15,16,16)) +
  scale_color_manual(values=c("#663333","#6600FF","#6600FF","#663333","#339999",
                              "#6600FF","#339999","#663333","#663333","#663333",
                              "#663333","#339999","#6600FF","#6600FF")) +
  scale_y_reverse() +
  scale_x_reverse() +
  theme(legend.position="none", text=element_text(size=14),
        axis.text=element_text(size=14), axis.line = element_line(color = 'grey'),
        panel.background = element_rect(fill = "white", colour = "grey"))
a

ggsave(a, filename = "Climate_PCA.png", width = 5, height = 5)



#### Mean genetic divergence (Dxy) ####

d_xy <- read.table("Dxy_intergenic_intronic_coding_allopatric_sympatric.64bp.txt", header=T)


d_xy$Geographic_genetic_comparison <- factor(d_xy$Geographic_genetic_comparison, levels = c("Sympatry_intergenic",
                                                                                            "Allopatry_intergenic",
                                                                                            "Sympatry_intronic",
                                                                                            "Allopatry_intronic",
                                                                                            "Sympatry_coding",
                                                                                            "Allopatry_coding"))
levels(d_xy$Geographic_genetic_comparison)[levels(d_xy$Geographic_genetic_comparison)=="Allopatry_intergenic"] <- "Allopatry, intergenic"
levels(d_xy$Geographic_genetic_comparison)[levels(d_xy$Geographic_genetic_comparison)=="Sympatry_intergenic"] <- "Sympatry, intergenic"
levels(d_xy$Geographic_genetic_comparison)[levels(d_xy$Geographic_genetic_comparison)=="Allopatry_intronic"] <- "Allopatry, intronic"
levels(d_xy$Geographic_genetic_comparison)[levels(d_xy$Geographic_genetic_comparison)=="Sympatry_intronic"] <- "Sympatry, intronic"
levels(d_xy$Geographic_genetic_comparison)[levels(d_xy$Geographic_genetic_comparison)=="Allopatry_coding"] <- "Allopatry, coding"
levels(d_xy$Geographic_genetic_comparison)[levels(d_xy$Geographic_genetic_comparison)=="Sympatry_coding"] <- "Sympatry, coding"

d_xy$Chromosomal_region <- factor(d_xy$Chromosomal_region, levels = c("colinear_autosomes",
                                                                      "4_inversion",
                                                                      "5_inversion",
                                                                      "colinear_X",
                                                                      "X_inversion"))
levels(d_xy$Chromosomal_region)[levels(d_xy$Chromosomal_region)=="colinear_autosomes"] <- "COL"
levels(d_xy$Chromosomal_region)[levels(d_xy$Chromosomal_region)=="4_inversion"] <- "4 INV"
levels(d_xy$Chromosomal_region)[levels(d_xy$Chromosomal_region)=="5_inversion"] <- "5 INV"
levels(d_xy$Chromosomal_region)[levels(d_xy$Chromosomal_region)=="colinear_X"] <- "COL "
levels(d_xy$Chromosomal_region)[levels(d_xy$Chromosomal_region)=="X_inversion"] <- "X INV"

d_xy$Chromosome <- factor(d_xy$Chromosome, levels = c("Autosomes", "X"))

text1<-data.frame(x=1.8,y=0.04, Chromosome="Autosomes")
text2<-data.frame(x=2.2,y=0.04, Chromosome="Autosomes")

text3<-data.frame(x=2.8,y=0.04, Chromosome="Autosomes")
text4<-data.frame(x=3.2,y=0.04, Chromosome="Autosomes")

text5<-data.frame(x=1.8,y=0.04, Chromosome="X")
text6<-data.frame(x=2.2,y=0.04, Chromosome="X")

text7<-data.frame(x=1.8,y=0.015, Chromosome="Autosomes")
text8<-data.frame(x=2.2,y=0.015, Chromosome="Autosomes")

text9<-data.frame(x=1.8,y=0.015, Chromosome="X")
text10<-data.frame(x=2.2,y=0.015, Chromosome="X")

legend_title<-"Comparison"

a<-ggplot(d_xy, aes(y=DXY, x=Chromosomal_region)) + 
  labs(x="", y="Mean genetic divergence\n(dxy)") +
  geom_point(aes(color=Geographic_genetic_comparison, shape=Geographic_genetic_comparison, fill=Geographic_genetic_comparison), 
             alpha=1, size=3, position = position_dodge(width = 0.4)) +
  facet_grid(cols=vars(Chromosome), scales="free_x", space = "free") +
  scale_color_manual(legend_title, values = c("black", "black", "black", "black", "black", "black")) +
  scale_fill_manual(legend_title, values = c("#6CD1AC", "#7D0DDB", "#6CD1AC", "#7D0DDB", "#6CD1AC", "#7D0DDB")) +
  scale_shape_manual(legend_title, values = c(21, 21, 24, 24, 22, 22)) +
  scale_y_continuous(limits = c(0.015, 0.04),
                     breaks=c(0.02, 0.03, 0.04)) +
  geom_text(data = text1, aes(x=x, y=y), label="***", size = 5, color="#009966") +
  geom_text(data = text2, aes(x=x, y=y), label="***", size = 5, color="#7D0DDB") +
  geom_text(data = text3, aes(x=x, y=y), label="***", size = 5, color="#009966") +
  geom_text(data = text4, aes(x=x, y=y), label="***", size = 5, color="#7D0DDB") +
  geom_text(data = text5, aes(x=x, y=y), label="***", size = 5, color="#009966") +
  geom_text(data = text6, aes(x=x, y=y), label="***", size = 5, color="#7D0DDB") +
  geom_text(data = text7, aes(x=x, y=y), label="***", size = 5, color="#009966") +
  geom_text(data = text8, aes(x=x, y=y), label="***", size = 5, color="#7D0DDB") +
  geom_text(data = text9, aes(x=x, y=y), label="***", size = 5, color="#009966") +
  geom_text(data = text10, aes(x=x, y=y), label="***", size = 5, color="#7D0DDB") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(legend.text=element_text(size=9), 
        axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(size=9),
        legend.title = element_text(size=9),
        strip.text.x = element_text(size=8))
a

ggsave(a, filename = "C:/Users/Noora/Desktop/mean_Dxy.png", height = 3, width = 6)




#### Mean genetic diversity (Pi) ####

d_xy <- read.table("Pi_intergenic_intronic_coding_allopatric_sympatric.64bp.txt", header=T)

d_xy$Genetic_region <- factor(d_xy$Genetic_region, levels = c("Intergenic", "Intronic", "Coding"))

d_xy$Population_species <- factor(d_xy$Population_species, levels = c("mon_Allopatry", "mon_Sympatry", "fla_Sympatry"))

d_xy$Genetic_region_chromosomal_region <- factor(d_xy$Genetic_region_chromosomal_region, 
                                                 levels = c("Intergenic_autosomes", 
                                                            "Intergenic_X", 
                                                            "Intronic_autosomes",
                                                            "Intronic_X",
                                                            "Coding_autosomes",
                                                            "Coding_X"))

levels(d_xy$Genetic_region_chromosomal_region)[levels(d_xy$Genetic_region_chromosomal_region)=="Intergenic_autosomes"] <- "Autosomes,\nintergenic"
levels(d_xy$Genetic_region_chromosomal_region)[levels(d_xy$Genetic_region_chromosomal_region)=="Intergenic_X"] <- "X,\nintergenic"
levels(d_xy$Genetic_region_chromosomal_region)[levels(d_xy$Genetic_region_chromosomal_region)=="Intronic_autosomes"] <- "Autosomes,\nintrons"
levels(d_xy$Genetic_region_chromosomal_region)[levels(d_xy$Genetic_region_chromosomal_region)=="Intronic_X"] <- "X,\nintrons"
levels(d_xy$Genetic_region_chromosomal_region)[levels(d_xy$Genetic_region_chromosomal_region)=="Coding_autosomes"] <- "Autosomes,\ncoding"
levels(d_xy$Genetic_region_chromosomal_region)[levels(d_xy$Genetic_region_chromosomal_region)=="Coding_X"] <- "X,\ncoding "

d_xy$Chromosomal_region <- factor(d_xy$Chromosomal_region, levels = c("colinear_autosomes",                                                                 "2Lb_inversion",
                                                                      "4_inversion",
                                                                      "5_inversion",
                                                                      "colinear_X",
                                                                      "X_inversion"))
levels(d_xy$Chromosomal_region)[levels(d_xy$Chromosomal_region)=="colinear_autosomes"] <- "COL"
levels(d_xy$Chromosomal_region)[levels(d_xy$Chromosomal_region)=="4_inversion"] <- "4 INV"
levels(d_xy$Chromosomal_region)[levels(d_xy$Chromosomal_region)=="5_inversion"] <- "5 INV"
levels(d_xy$Chromosomal_region)[levels(d_xy$Chromosomal_region)=="colinear_X"] <- "COL "
levels(d_xy$Chromosomal_region)[levels(d_xy$Chromosomal_region)=="X_inversion"] <- "X INV"

a<-ggplot(d_xy, aes(y=Pi, x=Chromosomal_region)) + 
  labs(x="", y="Mean genetic diversity\n(π)", color="Population", fill="Population", shape="Population") +
  geom_point(aes(color=Population_species, fill=Population_species, shape=Population_species), alpha=1, size=2) + #position = position_dodge(width = 0.7)
  facet_grid(cols=vars(Genetic_region_chromosomal_region), scales="free_x", space = "free") +
  scale_color_manual(values = c("black", "black", "black"), labels = c(expression('Allopatry'~italic("D. mon")~''),
                                                                       expression('Sympatry'~italic("D. mon")~''),
                                                                       expression('Sympatry'~italic("D. fla")~''))) +
  scale_fill_manual(values = c("#9999FF", "#9999FF", "#FF9933"), labels = c(expression('Allopatry'~italic("D. mon")~''),
                                                                            expression('Sympatry'~italic("D. mon")~''),
                                                                            expression('Sympatry'~italic("D. fla")~''))) +
  scale_shape_manual(values = c(21, 23, 21), labels = c(expression('Allopatry'~italic("D. mon")~''),
                                                        expression('Sympatry'~italic("D. mon")~''),
                                                        expression('Sympatry'~italic("D. fla")~''))) +
  scale_y_continuous(limits = c(0, 0.014),
                     breaks=c(0, 0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(legend.text=element_text(size=8), 
        axis.text=element_text(size=8),
        axis.title = element_text(size=8),
        plot.title = element_text(size=8),
        legend.title = element_text(size=8),
        strip.text.x = element_text(size=8))
a

ggsave(a, filename = "C:/Users/Noora/Desktop/mean_Pi.png", height = 2.5, width = 7)


#### Genome-wide dxy plot ####

# Read the data frames
dxy1 <- read.table("2montana.flavomontana.allopatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header = TRUE)
dxy2 <- read.table("2montana.flavomontana.sympatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header = TRUE)
dxy3 <- read.table("2montana.allopatric.sympatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header = TRUE)

# Scale the start column
dxy1$start_NEW <- dxy1$start / 1000000
dxy2$start_NEW <- dxy2$start / 1000000
dxy3$start_NEW <- dxy3$start / 1000000

# Replace sequence values
dxy_replace_sequence <- function(df) {
  df$sequence[df$sequence == 'X_chromosome'] <- 'X'
  df$sequence[df$sequence == '2L_chromosome'] <- '2L'
  df$sequence[df$sequence == '2R_chromosome'] <- '2R'
  df$sequence[df$sequence == '3_chromosome'] <- '3'
  df$sequence[df$sequence == '4_chromosome'] <- '4'
  df$sequence[df$sequence == '5_chromosome'] <- '5'
  df$sequence <- factor(df$sequence, levels = c("X", "2L", "2R", "3", "4", "5"))
  df$sequence_NEW <- factor(df$sequence, levels = c("X", "2L", "2R", "3", "4", "5"), labels = c("X", "2L", "2R", "3", "4", "5"))
  return(df)
}

dxy1 <- dxy_replace_sequence(dxy1)
dxy2 <- dxy_replace_sequence(dxy2)
dxy3 <- dxy_replace_sequence(dxy3)

# Create segment data frame
data.segmX1 <- data.frame(x = 17.044654, y = 0.00, xend = 17.044654, yend = 0.045, sequence_NEW = "X")
data.segmX2 <- data.frame(x = 20.200540, y = 0.00, xend = 20.200540, yend = 0.045, sequence_NEW = "X")
data.segmX3 <- data.frame(x = 3.993915, y = 0.00, xend = 3.993915, yend = 0.045, sequence_NEW = "X")
data.segmX4 <- data.frame(x = 11.174443, y = 0.00, xend = 11.174443, yend = 0.045, sequence_NEW = "X")
data.segmX5 <- data.frame(x = 6.042373, y = 0.00, xend = 6.042373, yend = 0.045, sequence_NEW = "X")
data.segmX6 <- data.frame(x = 20.201708, y = 0.00, xend = 20.201708, yend = 0.045, sequence_NEW = "X")
data.segm41 <- data.frame(x = 3.246708, y = 0.00, xend = 3.246708, yend = 0.045, sequence_NEW = "4")
data.segm42 <- data.frame(x = 17.382172, y = 0.00, xend = 17.382172, yend = 0.045, sequence_NEW = "4")
data.segm51 <- data.frame(x = 9.422591, y = 0.00, xend = 9.422591, yend = 0.045, sequence_NEW = "5")
data.segm52 <- data.frame(x = 18.649749, y = 0.00, xend = 18.649749, yend = 0.045, sequence_NEW = "5")

# Reorder the sequence_NEW factor
dxy1$sequence_NEW <- factor(dxy1$sequence_NEW, levels = c("X", "2L", "2R", "3", "4", "5"))
dxy2$sequence_NEW <- factor(dxy2$sequence_NEW, levels = c("X", "2L", "2R", "3", "4", "5"))
dxy3$sequence_NEW <- factor(dxy3$sequence_NEW, levels = c("X", "2L", "2R", "3", "4", "5"))

# Combine the segment data frames with the main data frame
dxy1 <- cbind(dxy1, sequence_NEW_order = dxy1$sequence_NEW)
dxy2 <- cbind(dxy2, sequence_NEW_order = dxy2$sequence_NEW)
dxy3 <- cbind(dxy3, sequence_NEW_order = dxy3$sequence_NEW)
data.segm <- rbind(data.segmX1, data.segmX2, data.segmX3, data.segmX4, data.segmX5, data.segmX6, data.segm41, data.segm42, data.segm51, data.segm52)

data.segm$sequence_NEW <- fct_relevel(data.segm$sequence_NEW, "X", "2L", "2R", "3", "4", "5")

#NOORA tässä vaihdetaan värejä rivinumeron perusteella. Eli katso data.segm dataframesta minkä viivan haluat muuttaa
data.segm$dashGroup <- "1" #the dash group defines if a line is dashed or not (1 = solid, 2 = dashed)
data.segm[1, "dashGroup"] <- "3" #make a specific line dashed by changing their group. Change to whatever row in the dataframe needs to be dashed
data.segm[2, "dashGroup"] <- "3" #make a specific line dashed by changing their group. Change to whatever row in the dataframe needs to be dashed

data.segm$colorGroup <- "1" #the color group defines the color
data.segm[1, "colorGroup"] <- "1" #same logic as with dashGroup
data.segm[2, "colorGroup"] <- "3"
data.segm[3, "colorGroup"] <- "2"
data.segm[4, "colorGroup"] <- "2"
data.segm[5, "colorGroup"] <- "1"
data.segm[6, "colorGroup"] <- "1"
data.segm[7, "colorGroup"] <- "1"
data.segm[8, "colorGroup"] <- "1"
data.segm[9, "colorGroup"] <- "1"
data.segm[10, "colorGroup"] <- "1"

global_size <- 0.5
global_size2 <- 0.5
global_alpha <- 0.6

a <- ggplot(data = dxy1, aes(x = start_NEW, y = d_xy)) +
  labs(x = "Position (Mb)", y = "Genetic divergence (dxy)", tag = "") +
  ggtitle(label = 'Inter- and intraspecific comparisons') +
  facet_grid(cols = vars(sequence_NEW), scales = "free_x", space = "free") +
  stat_summary(data = dxy1, geom = "line", size = global_size, alpha = 1, color = "#7D0DDB") +
  stat_summary(data = dxy2, geom = "line", size = global_size, alpha = 1, color = "#6CD1AC") +
  stat_summary(data = dxy3, geom = "line", size = global_size, alpha = 1, color = "#3394A8") +
  coord_cartesian(ylim = c(0.00, 0.044)) +
  scale_x_continuous(labels = scaleFUN, breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02, 0.03, 0.04)) +
  geom_segment(data = data.segm, aes(x = x, y = y, yend = yend, xend = xend, linetype = dashGroup, color = colorGroup), 
               size = global_size2, alpha = global_alpha) +
  scale_color_manual(values = c("#CC6600", "#3300FF", "#993300"),
                     breaks = c("1", "2", "3")) + #these refer to colorGroups
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="bold",hjust = 0))
a

Fig5 <- a

ggsave(Fig5, filename = "C:/Users/Noora/Desktop/Fig5.png", height = 3, width = 8)



#### Genome-wide Fst plot ####

# Read the data frames
dxy1 <- read.table("2montana.flavomontana.allopatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header = TRUE)
dxy2 <- read.table("2montana.flavomontana.sympatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header = TRUE)
dxy3 <- read.table("2montana.allopatric.sympatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header = TRUE)

# Scale the start column
dxy1$start_NEW <- dxy1$start / 1000000
dxy2$start_NEW <- dxy2$start / 1000000
dxy3$start_NEW <- dxy3$start / 1000000

# Replace sequence values
dxy_replace_sequence <- function(df) {
  df$sequence[df$sequence == 'X_chromosome'] <- 'X'
  df$sequence[df$sequence == '2L_chromosome'] <- '2L'
  df$sequence[df$sequence == '2R_chromosome'] <- '2R'
  df$sequence[df$sequence == '3_chromosome'] <- '3'
  df$sequence[df$sequence == '4_chromosome'] <- '4'
  df$sequence[df$sequence == '5_chromosome'] <- '5'
  df$sequence <- factor(df$sequence, levels = c("X", "2L", "2R", "3", "4", "5"))
  df$sequence_NEW <- factor(df$sequence, levels = c("X", "2L", "2R", "3", "4", "5"), labels = c("X", "2L", "2R", "3", "4", "5"))
  return(df)
}

dxy1 <- dxy_replace_sequence(dxy1)
dxy2 <- dxy_replace_sequence(dxy2)
dxy3 <- dxy_replace_sequence(dxy3)

# Create segment data frame
data.segmX1 <- data.frame(x = 17.044654, y = 0.00, xend = 17.044654, yend = 0.9, sequence_NEW = "X")
data.segmX2 <- data.frame(x = 20.200540, y = 0.00, xend = 20.200540, yend = 0.9, sequence_NEW = "X")
data.segmX3 <- data.frame(x = 3.993915, y = 0.00, xend = 3.993915, yend = 0.9, sequence_NEW = "X")
data.segmX4 <- data.frame(x = 11.174443, y = 0.00, xend = 11.174443, yend = 0.9, sequence_NEW = "X")
data.segmX5 <- data.frame(x = 6.042373, y = 0.00, xend = 6.042373, yend = 0.9, sequence_NEW = "X")
data.segmX6 <- data.frame(x = 20.201708, y = 0.00, xend = 20.201708, yend = 0.9, sequence_NEW = "X")
data.segm41 <- data.frame(x = 3.246708, y = 0.00, xend = 3.246708, yend = 0.9, sequence_NEW = "4")
data.segm42 <- data.frame(x = 17.382172, y = 0.00, xend = 17.382172, yend = 0.9, sequence_NEW = "4")
data.segm51 <- data.frame(x = 9.422591, y = 0.00, xend = 9.422591, yend = 0.9, sequence_NEW = "5")
data.segm52 <- data.frame(x = 18.649749, y = 0.00, xend = 18.649749, yend = 0.9, sequence_NEW = "5")

# Reorder the sequence_NEW factor
dxy1$sequence_NEW <- factor(dxy1$sequence_NEW, levels = c("X", "2L", "2R", "3", "4", "5"))
dxy2$sequence_NEW <- factor(dxy2$sequence_NEW, levels = c("X", "2L", "2R", "3", "4", "5"))
dxy3$sequence_NEW <- factor(dxy3$sequence_NEW, levels = c("X", "2L", "2R", "3", "4", "5"))

# Combine the segment data frames with the main data frame
dxy1 <- cbind(dxy1, sequence_NEW_order = dxy1$sequence_NEW)
dxy2 <- cbind(dxy2, sequence_NEW_order = dxy2$sequence_NEW)
dxy3 <- cbind(dxy3, sequence_NEW_order = dxy3$sequence_NEW)
data.segm <- rbind(data.segmX1, data.segmX2, data.segmX3, data.segmX4, data.segmX5, data.segmX6, data.segm41, data.segm42, data.segm51, data.segm52)

data.segm$sequence_NEW <- fct_relevel(data.segm$sequence_NEW, "X", "2L", "2R", "3", "4", "5")

data.segm$dashGroup <- "1" #the dash group defines if a line is dashed or not (1 = solid, 2 = dashed)
data.segm[1, "dashGroup"] <- "3" #make a specific line dashed by changing their group. Change to whatever row in the dataframe needs to be dashed
data.segm[2, "dashGroup"] <- "3" #make a specific line dashed by changing their group. Change to whatever row in the dataframe needs to be dashed

data.segm$colorGroup <- "1" #the color group defines the color
data.segm[1, "colorGroup"] <- "1" #same logic as with dashGroup
data.segm[2, "colorGroup"] <- "3"
data.segm[3, "colorGroup"] <- "2"
data.segm[4, "colorGroup"] <- "2"
data.segm[5, "colorGroup"] <- "1"
data.segm[6, "colorGroup"] <- "1"
data.segm[7, "colorGroup"] <- "1"
data.segm[8, "colorGroup"] <- "1"
data.segm[9, "colorGroup"] <- "1"
data.segm[10, "colorGroup"] <- "1"

global_size <- 0.5
global_size2 <- 0.5
global_alpha <- 0.6

a <- ggplot(data = dxy1, aes(x = start_NEW, y = f_st)) +
  labs(x = "Position (Mb)", y = "Genetic differentiation (Fst)", tag = "") +
  ggtitle(label = 'Inter- and intraspecific comparisons') +
  facet_grid(cols = vars(sequence_NEW), scales = "free_x", space = "free") +
  stat_summary(data = dxy1, geom = "line", size = global_size, alpha = 1, color = "#7D0DDB") +
  stat_summary(data = dxy2, geom = "line", size = global_size, alpha = 1, color = "#6CD1AC") +
  stat_summary(data = dxy3, geom = "line", size = global_size, alpha = 1, color = "#3394A8") +
  coord_cartesian(ylim = c(0.0, 0.9)) +
  scale_x_continuous(labels = scaleFUN, breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) +
  geom_segment(data = data.segm, aes(x = x, y = y, yend = yend, xend = xend, linetype = dashGroup, color = colorGroup), 
               size = global_size2, alpha = global_alpha) +
  scale_color_manual(values = c("#CC6600", "#3300FF", "#999999", "#993300"),
                     breaks = c("1", "2", "3")) + #these refer to colorGroups
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="bold",hjust = 0))
a

Fig_fst <- a

ggsave(Fig_fst, filename = "C:/Users/Noora/Desktop/Fig_genomewide_fst.png", height = 3, width = 8)


#### Genome-wide heterozygosity plot ####

# Read the data frames
dxy1 <- read.table("2montana.flavomontana.allopatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header = TRUE)
dxy2 <- read.table("2montana.flavomontana.sympatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header = TRUE)
dxy3 <- read.table("2montana.allopatric.sympatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header = TRUE)
dxy4 <- read.table("2montana.flavomontana.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header = TRUE)

# Scale the start column
dxy1$start_NEW <- dxy1$start / 1000000
dxy2$start_NEW <- dxy2$start / 1000000
dxy3$start_NEW <- dxy3$start / 1000000
dxy4$start_NEW <- dxy4$start / 1000000

# Replace sequence values
dxy_replace_sequence <- function(df) {
  df$sequence[df$sequence == 'X_chromosome'] <- 'X'
  df$sequence[df$sequence == '2L_chromosome'] <- '2L'
  df$sequence[df$sequence == '2R_chromosome'] <- '2R'
  df$sequence[df$sequence == '3_chromosome'] <- '3'
  df$sequence[df$sequence == '4_chromosome'] <- '4'
  df$sequence[df$sequence == '5_chromosome'] <- '5'
  df$sequence <- factor(df$sequence, levels = c("X", "2L", "2R", "3", "4", "5"))
  df$sequence_NEW <- factor(df$sequence, levels = c("X", "2L", "2R", "3", "4", "5"), labels = c("X", "2L", "2R", "3", "4", "5"))
  return(df)
}

dxy1 <- dxy_replace_sequence(dxy1)
dxy2 <- dxy_replace_sequence(dxy2)
dxy3 <- dxy_replace_sequence(dxy3)
dxy4 <- dxy_replace_sequence(dxy4)

# Create segment data frame
data.segmX1 <- data.frame(x = 17.044654, y = 0.00, xend = 17.044654, yend = 0.02, sequence_NEW = "X")
data.segmX2 <- data.frame(x = 20.200540, y = 0.00, xend = 20.200540, yend = 0.02, sequence_NEW = "X")
data.segmX3 <- data.frame(x = 3.993915, y = 0.00, xend = 3.993915, yend = 0.02, sequence_NEW = "X")
data.segmX4 <- data.frame(x = 11.174443, y = 0.00, xend = 11.174443, yend = 0.02, sequence_NEW = "X")
data.segmX5 <- data.frame(x = 6.042373, y = 0.00, xend = 6.042373, yend = 0.02, sequence_NEW = "X")
data.segmX6 <- data.frame(x = 20.201708, y = 0.00, xend = 20.201708, yend = 0.02, sequence_NEW = "X")
data.segm41 <- data.frame(x = 3.246708, y = 0.00, xend = 3.246708, yend = 0.02, sequence_NEW = "4")
data.segm42 <- data.frame(x = 17.382172, y = 0.00, xend = 17.382172, yend = 0.02, sequence_NEW = "4")
data.segm51 <- data.frame(x = 9.422591, y = 0.00, xend = 9.422591, yend = 0.02, sequence_NEW = "5")
data.segm52 <- data.frame(x = 18.649749, y = 0.00, xend = 18.649749, yend = 0.02, sequence_NEW = "5")

# Reorder the sequence_NEW factor
dxy1$sequence_NEW <- factor(dxy1$sequence_NEW, levels = c("X", "2L", "2R", "3", "4", "5"))
dxy2$sequence_NEW <- factor(dxy2$sequence_NEW, levels = c("X", "2L", "2R", "3", "4", "5"))
dxy3$sequence_NEW <- factor(dxy3$sequence_NEW, levels = c("X", "2L", "2R", "3", "4", "5"))

# Combine the segment data frames with the main data frame
dxy1 <- cbind(dxy1, sequence_NEW_order = dxy1$sequence_NEW)
dxy2 <- cbind(dxy2, sequence_NEW_order = dxy2$sequence_NEW)
dxy3 <- cbind(dxy3, sequence_NEW_order = dxy3$sequence_NEW)
data.segm <- rbind(data.segmX1, data.segmX2, data.segmX3, data.segmX4, data.segmX5, data.segmX6, data.segm41, data.segm42, data.segm51, data.segm52)

data.segm$sequence_NEW <- fct_relevel(data.segm$sequence_NEW, "X", "2L", "2R", "3", "4", "5")

data.segm$dashGroup <- "1" #the dash group defines if a line is dashed or not (1 = solid, 2 = dashed)
data.segm[1, "dashGroup"] <- "3" #make a specific line dashed by changing their group. Change to whatever row in the dataframe needs to be dashed
data.segm[2, "dashGroup"] <- "3" #make a specific line dashed by changing their group. Change to whatever row in the dataframe needs to be dashed

data.segm$colorGroup <- "1" #the color group defines the color
data.segm[1, "colorGroup"] <- "1" #same logic as with dashGroup
data.segm[2, "colorGroup"] <- "3"
data.segm[3, "colorGroup"] <- "2"
data.segm[4, "colorGroup"] <- "2"
data.segm[5, "colorGroup"] <- "1"
data.segm[6, "colorGroup"] <- "1"
data.segm[7, "colorGroup"] <- "1"
data.segm[8, "colorGroup"] <- "1"
data.segm[9, "colorGroup"] <- "1"
data.segm[10, "colorGroup"] <- "1"

global_size <- 0.3
global_size2 <- 0.5
global_alpha <- 0.6

a <- ggplot(data = dxy1, aes(x = start_NEW)) +
  labs(x = "Position (Mb)", y = "Heterozygosity", tag = "") +
  ggtitle(label = '') +
  facet_grid(cols = vars(sequence_NEW), scales = "free_x", space = "free") +
  stat_summary(data = dxy1, aes(x=start_NEW, y=heterozygosity_A), geom = "line", size = global_size, alpha = 1, color = "#E69F00") +
  stat_summary(data = dxy2, aes(x=start_NEW, y=heterozygosity_B), geom = "line", size = global_size, alpha = 1, color = "#6CD1AC") +
  stat_summary(data = dxy3, aes(x=start_NEW, y=heterozygosity_A), geom = "line", size = global_size, alpha = 1, color = "#7D0DDB") +
  stat_summary(data = dxy4, aes(x=start_NEW, y=heterozygosity_B), geom = "line", size = global_size, alpha = 1, color = "#000066") +
  coord_cartesian(ylim = c(0.00, 0.02)) +
  scale_x_continuous(labels = scaleFUN, breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  geom_segment(data = data.segm, aes(x = x, y = y, yend = yend, xend = xend, linetype = dashGroup, color = colorGroup), 
               size = global_size2, alpha = global_alpha) +
  scale_color_manual(values = c("#CC6600", "#3300FF", "#993300"),
                     breaks = c("1", "2", "3")) + #these refer to colorGroups
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="bold",hjust = 0))
a

Fig_het <- a

ggsave(Fig_het, filename = "C:/Users/Noora/Desktop/Fig_genomewide_het.png", height = 3, width = 8)



#### Genome-wide dxy plot; chromosomes zoomed in ####

dxy1<-read.table("2montana.flavomontana.allopatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header=T)
dxy2<-read.table("2montana.flavomontana.sympatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header=T)
dxy3<-read.table("2montana.allopatric.sympatric.intergenic_regions_for_plotting.64bp.windows.popgen.bed", header=T)

dxy1$start_NEW <- dxy1$start / 1000000
dxy2$start_NEW <- dxy2$start / 1000000
dxy3$start_NEW <- dxy3$start / 1000000

dxy1$sequence[dxy1$sequence == 'X_chromosome'] <- 'X'
dxy1$sequence[dxy1$sequence == '2L_chromosome'] <- '2L'
dxy1$sequence[dxy1$sequence == '2R_chromosome'] <- '2R'
dxy1$sequence[dxy1$sequence == '3_chromosome'] <- '3'
dxy1$sequence[dxy1$sequence == '4_chromosome'] <- '4'
dxy1$sequence[dxy1$sequence == '5_chromosome'] <- '5'

dxy2$sequence[dxy2$sequence == 'X_chromosome'] <- 'X'
dxy2$sequence[dxy2$sequence == '2L_chromosome'] <- '2L'
dxy2$sequence[dxy2$sequence == '2R_chromosome'] <- '2R'
dxy2$sequence[dxy2$sequence == '3_chromosome'] <- '3'
dxy2$sequence[dxy2$sequence == '4_chromosome'] <- '4'
dxy2$sequence[dxy2$sequence == '5_chromosome'] <- '5'

dxy3$sequence[dxy3$sequence == 'X_chromosome'] <- 'X'
dxy3$sequence[dxy3$sequence == '2L_chromosome'] <- '2L'
dxy3$sequence[dxy3$sequence == '2R_chromosome'] <- '2R'
dxy3$sequence[dxy3$sequence == '3_chromosome'] <- '3'
dxy3$sequence[dxy3$sequence == '4_chromosome'] <- '4'
dxy3$sequence[dxy3$sequence == '5_chromosome'] <- '5'

data.segmX1<-data.frame(x=17.044654,y=0.00,xend=17.044654,yend=0.045, sequence="X")
data.segmX2<-data.frame(x=20.200540,y=0.00,xend=20.200540,yend=0.045, sequence="X")
data.segmX3<-data.frame(x=3.993915,y=0.00,xend=3.993915,yend=0.045, sequence="X")
data.segmX4<-data.frame(x=11.174443,y=0.00,xend=11.174443,yend=0.045, sequence="X")
data.segmX5<-data.frame(x=6.042373,y=0.00,xend=6.042373,yend=0.045, sequence="X")
data.segmX6<-data.frame(x=20.201708,y=0.00,xend=20.201708,yend=0.045, sequence="X")

data.segm41<-data.frame(x=3.246708,y=0.00,xend=3.246708,yend=0.045, sequence="4")
data.segm42<-data.frame(x=17.382172,y=0.00,xend=17.382172,yend=0.045, sequence="4")

data.segm51<-data.frame(x=9.422591,y=0.00,xend=9.422591,yend=0.045, sequence="5")
data.segm52<-data.frame(x=18.649749,y=0.00,xend=18.649749,yend=0.045, sequence="5")

scaleFUN <- function(x) sprintf("%.0f", x)

global_size <- 0.5
global_size2 <- 0.5
global_alpha <- 0.5

dxy1_X <- subset(dxy1, sequence=="X")
dxy2_X <- subset(dxy2, sequence=="X")
dxy3_X <- subset(dxy3, sequence=="X")

a<-ggplot(data = dxy1_X, aes(x=start_NEW, y=d_xy)) + 
  labs(x="Position (Mb)", y="Genetic divergence (dxy)", tag="") +
  ggtitle(label = '') +
  facet_grid(cols=vars(sequence), scales="free_x", space = "free") +
  stat_summary(data=dxy1_X, geom = "line", size = global_size, alpha=1, color="#7D0DDB") + 
  stat_summary(data=dxy2_X, geom = "line", size = global_size, alpha=1, color="#6CD1AC") + 
  stat_summary(data=dxy3_X, geom = "line", size = global_size, alpha=1, color="#3394A8") +
  coord_cartesian(ylim=c(0.00,0.044)) +
  scale_x_continuous(labels=scaleFUN, breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02, 0.03, 0.04)) +
  geom_segment(data=data.segmX1,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, 
               color="#CC6600", linetype=2, linewidth=global_size2, alpha=global_alpha) +
  geom_segment(data=data.segmX2,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, 
               color="#993300", linetype=2, linewidth=global_size2, alpha=global_alpha) +
  geom_segment(data=data.segmX3,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, 
               color="#3300FF", linetype=1, linewidth=global_size2, alpha=global_alpha) +
  geom_segment(data=data.segmX4,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, 
               color="#3300FF", linetype=1, linewidth=global_size2, alpha=global_alpha) +
  geom_segment(data=data.segmX5,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, 
               color="#CC6600", linetype=1, linewidth=global_size2, alpha=global_alpha) +
  geom_segment(data=data.segmX6,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, 
               color="#CC6600", linetype=1, linewidth=global_size2, alpha=global_alpha) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
a

dxy1_2L <- subset(dxy1, sequence=="2L")
dxy2_2L <- subset(dxy2, sequence=="2L")
dxy3_2L <- subset(dxy3, sequence=="2L")

b<-ggplot(data = dxy1_2L, aes(x=start_NEW, y=d_xy)) + 
  labs(x="Position (Mb)", y="Genetic divergence (dxy)", tag="") +
  ggtitle(label = '') +
  facet_grid(cols=vars(sequence), scales="free_x", space = "free") +
  stat_summary(data=dxy1_2L, geom = "line", size = global_size, alpha=1, color="#7D0DDB") + 
  stat_summary(data=dxy2_2L, geom = "line", size = global_size, alpha=1, color="#6CD1AC") + 
  stat_summary(data=dxy3_2L, geom = "line", size = global_size, alpha=1, color="#3394A8") +
  coord_cartesian(ylim=c(0.00,0.044)) +
  scale_x_continuous(labels=scaleFUN, breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02, 0.03, 0.04)) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
b


dxy1_2R <- subset(dxy1, sequence=="2R")
dxy2_2R <- subset(dxy2, sequence=="2R")
dxy3_2R <- subset(dxy3, sequence=="2R")

c<-ggplot(data = dxy1_2R, aes(x=start_NEW, y=d_xy)) + 
  labs(x="Position (Mb)", y="Genetic divergence (dxy)", tag="") +
  ggtitle(label = '') +
  facet_grid(cols=vars(sequence), scales="free_x", space = "free") +
  stat_summary(data=dxy1_2R, geom = "line", size = global_size, alpha=1, color="#7D0DDB") + 
  stat_summary(data=dxy2_2R, geom = "line", size = global_size, alpha=1, color="#6CD1AC") + 
  stat_summary(data=dxy3_2R, geom = "line", size = global_size, alpha=1, color="#3394A8") +
  coord_cartesian(ylim=c(0.00,0.044)) +
  scale_x_continuous(labels=scaleFUN, breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02, 0.03, 0.04)) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
c


dxy1_3 <- subset(dxy1, sequence=="3")
dxy2_3 <- subset(dxy2, sequence=="3")
dxy3_3 <- subset(dxy3, sequence=="3")

d<-ggplot(data = dxy1_3, aes(x=start_NEW, y=d_xy)) + 
  labs(x="Position (Mb)", y="Genetic divergence (dxy)", tag="") +
  ggtitle(label = '') +
  facet_grid(cols=vars(sequence), scales="free_x", space = "free") +
  stat_summary(data=dxy1_3, geom = "line", size = global_size, alpha=1, color="#7D0DDB") + 
  stat_summary(data=dxy2_3, geom = "line", size = global_size, alpha=1, color="#6CD1AC") + 
  stat_summary(data=dxy3_3, geom = "line", size = global_size, alpha=1, color="#3394A8") +
  coord_cartesian(ylim=c(0.00,0.044)) +
  scale_x_continuous(labels=scaleFUN, breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02, 0.03, 0.04)) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
d


dxy1_4 <- subset(dxy1, sequence=="4")
dxy2_4 <- subset(dxy2, sequence=="4")
dxy3_4 <- subset(dxy3, sequence=="4")

e<-ggplot(data = dxy1_4, aes(x=start_NEW, y=d_xy)) + 
  labs(x="Position (Mb)", y="Genetic divergence (dxy)", tag="") +
  ggtitle(label = '') +
  facet_grid(cols=vars(sequence), scales="free_x", space = "free") +
  stat_summary(data=dxy1_4, geom = "line", size = global_size, alpha=1, color="#7D0DDB") + 
  stat_summary(data=dxy2_4, geom = "line", size = global_size, alpha=1, color="#6CD1AC") + 
  stat_summary(data=dxy3_4, geom = "line", size = global_size, alpha=1, color="#3394A8") +
  coord_cartesian(ylim=c(0.00,0.044)) +
  scale_x_continuous(labels=scaleFUN, breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02, 0.03, 0.04)) +
  geom_segment(data=data.segm41,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, 
               color="#CC6600", linewidth=global_size2, alpha=global_alpha) +
  geom_segment(data=data.segm42,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, 
               color="#CC6600", linewidth=global_size2, alpha=global_alpha) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
e



dxy1_5 <- subset(dxy1, sequence=="5")
dxy2_5 <- subset(dxy2, sequence=="5")
dxy3_5 <- subset(dxy3, sequence=="5")

f<-ggplot(data = dxy1_5, aes(x=start_NEW, y=d_xy)) + 
  labs(x="Position (Mb)", y="Genetic divergence (dxy)", tag="") +
  ggtitle(label = '') +
  facet_grid(cols=vars(sequence), scales="free_x", space = "free") +
  stat_summary(data=dxy1_5, geom = "line", size = global_size, alpha=1, color="#7D0DDB") + 
  stat_summary(data=dxy2_5, geom = "line", size = global_size, alpha=1, color="#6CD1AC") + 
  stat_summary(data=dxy3_5, geom = "line", size = global_size, alpha=1, color="#3394A8") +
  coord_cartesian(ylim=c(0.00,0.044)) +
  scale_x_continuous(labels=scaleFUN, breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(breaks = c(0.00, 0.01, 0.02, 0.03, 0.04)) +
  geom_segment(data=data.segm51,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, 
               color="#CC6600", linewidth=global_size2, alpha=global_alpha) +
  geom_segment(data=data.segm52,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, 
               color="#CC6600", size=global_size2, alpha=global_alpha) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
f


Fig<-(a|b|c)/(d|e|f)

ggsave(Fig, filename = "C:/Users/Noora/Desktop/Fig_dxy_chroms_separately.png", height = 6, width = 8)





#### Repeats along the genome ####

# Read the data frames
dxy1 <- read.table("D_montana.ref_genome.unmasked.repeats_for_plotting.txt", header = TRUE)

# Scale the start column
dxy1$start_NEW <- dxy1$start / 1000000

# Replace sequence values
dxy_replace_sequence <- function(df) {
  df$sequence[df$sequence == 'X_chromosome'] <- 'X'
  df$sequence[df$sequence == '2L_chromosome'] <- '2L'
  df$sequence[df$sequence == '2R_chromosome'] <- '2R'
  df$sequence[df$sequence == '3_chromosome'] <- '3'
  df$sequence[df$sequence == '4_chromosome'] <- '4'
  df$sequence[df$sequence == '5_chromosome'] <- '5'
  df$sequence <- factor(df$sequence, levels = c("X", "2L", "2R", "3", "4", "5"))
  df$sequence_NEW <- factor(df$sequence, levels = c("X", "2L", "2R", "3", "4", "5"), labels = c("X", "2L", "2R", "3", "4", "5"))
  return(df)
}

dxy1 <- dxy_replace_sequence(dxy1)

# Create segment data frame
data.segmX1 <- data.frame(x = 17.044654, y = 0.00, xend = 17.044654, yend = 1800, sequence_NEW = "X")
data.segmX2 <- data.frame(x = 20.200540, y = 0.00, xend = 20.200540, yend = 1800, sequence_NEW = "X")
data.segmX3 <- data.frame(x = 3.993915, y = 0.00, xend = 3.993915, yend = 1800, sequence_NEW = "X")
data.segmX4 <- data.frame(x = 11.174443, y = 0.00, xend = 11.174443, yend = 1800, sequence_NEW = "X")
data.segmX5 <- data.frame(x = 6.042373, y = 0.00, xend = 6.042373, yend = 1800, sequence_NEW = "X")
data.segmX6 <- data.frame(x = 20.201708, y = 0.00, xend = 20.201708, yend = 1800, sequence_NEW = "X")
data.segm41 <- data.frame(x = 3.246708, y = 0.00, xend = 3.246708, yend = 1800, sequence_NEW = "4")
data.segm42 <- data.frame(x = 17.382172, y = 0.00, xend = 17.382172, yend = 1800, sequence_NEW = "4")
data.segm51 <- data.frame(x = 9.422591, y = 0.00, xend = 9.422591, yend = 1800, sequence_NEW = "5")
data.segm52 <- data.frame(x = 18.649749, y = 0.00, xend = 18.649749, yend = 1800, sequence_NEW = "5")

data.segm_a <- data.frame(x = 0, y = 1314.167, xend = 29.14, yend = 1314.167, sequence_NEW = "X")
data.segm_b <- data.frame(x = 0, y = 849.7815, xend = 20.25, yend = 849.7815, sequence_NEW = "2L")
data.segm_c <- data.frame(x = 0, y = 849.7815, xend = 11.00, yend = 849.7815, sequence_NEW = "2R")
data.segm_d <- data.frame(x = 0, y = 849.7815, xend = 26.02, yend = 849.7815, sequence_NEW = "3")
data.segm_e <- data.frame(x = 0, y = 849.7815, xend = 32.54, yend = 849.7815, sequence_NEW = "4")
data.segm_f <- data.frame(x = 0, y = 849.7815, xend = 26.51, yend = 849.7815, sequence_NEW = "5")

# Reorder the sequence_NEW factor
dxy1$sequence_NEW <- factor(dxy1$sequence_NEW, levels = c("X", "2L", "2R", "3", "4", "5"))

# Combine the segment data frames with the main data frame
dxy1 <- cbind(dxy1, sequence_NEW_order = dxy1$sequence_NEW)
data.segm <- rbind(data.segmX1, data.segmX2, data.segmX3, data.segmX4, data.segmX5, data.segmX6, data.segm41, data.segm42, data.segm51, data.segm52, data.segm_a, data.segm_b, data.segm_c, data.segm_d, data.segm_e, data.segm_f)

data.segm$sequence_NEW <- fct_relevel(data.segm$sequence_NEW, "X", "2L", "2R", "3", "4", "5")

data.segm$dashGroup <- "1" #the dash group defines if a line is dashed or not (1 = solid, 2 = dashed)
data.segm[1, "dashGroup"] <- "3" #make a specific line dashed by changing their group. Change to whatever row in the dataframe needs to be dashed
data.segm[2, "dashGroup"] <- "3" #make a specific line dashed by changing their group. Change to whatever row in the dataframe needs to be dashed

data.segm$colorGroup <- "1" #the color group defines the color
data.segm[1, "colorGroup"] <- "1" #same logic as with dashGroup
data.segm[2, "colorGroup"] <- "3"
data.segm[3, "colorGroup"] <- "2"
data.segm[4, "colorGroup"] <- "2"
data.segm[5, "colorGroup"] <- "1"
data.segm[6, "colorGroup"] <- "1"
data.segm[7, "colorGroup"] <- "1"
data.segm[8, "colorGroup"] <- "1"
data.segm[9, "colorGroup"] <- "1"
data.segm[10, "colorGroup"] <- "1"

data.segm[11, "colorGroup"] <- "4"
data.segm[12, "colorGroup"] <- "4"
data.segm[13, "colorGroup"] <- "4"
data.segm[14, "colorGroup"] <- "4"
data.segm[15, "colorGroup"] <- "4"
data.segm[16, "colorGroup"] <- "4"

scaleFUN <- function(x) sprintf("%.0f", x)

global_size <- 0.8
global_size2 <- 0.8
global_alpha <- 0.6

a <- ggplot(dxy1) +
  geom_density(aes(x=start_NEW, y=after_stat(count)), adjust=0.7) +
  facet_grid(cols=vars(sequence_NEW), scales="free_x", space = "free") +
  labs(x="Chromosome", y="Repeat density") +
  geom_segment(data = data.segm, aes(x = x, y = y, yend = yend, xend = xend, linetype = dashGroup, color = colorGroup), 
               linewidth = global_size2, alpha = global_alpha) +
  scale_color_manual(values = c("#CC6600", "#3300FF", "#993300", "black"),
                     breaks = c("1", "2", "3", "4")) + #these refer to colorGroups
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="bold",hjust = 0))
a

Fig_repeats <- a

ggsave(Fig_repeats, filename = "C:/Users/Noora/Desktop/Fig_genomewide_repeats.png", height = 3, width = 8)





#### Simulations: montana as reference genome ####

sim_allop_64bp<-read.table("Difference_in_likelihood_under_correct_and_wrong_model_autosomes_collinear_allopatry_64bp.monref.txt", header=T)

a<-ggplot(sim_allop_64bp, aes(x=Difference_in_likelihood, width=0.8)) + 
  labs(x="Difference in Likelihood", y="Frequency", fill='', tag="A") +
  ggtitle("Colinear autosomes, \n allopatric comparison, \n 64bp blocks") +
  geom_density(aes(Difference_in_likelihood), linetype=1, fill="#7D0DDB", alpha=0.2) +
  geom_vline(xintercept = 12868.67568, color="#7D0DDB", size=1.5, linetype=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(8000,15000)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.00075), labels = scales::scientific) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
a


sim_symp_64bp<-read.table("Difference_in_likelihood_under_correct_and_wrong_model_autosomes_collinear_sympatry_64bp.monref.txt", header=T)

c<-ggplot(sim_symp_64bp, aes(x=Difference_in_likelihood, width=0.8)) + 
  labs(x="Difference in Likelihood", y="Frequency", fill='', tag="C") +
  ggtitle("Colinear autosomes, \n sympatric comparison, \n 64bp blocks") +
  geom_density(aes(Difference_in_likelihood), linetype=1, fill="#6CD1AC", alpha=0.2) +
  geom_vline(xintercept = 47789.95747, color="#6CD1AC", size=1.5, linetype=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(33000,52000)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.00028), labels = scales::scientific) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
c



sim_allop_128bp<-read.table("Difference_in_likelihood_under_correct_and_wrong_model_autosomes_collinear_allopatry_128bp.monref.txt", header=T)

b<-ggplot(sim_allop_128bp, aes(x=Difference_in_likelihood, width=0.8)) + 
  labs(x="Difference in Likelihood", y="Frequency", fill='', tag="B") +
  ggtitle("Colinear autosomes, \n allopatric comparison, \n 128bp blocks") +
  geom_density(aes(Difference_in_likelihood), linetype=1, fill="#7D0DDB", alpha=0.2) +
  geom_vline(xintercept = 15401.17199, color="#7D0DDB", size=1.5, linetype=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,17000)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.00125), labels = scales::scientific) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
b


sim_symp_128bp<-read.table("Difference_in_likelihood_under_correct_and_wrong_model_autosomes_collinear_sympatry_128bp.monref.txt", header=T)

d<-ggplot(sim_symp_128bp, aes(x=Difference_in_likelihood, width=0.8)) + 
  labs(x="Difference in Likelihood", y="Frequency", fill='', tag="D") +
  ggtitle("Colinear autosomes, \n sympatric comparison, \n 128bp blocks") +
  geom_density(aes(Difference_in_likelihood), linetype=1, fill="#6CD1AC", alpha=0.2) +
  geom_vline(xintercept = 56289.81669, color="#6CD1AC", size=1.5, linetype=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(10000,60000)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.00043), labels = scales::scientific) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
d


sim_intra_128bp<-read.table("Difference_in_likelihood_under_correct_and_wrong_model_autosomes_collinear_intraspecific_128bp.monref.txt", header=T)

e<-ggplot(sim_intra_128bp, aes(x=Difference_in_likelihood, width=0.8)) + 
  labs(x="Difference in Likelihood", y="Frequency", fill='', tag="E") +
  ggtitle("Colinear autosomes, \n intraspecific comparison, \n 128bp blocks") +
  geom_density(aes(Difference_in_likelihood), linetype=1, fill="#3394A8", alpha=0.2) +
  geom_vline(xintercept = 14069.83183, color="#3394A8", size=1.5, linetype=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(10000,90000)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.00038), labels = scales::scientific) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
e


fig<-(a|b)/(c|d)/(e + plot_spacer())
ggsave(fig, filename = "Fig_S5.png", width = 6, height = 8)


#### Simulations: flavomontana as reference genome ####

sim_allop_64bp<-read.table("Difference_in_likelihood_under_correct_and_wrong_model_autosomes_collinear_allopatry_64bp.flaref.txt", header=T)

a<-ggplot(sim_allop_64bp, aes(x=Difference_in_likelihood, width=0.8)) + 
  labs(x="Difference in Likelihood", y="Frequency", fill='', tag="A") +
  ggtitle("Colinear autosomes, \n allopatric comparison, \n 64bp blocks") +
  geom_density(aes(Difference_in_likelihood), linetype=1, fill="#7D0DDB", alpha=0.2) +
  geom_vline(xintercept = 9057.080976, color="#7D0DDB", size=1.5, linetype=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(7000,14000)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.00064), labels = scales::scientific) + 
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
a


sim_symp_64bp<-read.table("Difference_in_likelihood_under_correct_and_wrong_model_autosomes_collinear_sympatry_64bp.flaref.txt", header=T)

c<-ggplot(sim_symp_64bp, aes(x=Difference_in_likelihood, width=0.8)) + 
  labs(x="Difference in Likelihood", y="Frequency", fill='', tag="C") +
  ggtitle("Colinear autosomes, \n sympatric comparison, \n 64bp blocks") +
  geom_density(aes(Difference_in_likelihood), linetype=1, fill="#6CD1AC", alpha=0.2) +
  geom_vline(xintercept = 35846.52767, color="#6CD1AC", size=1.5, linetype=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(25000,50000)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.00023), labels = scales::scientific) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
c



sim_allop_128bp<-read.table("Difference_in_likelihood_under_correct_and_wrong_model_autosomes_collinear_allopatry_128bp.flaref.txt", header=T)

b<-ggplot(sim_allop_128bp, aes(x=Difference_in_likelihood, width=0.8)) + 
  labs(x="Difference in Likelihood", y="Frequency", fill='', tag="B") +
  ggtitle("Colinear autosomes, \n allopatric comparison, \n 128bp blocks") +
  geom_density(aes(Difference_in_likelihood), linetype=1, fill="#7D0DDB", alpha=0.2) +
  geom_vline(xintercept = 11739.32471, color="#7D0DDB", size=1.5, linetype=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,13000)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.00141), labels = scales::scientific) + 
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
b


sim_symp_128bp<-read.table("Difference_in_likelihood_under_correct_and_wrong_model_autosomes_collinear_sympatry_128bp.flaref.txt", header=T)

d<-ggplot(sim_symp_128bp, aes(x=Difference_in_likelihood, width=0.8)) + 
  labs(x="Difference in Likelihood", y="Frequency", fill='', tag="D") +
  ggtitle("Colinear autosomes, \n sympatric comparison, \n 128bp blocks") +
  geom_density(aes(Difference_in_likelihood), linetype=1, fill="#6CD1AC", alpha=0.2) +
  geom_vline(xintercept = 44004.21189, color="#6CD1AC", size=1.5, linetype=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(5000,47000)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.00045), labels = scales::scientific) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
d


sim_intra_128bp<-read.table("Difference_in_likelihood_under_correct_and_wrong_model_autosomes_collinear_intraspecific_128bp.flaref.txt", header=T)

e<-ggplot(sim_intra_128bp, aes(x=Difference_in_likelihood, width=0.8)) + 
  labs(x="Difference in Likelihood", y="Frequency", fill='', tag="E") +
  ggtitle("Colinear autosomes, \n intraspecific comparison, \n 128bp blocks") +
  geom_density(aes(Difference_in_likelihood), linetype=1, fill="#3394A8", alpha=0.2) +
  geom_vline(xintercept = 17638.01461, color="#3394A8", size=1.5, linetype=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(10000,90000)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.0003), labels = scales::scientific) +
  theme(panel.background = element_rect(fill='white', colour='white')
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'), axis.text=element_text(size=9),
        axis.title = element_text(size=9),
        plot.title = element_text(color="black", size=9, face="italic",hjust = 0)) +
  theme(legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=9)) +
  theme(plot.title = element_text(color = 'black', size=9, face="bold"))
e


fig<-(a|b)/(c|d)/(e + plot_spacer())
ggsave(fig, filename = "C:/Users/Noora/Desktop/Fig_S6.png", width = 6, height = 8)



#### Divergence time and migration rate ####

test<-read.table("divergence_time_migration_rate_results.txt", header=T)

test$Genomic_region <- factor(test$Genomic_region, levels = c("autosomes_colinear", "4_inverted", "5_inverted",
                                                              "X_colinear", "X_inverted"))

test$Geography <- factor(test$Geography, levels = c("Allopatry", "Sympatry"))

levels(test$Genomic_region)[levels(test$Genomic_region)=="autosomes_colinear"] <- "COL"
levels(test$Genomic_region)[levels(test$Genomic_region)=="4_inverted"] <- "INV 4"
levels(test$Genomic_region)[levels(test$Genomic_region)=="5_inverted"] <- "INV 5"
levels(test$Genomic_region)[levels(test$Genomic_region)=="X_colinear"] <- "COL "
levels(test$Genomic_region)[levels(test$Genomic_region)=="X_inverted"] <- "INV X"

test1<-subset(test, Chromosome_type=="Autosomes")
test2<-subset(test, Chromosome_type=="X")

a<-ggplot(test) +
  aes(y = T_scaled_Mya, 
      x = Genomic_region,
      color = Geography) +
  geom_point(position=position_dodge(0.3), size=1.5)+
  geom_errorbar(aes(ymin=T_scaled_Mya-SD2_T_Mya, ymax=T_scaled_Mya+SD2_T_Mya), width=.2,
                position=position_dodge(0.3), size=0.6) +
  facet_grid(cols=vars(Chromosome_type), scales="free_x", space = "free") +
  scale_color_manual(values = c("#7D0DDB", "#5BCFA7")) +
  scale_y_continuous(limits = c(2.4, 3.5),
                     breaks=c(2.4, 2.6, 2.8, 3.0, 3.2, 3.4)) +
  labs(x = "",
       y = "Divergence time (Mya)",
       fill = "Geography",
       tag="A") +
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(size=13),
        plot.title = element_text(size = 13, face = "bold"),
        legend.text=element_text(size=13), 
        legend.title=element_text(size=13),
        axis.text=element_text(size=13, color="black"),
        strip.text.x = element_text(size = 13)) +
  theme(axis.text.x = element_text(angle = 35, vjust = 1.05, hjust=1))
a

b<-ggplot(test) +
  aes(y = me_unscaled, 
      x = Genomic_region,
      color = Geography) +
  geom_point(position=position_dodge(0.3), size=1.5)+
  geom_errorbar(aes(ymin=me_unscaled-SD2_M, ymax=me_unscaled+SD2_M), width=.2,
                position=position_dodge(0.3), size=0.6) +
  facet_grid(cols=vars(Chromosome_type), scales="free_x", space = "free") +
  scale_color_manual(values = c("#7D0DDB", "#5BCFA7")) +
  scale_y_continuous(limits = c(0.0014, 0.014),
                     breaks=c(0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014)) +
  labs(x = "",
       y = "Migration rate \n(migrants / generation)",
       fill = "Geography",
       tag="B") +
  theme_bw() +
  theme(legend.position="right",
        axis.title = element_text(size=13),
        plot.title = element_text(size = 13, face = "bold"),
        legend.text=element_text(size=13), 
        legend.title=element_text(size=13),
        axis.text=element_text(size=13, color="black"),
        strip.text.x = element_text(size = 13)) +
  theme(axis.text.x = element_text(angle = 35, vjust = 1.05, hjust=1))
b

fig<-a|b

ggsave(fig, filename = "C:/Users/Noora/Desktop/Fig6.png", width = 9, height = 3.5)




#### Half violin plot ####

test<-read.table("results_T_M_Ne_64bp.txt", header=T)

test$Genomic_region <- factor(test$Genomic_region, levels = c("autosomes_colinear", "4_inverted", "5_inverted",
                                                              "X_colinear", "X_inverted"))

test$Geography <- factor(test$Geography, levels = c("Allopatry", "Sympatry"))

levels(test$Genomic_region)[levels(test$Genomic_region)=="autosomes_colinear"] <- "COL"
levels(test$Genomic_region)[levels(test$Genomic_region)=="4_inverted"] <- "INV 4"
levels(test$Genomic_region)[levels(test$Genomic_region)=="5_inverted"] <- "INV 5"
levels(test$Genomic_region)[levels(test$Genomic_region)=="X_colinear"] <- "COL "
levels(test$Genomic_region)[levels(test$Genomic_region)=="X_inverted"] <- "INV X"

test1<-subset(test, Chromosome_type=="Autosomes")
test2<-subset(test, Chromosome_type=="X")

text1<-data.frame(x=1.75,y=4.3, Chromosome_type="Autosomes")
text2<-data.frame(x=2.25,y=4.3, Chromosome_type="Autosomes")

text3<-data.frame(x=2.75,y=4.3, Chromosome_type="Autosomes")
text4<-data.frame(x=3.25,y=4.3, Chromosome_type="Autosomes")

text5<-data.frame(x=1.75,y=4.3, Chromosome_type="X")
text6<-data.frame(x=2.25,y=4.3, Chromosome_type="X")

a<-ggplot(test) +
  aes(y = T_scaled_Mya, 
      x = Genomic_region,
      fill = Geography) +
  introdataviz::geom_split_violin(data=test1, alpha = .2, trim = FALSE, width=2, color="#CCCCCC") +
  introdataviz::geom_split_violin(data=test2, alpha = .2, trim = FALSE, width=0.75, color="#CCCCCC") +
  facet_grid(cols=vars(Chromosome_type), scales="free_x", space = "free") +
  geom_boxplot(width = .5, 
               outlier.shape = NA,
               alpha = 0.2, size=0.5) +
  scale_fill_manual(values = c("#7D0DDB", "#6CD1AC")) +
  scale_y_continuous(limits = c(3.1, 4.3),
                     breaks=c(3.00, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2)) +
  labs(x = "",
       y = "Divergence time (Mya)",
       fill = "Geography",
       tag="A") +
  theme_bw() +
  theme(legend.position="none",
        axis.title = element_text(size=13),
        plot.title = element_text(size = 13, face = "bold"),
        legend.text=element_text(size=13), 
        legend.title=element_text(size=13),
        axis.text=element_text(size=13, color="black"),
        strip.text.x = element_text(size = 13)) +
  theme(axis.text.x = element_text(angle = 35, vjust = 1.05, hjust=1))
a


text1<-data.frame(x=1.8,y=0.0175, Chromosome_type="Autosomes")
text2<-data.frame(x=2.2,y=0.0175, Chromosome_type="Autosomes")

text3<-data.frame(x=2.8,y=0.0175, Chromosome_type="Autosomes")
text4<-data.frame(x=3.2,y=0.0175, Chromosome_type="Autosomes")

text5<-data.frame(x=1.8,y=0.0175, Chromosome_type="X")
text6<-data.frame(x=2.2,y=0.0175, Chromosome_type="X")

data.segm1<-data.frame(x=0.8, y=0.01, xend=1.2, yend=0.01, Chromosome_type="Autosomes")
data.segm2<-data.frame(x=1.8, y=0.0053, xend=2.2, yend=0.0053, Chromosome_type="Autosomes")
data.segm3<-data.frame(x=2.8, y=0.0047, xend=3.2, yend=0.0047, Chromosome_type="Autosomes")
data.segm4<-data.frame(x=0.8, y=0.0026, xend=1.2, yend=0.0026, Chromosome_type="X")
data.segm5<-data.frame(x=1.8, y=0.0022, xend=2.2, yend=0.0022, Chromosome_type="X")

text1b<-data.frame(x=1, y=0.0091, Chromosome_type="Autosomes")
text2b<-data.frame(x=2, y=0.0044, Chromosome_type="Autosomes")
text3b<-data.frame(x=3, y=0.0038, Chromosome_type="Autosomes")
text4b<-data.frame(x=1, y=0.0016, Chromosome_type="X")
text5b<-data.frame(x=2, y=0.0012, Chromosome_type="X")

b<-ggplot(test) +
  aes(y = me_unscaled, 
      x = Genomic_region,
      fill = Geography) +
  introdataviz::geom_split_violin(data=test1, alpha = .2, trim = FALSE, width=1.5, color="#CCCCCC") +
  introdataviz::geom_split_violin(data=test2, alpha = .2, trim = FALSE, width=0.9, color="#CCCCCC") +
  facet_grid(cols=vars(Chromosome_type), scales="free_x", space = "free") +
  geom_boxplot(width = .4, 
               outlier.shape = NA,
               alpha = 0.2, size=0.5) +
  scale_fill_manual(values = c("#7D0DDB", "#6CD1AC")) +
  scale_color_manual(values = c("#7D0DDB", "#6CD1AC")) +
labs(x = "",
     y = "Migration rate (migrants / generation)",
     fill = "Geography",
     tag="B") +
  theme_bw() +
  theme(legend.position="right",
        axis.title = element_text(size=13),
        plot.title = element_text(size = 13, face = "bold"),
        legend.text=element_text(size=13), 
        legend.title=element_text(size=13),
        axis.text=element_text(size=13, color="black"),
        strip.text.x = element_text(size = 13)) +
  theme(axis.text.x = element_text(angle = 35, vjust = 1.05, hjust=1))
b

fig<-a|b

ggsave(fig, filename = "C:/Users/Noora/Desktop/half_violin_plot.png", width = 10, height = 4)


#### inter and intraspecific dxy; correlation ####

cor <- read.table("inter_intra_dxy_correlation.txt", header = TRUE)

legent_title<-"Region"

a<-ggplot(cor, aes(x=d_xy_inter, y=d_xy_intra, color=Region)) + 
  labs(x="dxy interspecific", y="dxy intraspecific") +
  geom_point(aes(color=Region, fill=Region, shape=Region)) +
  scale_fill_manual(legent_title, values = c("#FF6699", "#00CC99", "blue")) +
  scale_color_manual(legent_title, values = c("#FF6699", "#00CC99", "blue")) +
  scale_shape_manual(legent_title, values = c(21, 24, 22)) +
  geom_smooth(method='lm', formula= y~x, se=F) +
  theme_bw() +
  theme(legend.position="right",
        axis.title = element_text(size=11),
        plot.title = element_text(size = 11, face = "bold"),
        legend.text=element_text(size=11), 
        legend.title=element_text(size=11),
        axis.text=element_text(size=11, color="black"),
        strip.text.x = element_text(size = 11)) +
  annotate("text", x=0.0235, y=0.019, label="col: r = 0.276, P = <0.001", size=4, color="#FF6699") +
  annotate("text", x=0.0236, y=0.018, label="inv4: r = 0.647, P = <0.001", size=4, color="#009966") +
  annotate("text", x=0.0235, y=0.017, label="inv5: r = 0.140, P = 0.144", size=4, color="blue")
a

cor1<-subset(cor, Region=="col")
cor1<-subset(cor, Region=="inv4")
cor1<-subset(cor, Region=="inv5")

correlation<-cor.test(x=cor1$d_xy_inter, cor1$d_xy_intra, method="pearson")
print(correlation)
print(correlation$p.value)


