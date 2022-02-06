################################# Coverage plotter for Project - Clinical data ###################################
###################################### Damilola Oresegun ##########################################

# Clear workspace and close figures
rm(list=ls())
graphics.off()
##Get libraries
library(ggplot2)
setwd("~/")
source("multiplot.R")
########### Coverage plot for 3D7 ##########
pdf("testimage.pdf")
# Set working directory
setwd("~/Testing/Mapping/3D7_genome")

# Read in the coverage file
s3D7_data<-read.csv("3D7.cov",sep="\t",header=T)

# Extract features of the file which fall within the range of Chromosome 12
s3D7_12_var<-s3D7_data[Reduce(intersect,list(which(s3D7_data[,1]=="Pf3D7_12_v3"),which(s3D7_data[,2]>=40000),which(s3D7_data[,2]<=65000))),]

# Calculate the mean coverage of the entire chromosome 12
s3D7_mean<-mean(s3D7_data[which(s3D7_data[,1]=="Pf3D7_12_v3"),5])

# Set the names of the axes
Position=s3D7_12_var[,2]+250
Coverage=s3D7_12_var[,5]

# Plot the coverage of the range
s1<-ggplot(s3D7_12_var,aes(x=Position,y=Coverage))+
  geom_area(fill="orange",alpha=.2)+
  #use the mean coverage as a threshold in red
  geom_hline(mapping = NULL, data=NULL,color="red", yintercept=s3D7_mean, na.rm = FALSE,show.legend = TRUE)+
  # set the starting point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=46788, na.rm = FALSE, show.legend = TRUE)+
  # set the starting point of the inter-exon intron
  geom_vline(mapping = NULL, data =NULL,color="forest green",xintercept=47955, na.rm = FALSE, show.legend = TRUE)+
  # set the end point of the intron
  geom_vline(mapping = NULL, data =NULL,color="forest green",xintercept=48803, na.rm = FALSE, show.legend = TRUE)+
  #set the end point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=56805, na.rm = FALSE, show.legend = TRUE)+
  geom_line()+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  ggtitle("3D7")

#################################### Coverage plot for 7G8 ##################################
# Set working directory
setwd("~/Testing/Mapping/7G8_genome")

# Read in the coverage file
s7G8_data<-read.csv("7G8.cov",sep="\t",header=T)

# Extract features of the file which fall within the range of Chromosome 12
s7G8_12_var<-s7G8_data[Reduce(intersect,list(which(s7G8_data[,1]=="Pf3D7_12_v3"),which(s7G8_data[,2]>=40000),which(s7G8_data[,2]<=65000))),]

# Calculate the mean coverage of the entire chromosome 12
s7G8_mean<-mean(s7G8_data[which(s7G8_data[,1]=="Pf3D7_12_v3"),5])

# Set the names of the axes
Position=s7G8_12_var[,2]+250
Coverage=s7G8_12_var[,5]

# Plot the coverage of the range
s2<-ggplot(s7G8_12_var,aes(x=Position,y=Coverage))+
  geom_area(fill="green",alpha=.2)+
  #use the mean coverage as a threshold in red
  geom_hline(mapping = NULL, data=NULL,color="red", yintercept=s7G8_mean, na.rm = FALSE,show.legend = TRUE)+
  # set the starting point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=46788, na.rm = FALSE, show.legend = TRUE)+
  # set the starting point of the inter-exon intron
  geom_vline(mapping = NULL, data =NULL,color="pink",xintercept=47955, na.rm = FALSE, show.legend = TRUE)+
  # set the end point of the intron
  geom_vline(mapping = NULL, data =NULL,color="pink",xintercept=48803, na.rm = FALSE, show.legend = TRUE)+
  #set the end point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=56805, na.rm = FALSE, show.legend = TRUE)+
  geom_line()+
  ggtitle("7G8")

########### Coverage plot for DD2 ##########

# Set working directory
setwd("~/Testing/Mapping/DD2_genome")

# Read in the coverage file
sDD2_data<-read.csv("DD2.cov",sep="\t",header=T)

# Extract features of the file which fall within the range of Chromosome 12
sDD2_12_var<-sDD2_data[Reduce(intersect,list(which(sDD2_data[,1]=="Pf3D7_12_v3"),which(sDD2_data[,2]>=40000),which(sDD2_data[,2]<=65000))),]

# Calculate the mean coverage of the entire chromosome 12
sDD2_mean<-mean(sDD2_data[which(sDD2_data[,1]=="Pf3D7_12_v3"),5])
# Set the names of the axes
Position=sDD2_12_var[,2]+250
Coverage=sDD2_12_var[,5]

# Plot the coverage of the range
s3<-ggplot(sDD2_12_var,aes(x=Position,y=Coverage))+
  geom_area(fill="pink",alpha=.2)+
  #use the mean coverage as a threshold in red
  geom_hline(mapping = NULL, data=NULL,color="red", yintercept=sDD2_mean, na.rm = FALSE,show.legend = TRUE)+
  # set the starting point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=46788, na.rm = FALSE, show.legend = TRUE)+
  # set the starting point of the inter-exon intron
  geom_vline(mapping = NULL, data =NULL,color="forest green",xintercept=47955, na.rm = FALSE, show.legend = TRUE)+
  # set the end point of the intron
  geom_vline(mapping = NULL, data =NULL,color="forest green",xintercept=48803, na.rm = FALSE, show.legend = TRUE)+
  #set the end point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=56805, na.rm = FALSE, show.legend = TRUE)+
  geom_line()+
  ggtitle("DD2")

########### Coverage plot for GB4 ##########

# Set working directory
setwd("~/Testing/Mapping/GB4_genome")

# Read in the coverage file
sGB4_data<-read.csv("GB4.cov",sep="\t",header=T)

# Extract features of the file which fall within the range of Chromosome 12
sGB4_12_var<-sGB4_data[Reduce(intersect,list(which(sGB4_data[,1]=="Pf3D7_12_v3"),which(sGB4_data[,2]>=40000),which(sGB4_data[,2]<=65000))),]

# Calculate the mean coverage of the entire chromosome 12
sGB4_mean<-mean(sGB4_data[which(sGB4_data[,1]=="Pf3D7_12_v3"),5])
# Set the names of the axes
Position=sGB4_12_var[,2]+250
Coverage=sGB4_12_var[,5]

# Plot the coverage of the range
s4<-ggplot(sGB4_12_var,aes(x=Position,y=Coverage))+
  geom_area(fill="gold",alpha=.2)+
  #use the mean coverage as a threshold in red
  geom_hline(mapping = NULL, data=NULL,color="red", yintercept=sGB4_mean, na.rm = FALSE,show.legend = TRUE)+
  # set the starting point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=46788, na.rm = FALSE, show.legend = TRUE)+
  # set the starting point of the inter-exon intron
  geom_vline(mapping = NULL, data =NULL,color="forest green",xintercept=47955, na.rm = FALSE, show.legend = TRUE)+
  # set the end point of the intron
  geom_vline(mapping = NULL, data =NULL,color="forest green",xintercept=48803, na.rm = FALSE, show.legend = TRUE)+
  #set the end point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=56805, na.rm = FALSE, show.legend = TRUE)+
  geom_line()+
  ggtitle("GB4")

#################################### Coverage plot for HB3 ##################################
# Set working directory
setwd("~/Testing/Mapping/HB3_genome")

# Read in the coverage file
sHB3_data<-read.csv("HB3.cov",sep="\t",header=T)

# Extract features of the file which fall within the range of Chromosome 12
sHB3_12_var<-sHB3_data[Reduce(intersect,list(which(sHB3_data[,1]=="Pf3D7_12_v3"),which(sHB3_data[,2]>=40000),which(sHB3_data[,2]<=65000))),]

# Calculate the mean coverage of the entire chromosome 12
sHB3_mean<-mean(sHB3_data[which(sHB3_data[,1]=="Pf3D7_12_v3"),5])
# Set the names of the axes
Position=sHB3_12_var[,2]+250
Coverage=sHB3_12_var[,5]

# Plot the coverage of the range
s5<-ggplot(sHB3_12_var,aes(x=Position,y=Coverage))+
  geom_area(fill="purple",alpha=.2)+
  #use the mean coverage as a threshold in red
  geom_hline(mapping = NULL, data=NULL,color="red", yintercept=sHB3_mean, na.rm = FALSE,show.legend = TRUE)+
  # set the starting point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=46788, na.rm = FALSE, show.legend = TRUE)+
  # set the starting point of the inter-exon intron
  geom_vline(mapping = NULL, data =NULL,color="pink",xintercept=47955, na.rm = FALSE, show.legend = TRUE)+
  # set the end point of the intron
  geom_vline(mapping = NULL, data =NULL,color="pink",xintercept=48803, na.rm = FALSE, show.legend = TRUE)+
  #set the end point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=56805, na.rm = FALSE, show.legend = TRUE)+
  geom_line()+
  ggtitle("HB3")

#################################### Coverage plot for ERR022699 ##################################
# Set working directory
#setwd("~/Pictures/ERR022699_genome")

# Read in the coverage file
sERR022699_data<-read.csv("ERR022699.cov",sep="\t",header=T)

# Extract features of the file which fall within the range of Chromosome 12
sERR022699_12_var<-sERR022699_data[Reduce(intersect,list(which(sERR022699_data[,1]=="Pf3D7_12_v3"),which(sERR022699_data[,2]>=40000),which(sERR022699_data[,2]<=65000))),]

# Calculate the mean coverage of the entire chromosome 12
sERR022699_mean<-mean(sERR022699_data[which(sERR022699_data[,1]=="Pf3D7_12_v3"),5])
# Set the names of the axes
Position=sERR022699_12_var[,2]+250
Coverage=sERR022699_12_var[,5]

# Plot the coverage of the range
c1<-ggplot(sERR022699_12_var,aes(x=Position,y=Coverage))+
  geom_area(fill="red",alpha=.2)+
  #use the mean coverage as a threshold in red
  geom_hline(mapping = NULL, data=NULL,color="red", yintercept=sERR022699_mean, na.rm = FALSE,show.legend = TRUE)+
  # set the starting point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=46788, na.rm = FALSE, show.legend = TRUE)+
  # set the starting point of the inter-exon intron
  geom_vline(mapping = NULL, data =NULL,color="pink",xintercept=47955, na.rm = FALSE, show.legend = TRUE)+
  # set the end point of the intron
  geom_vline(mapping = NULL, data =NULL,color="pink",xintercept=48803, na.rm = FALSE, show.legend = TRUE)+
  #set the end point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",xintercept=56805, na.rm = FALSE, show.legend = TRUE)+
  geom_line()+
  ggtitle("Clinical 1")

#jpeg(filename = "Myplot.jpg", pointsize =12, quality = 200, bg = "white", res = NA, restoreConsole = TRUE)
multiplot(s1,s2, rows=2)

dev.off()