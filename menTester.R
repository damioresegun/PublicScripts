rm(list=ls())
graphics.off()
##Get libraries
#install.packages("ggplot2")
library(ggplot2)
setwd("~/")
source("multiplot.R")
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#################################### Coverage plot for HB3 ##################################
#isoPath<-args[1]
#isoName<-args[2]
chromo<-"Pf3D7_12_v3"
startPos<-46788
stopPos<-56805

midPosUp<-startPos+1200
midPosDown<-midPosUp+1000
print(midPosDown)

pdf(paste0("Test.pdf"))
# Set working directory
setwd("~/Testing/Mapping/HB3_genome")

# Read in the coverage file
ref1_data<-read.csv("HB3.cov",sep="\t",header=T)

# Extract features of the file which fall within the range of Chromosome 12
ref1_12_var<-ref1_data[Reduce(intersect,list(which(ref1_data[,1]==chromo),which(ref1_data[,2]>=30000),which(ref1_data[,2]<=75000))),]
# Extract sliding window before the start position
ref1_prior<-ref1_12_var[Reduce(intersect,list(which(ref1_12_var[,2]<midPosDown))),]
ref1_prior_last<-tail(ref1_prior,1)
# Extract the sliding window after the stop position
ref1_post<-ref1_12_var[Reduce(intersect,list(which(ref1_12_var[,2]>stopPos))),]
ref1_post_last<-head(ref1_post,1)
# Calculate the mean coverage of the entire chromosome 12
ref1_mean<-mean(ref1_data[which(ref1_data[,1]==chromo),5])
#### Calculate mean coverage of the gene of interest
ref1_gene<-ref1_12_var[Reduce(intersect,list(which(ref1_12_var[,2]>=midPosDown),which(ref1_12_var[,2]<=stopPos))),]
ref1_row<-nrow(ref1_gene)
ref1col<-colnames(ref1_gene[5])
ref1_covSum<-sum(ref1_gene[ref1col])
ref1_geneMean<-(ref1_covSum/ref1_row)

estimateCov<-round(ref1_geneMean/ref1_mean, digits=2)
positioner<-(ref1_geneMean+ref1_mean)/2
xPositioner<-(startPos+stopPos)/2

# Set the names of the axes
Positions=ref1_12_var[,2]+250
Coverages=ref1_12_var[,5]

# Plot the coverage of the range
s1<-ggplot(ref1_12_var,aes(x=Positions,y=Coverages))+
  geom_area(fill="purple",alpha=.2)+
  #use the mean coverage as a threshold in red
  geom_hline(mapping = NULL, data=NULL,color="red",size=1,yintercept=ref1_mean, na.rm = FALSE,show.legend = TRUE)+
  # Use mean gene coverage as a reference
  geom_segment(mapping = NULL, data=NULL,color="black",linetype="dotted",x=midPosDown-2000,xend=stopPos+2000,y=ref1_geneMean,yend=ref1_geneMean, na.rm = FALSE,show.legend = TRUE)+
  # set the starting point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",size=1.5,xintercept=startPos, na.rm = FALSE, show.legend = TRUE)+
  # set the starting point of the inter-exon intron
  geom_vline(mapping = NULL, data =NULL,color="green",size=1,xintercept=midPosUp, na.rm = FALSE, show.legend = TRUE)+
  # set the end point of the intron
  geom_vline(mapping = NULL, data =NULL,color="green",size=1,xintercept=midPosDown, na.rm = FALSE, show.legend = TRUE)+
  #set the end point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",size=1.5,xintercept=stopPos, na.rm = FALSE, show.legend = TRUE)+
  geom_text(x=xPositioner,y=positioner,label=estimateCov)+
  geom_line()+
  ggtitle("HB3")

#################################### Coverage plot for 3D7 ##################################

setwd("~/Testing/Mapping/3D7_genome")

# Read in the coverage file
ref2_data<-read.csv("3D7.cov",sep="\t",header=T)

# Extract features of the file which fall within the range of Chromosome 12
ref2_12_var<-ref2_data[Reduce(intersect,list(which(ref2_data[,1]==chromo),which(ref2_data[,2]>=30000),which(ref2_data[,2]<=75000))),]

# Extract sliding window before the start position
ref2_prior<-ref2_12_var[Reduce(intersect,list(which(ref2_12_var[,2]<midPosDown))),]
ref2_prior_last<-tail(ref2_prior,1)
# Extract the sliding window after the stop position
ref2_post<-ref2_12_var[Reduce(intersect,list(which(ref2_12_var[,2]>stopPos))),]
ref2_post_last<-head(ref2_post,1)
# Calculate the mean coverage of the entire chromosome 12
ref2_mean<-mean(ref2_data[which(ref2_data[,1]==chromo),5])
#### Calculate mean coverage of the gene of interest
ref2_gene<-ref2_12_var[Reduce(intersect,list(which(ref2_12_var[,2]>=midPosDown),which(ref2_12_var[,2]<=stopPos))),]
ref2_row<-nrow(ref2_gene)
ref2col<-colnames(ref2_gene[5])
ref2_covSum<-sum(ref2_gene[ref2col])
ref2_geneMean<-(ref2_covSum/ref2_row)

estimateCov<-round(ref2_geneMean/ref2_mean, digits=2)
positioner<-(ref2_geneMean+ref2_mean)/2
xPositioner<-(startPos+stopPos)/2

print(ref2col)

# Set the names of the axes
PositionX=ref2_12_var[,2]+250
CoverageX=ref2_12_var[,5]

# Plot the coverage of the range
s2<-ggplot(ref2_12_var,aes(x=PositionX,y=CoverageX))+
  geom_area(fill="grey",alpha=.2)+
  #use the mean coverage as a threshold in red
  geom_hline(mapping = NULL, data=NULL,color="red",size=1, yintercept=ref2_mean, na.rm = FALSE,show.legend = TRUE)+
  # Use mean gene coverage as a reference
  geom_segment(mapping = NULL, data=NULL,color="black",linetype="dotted",x=midPosDown-2000,xend=stopPos+2000,y=ref2_geneMean,yend=ref2_geneMean, na.rm = FALSE,show.legend = TRUE)+
  # set the starting point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",size=1.5,xintercept=startPos, na.rm = FALSE, show.legend = TRUE)+
  # set the starting point of the inter-exon intron
  geom_vline(mapping = NULL, data =NULL,color="forest green",size=1,xintercept=midPosUp, na.rm = FALSE, show.legend = TRUE)+
  # set the end point of the intron
  geom_vline(mapping = NULL, data =NULL,color="forest green",size=1,xintercept=midPosDown, na.rm = FALSE, show.legend = TRUE)+
  #set the end point of the gene
  geom_vline(mapping = NULL, data =NULL,color="blue",size=1.5,xintercept=stopPos, na.rm = FALSE, show.legend = TRUE)+
  geom_text(x=xPositioner,y=positioner,label=estimateCov)+
  geom_line()+
  ggtitle("3D7")
multiplot(s1,s2)
dev.off()

# Write out to csv 
setwd(Sys.getenv("HOME"))
ans<-matrix(, ncol = 2, nrow = 0)
ans<-rbind(ans,c("HB3",ref1_geneMean))
write.table(ans, file = "Total_Estimated_Copy_Number.csv",append=T,row.names=FALSE, na="",col.names=FALSE, sep=",")
