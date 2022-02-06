### first clear the workspace
rm(list = ls())
# install and load necessary packages
#BiocManager::install("karyoploteR") # uncomment this line to install
#BiocManager::install("rtracklayer") # uncomment this line to install
#BiocManager::install("formattable")
library(karyoploteR)
library(rtracklayer)
library(ggplotify)
library(cowplot)
library(formattable)

# load in your gff file
pknhGFF <- "Pknowlesi_pseudo.out.gff3"
sks047GFF <- "sks047_pseudo.out.gff3"
sks048GFF <- "sks048_pseudo.out.gff3"
cultGFF <- "Cultured_pseudo.out.gff3"
# read the first 'n' lines to get chromosome length information
pknhHEAD <- readLines(pknhGFF, n=36)
sks047HEAD <- readLines(sks047GFF, n=16)
sks048HEAD <- readLines(sks048GFF, n=16)
cultHEAD <- readLines(cultGFF, n=16)
# show the chromosome information
pknhHEAD
sks047HEAD
sks048HEAD
cultHEAD
# specifically get the lines that contain the chromosome length information
chroInfo047 <- sks047HEAD[grepl(sks047HEAD, pattern = "sequence-region   PKCLINC047")] # this will be hard to automate
chroInfo048 <- sks048HEAD[grepl(sks048HEAD, pattern = "sequence-region   PKCLINC048")]
chroInfoPKNH <- pknhHEAD[grepl(pknhHEAD, pattern = "sequence-region   PKNH")]
chroInfoCult <- cultHEAD[grepl(cultHEAD, pattern = "sequence-region   PKA1H1_STAND")]
# split this information and place them into a dataframe
ff047 <- data.frame(do.call(rbind, strsplit(chroInfo047, split = " ")))
ff048 <- data.frame(do.call(rbind, strsplit(chroInfo048, split = " ")))
ffPKNH <- data.frame(do.call(rbind, strsplit(chroInfoPKNH, split = " ")))
ffCULT <- data.frame(do.call(rbind, strsplit(chroInfoCult, split = " ")))
# Only place the chromosome name, start and stop positions into the dataframe
ff047[,5] <- as.numeric(as.character(ff047[,5]))
ff047[,6] <- as.numeric(as.character(ff047[,6]))
sks047 <- toGRanges(ff047[,c(4,5,6)])
ff048[,5] <- as.numeric(as.character(ff048[,5]))
ff048[,6] <- as.numeric(as.character(ff048[,6]))
sks048 <- toGRanges(ff048[,c(4,5,6)])
ffPKNH[,5] <- as.numeric(as.character(ffPKNH[,5]))
ffPKNH[,6] <- as.numeric(as.character(ffPKNH[,6]))
PKNHy <- toGRanges(ffPKNH[,c(4,5,6)])
ffCULT[,5] <- as.numeric(as.character(ffCULT[,5]))
ffCULT[,6] <- as.numeric(as.character(ffCULT[,6]))
CULT <- toGRanges(ffCULT[,c(4,5,6)])
# create features just for the gene  information
Feat047 <- import(sks047GFF)
Feat048 <- import(sks048GFF)
FeatPKNH <- import(pknhGFF)
FeatCULT <- import(cultGFF)
# show a table of the feature types present in the data
table(Feat047$type)
table(Feat048$type)
table(FeatPKNH$type)
table(FeatCULT$type)
# extract just the genes
Genes047 <- Feat047[Feat047$type=="gene"]
Genes048 <- Feat048[Feat048$type=="gene"]
GenesPKNH <- FeatPKNH[FeatPKNH$type=="gene"]
Genes0CULT <- FeatCULT[FeatCULT$type=="gene"]
Genes047
Genes048
GenesPKNH
Genes0CULT

seg.name <- paste("", 0:14, sep = "")
seg.name

#PKNH gene density plot
PKNHDens <- function() {
  #png("PKNH_genedensity.png", width = 1600, height = 800)
  pp <- getDefaultPlotParams(plot.type = 4)
  pp$data1inmargin <- 0
  pp$bottommargin <- 20
  kp <- plotKaryotype(genome = PKNHy, plot.type = 4, ideogram.plotter = NULL, labels.plotter = NULL, plot.params = pp)
  kpAddCytobandsAsLine(kp)
  kpAddChromosomeNames(kp, chr.names = seg.name, srt=0, cex=0.5)
  kp <- kpPlotDensity(kp, GenesPKNH, window.size = 1e05, col="sienna1")
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density,numticks = 8, cex=0.5)
  kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density)
  kpAddMainTitle(kp, "PKNH", cex=0.7)
  #dev.off()
}
CULTDens <- function() {
  #png("Cultured_genedensity.png", width = 1600, height = 800)
  pp <- getDefaultPlotParams(plot.type = 4)
  pp$data1inmargin <- 0
  pp$bottommargin <- 20
  kp <- plotKaryotype(genome = CULT, plot.type = 4, ideogram.plotter = NULL, labels.plotter = NULL, plot.params = pp)
  kpAddCytobandsAsLine(kp)
  kpAddChromosomeNames(kp, chr.names = seg.name, srt=0, cex=0.5)
  kp <- kpPlotDensity(kp, Genes0CULT, window.size = 1e05, col="seagreen")
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density,numticks = 8, cex=0.5)
  kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density)
  kpAddMainTitle(kp, "StAPkA1H1", cex=0.7)
  #dev.off()
}
#sks047 gene density plot
sks047Dens <- function() {
  #png("sks047_genedensity.png", width = 1600, height = 800)
  pp <- getDefaultPlotParams(plot.type = 4)
  pp$data1inmargin <- 0
  pp$bottommargin <- 20
  kp <- plotKaryotype(genome = sks047, plot.type = 4, ideogram.plotter = NULL, labels.plotter = NULL, plot.params = pp)
  kpAddCytobandsAsLine(kp)
  kpAddChromosomeNames(kp, chr.names = seg.name, srt=0, cex=0.5)
  kp <- kpPlotDensity(kp, Genes047, window.size = 1e05, col="skyblue")
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density,numticks = 8, cex=0.5)
  kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density)
  kpAddMainTitle(kp, "sks047", cex=0.7)
  #dev.off()
  }
#sks048 gene density plot
sks048Dens <- function() {
  #png("sks048_genedensity.png", width = 1600, height = 800)
  pp <- getDefaultPlotParams(plot.type = 4)
  pp$data1inmargin <- 0
  pp$bottommargin <- 20
  kp <- plotKaryotype(genome = sks048, plot.type = 4, ideogram.plotter = NULL, labels.plotter = NULL, plot.params = pp)
  kpAddCytobandsAsLine(kp)
  kpAddChromosomeNames(kp, chr.names = seg.name, srt=0, cex=0.5)
  kp <- kpPlotDensity(kp, Genes048, window.size = 1e05, col="coral2")
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density,numticks = 8, cex=0.5)
  kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density)
  kpAddMainTitle(kp, "sks048", cex=0.7)
  #dev.off()
}

#png("Isolate_Gene_Density.png", width = 2400, height = 1200)
tiff("Isolate_Gene_Density.tif", units="px", width=6000, height=3000,res = 600)
p1 <- as.ggplot(expression(PKNHDens()))
p2 <- as.ggplot(expression(CULTDens()))
p3 <- as.ggplot(expression(sks047Dens()))
p4 <- as.ggplot(expression(sks048Dens()))
plot_grid(p1, p2, p3, p4, ncol=2)
dev.off()
