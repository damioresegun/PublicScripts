#### Script to carry out visualisation of genomes after annotation by Companion
## Currently not automated but will attempt to implement this
## Requires a GFF file that you know the number of rows used BEFORE starting
## to state the gene information.
## This means, how many lines are used to described the chromosome information
## e.g. after "##gff-version 3", there will be multiple lines that look like this "##sequence-region   PKCLINC047_00 1 505490"
## how many lines like these are there BEFORE the first line that contains the information for each gene or contig/pesudogene
## that is what is needed here
## Input gene gffs also have to have the correct header row as follows: 
## "chromosome,unk1,feature,start,end,unk2,strand,unk3,gene,short,Attributes"
## The "unk#" columns are not used so can be named anything else.
## The "gene" column will have the common name for the gene identified to that region e.g calmodulin-binding protein
## The "short" column will have the short form of that gene name e.g. CBP  or SICAvar

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
setwd("C:/Users/dro/Dropbox/Reports/Thesis/ImageProcessing/Multigene Location/Gene_counts")
###### File load #############
# load in your gff file
pknhGFF <- "Pknowlesi_pseudo.out.gff3"
sks047GFF <- "sks047_pseudo.out.gff3"
sks048GFF <- "sks048_pseudo.out.gff3"
cultGFF <- "StAPkA1H1_pseudo.out.gff3"
sks050GFF <- "sks050_pseudo.out.gff3"
sks058GFF <- "sks058_pseudo.out.gff3"
sks070GFF <- "sks070_pseudo.out.gff3"
sks074GFF <- "sks074_pseudo.out.gff3"
sks078GFF <- "sks078_pseudo.out.gff3"
sks125GFF <- "sks125_pseudo.out.gff3"
sks325GFF <- "sks325_pseudo.out.gff3"
sks331GFF <- "sks331_pseudo.out.gff3"
sks333GFF <- "sks333_pseudo.out.gff3"
sks339GFF <- "sks339_pseudo.out.gff3"
sks344GFF <- "sks344_pseudo.out.gff3"
# load in the gene files
sks047IN_MSPs <- "sks047_MSPs.gff3"
sks048IN_MSPs <- "sks048_MSPs.gff3"
CultIN_MSPs <- "StAPkA1H1_MSPs.gff3"
pknhIN_MSPs <- "PKNH_MSPs.gff3"
sks050IN_MSPs <- "sks050_MSPs.gff3"
sks058IN_MSPs <- "sks058_MSPs.gff3"
sks070IN_MSPs <- "sks070_MSPs.gff3"
sks074IN_MSPs <- "sks074_MSPs.gff3"
sks078IN_MSPs <- "sks078_MSPs.gff3"
sks125IN_MSPs <- "sks125_MSPs.gff3"
sks325IN_MSPs <- "sks325_MSPs.gff3"
sks331IN_MSPs <- "sks331_MSPs.gff3"
sks333IN_MSPs <- "sks333_MSPs.gff3"
sks339IN_MSPs <- "sks339_MSPs.gff3"
sks344IN_MSPs <- "sks344_MSPs.gff3"
##
sks047IN_SPOR <- "sks047_sporo.gff3"
sks048IN_SPOR <- "sks048_sporo.gff3"
CultIN_SPOR <- "StAPkA1H1_sporo.gff3"
pknhIN_SPOR <- "PKNH_sporo.gff3"
sks050IN_SPOR <- "sks050_sporo.gff3"
sks058IN_SPOR <- "sks058_sporo.gff3"
sks070IN_SPOR <- "sks070_sporo.gff3"
sks074IN_SPOR <- "sks074_sporo.gff3"
sks078IN_SPOR <- "sks078_sporo.gff3"
sks125IN_SPOR <- "sks125_sporo.gff3"
sks325IN_SPOR <- "sks325_sporo.gff3"
sks331IN_SPOR <- "sks331_sporo.gff3"
sks333IN_SPOR <- "sks333_sporo.gff3"
sks339IN_SPOR <- "sks339_sporo.gff3"
sks344IN_SPOR <- "sks344_sporo.gff3"
##
sks047IN_pknbp <- "sks047_pknbp.gff3"
sks048IN_pknbp <- "sks048_pknbp.gff3"
CultIN_pknbp <- "StAPkA1H1_pknbp.gff3"
pknhIN_pknbp <- "PKNH_pknbp.gff3"
sks050IN_pknbp <- "sks050_pknbp.gff3"
sks058IN_pknbp <- "sks058_pknbp.gff3"
sks070IN_pknbp <- "sks070_pknbp.gff3"
sks074IN_pknbp <- "sks074_pknbp.gff3"
sks078IN_pknbp <- "sks078_pknbp.gff3"
sks125IN_pknbp <- "sks125_pknbp.gff3"
sks325IN_pknbp <- "sks325_pknbp.gff3"
sks331IN_pknbp <- "sks331_pknbp.gff3"
sks333IN_pknbp <- "sks333_pknbp.gff3"
sks339IN_pknbp <- "sks339_pknbp.gff3"
sks344IN_pknbp <- "sks344_pknbp.gff3"
##
sks047IN_kahrp <- "sks047_KAHRP.gff3"
sks048IN_kahrp <- "sks048_KAHRP.gff3"
CultIN_kahrp <- "StAPkA1H1_KAHRP.gff3"
pknhIN_kahrp <- "PKNH_KAHRP.gff3"
sks050IN_kahrp <- "sks050_KAHRP.gff3"
sks058IN_kahrp <- "sks058_KAHRP.gff3"
sks070IN_kahrp <- "sks070_KAHRP.gff3"
sks074IN_kahrp <- "sks074_KAHRP.gff3"
sks078IN_kahrp <- "sks078_KAHRP.gff3"
sks125IN_kahrp <- "sks125_KAHRP.gff3"
sks325IN_kahrp <- "sks325_KAHRP.gff3"
sks331IN_kahrp <- "sks331_KAHRP.gff3"
sks333IN_kahrp <- "sks333_KAHRP.gff3"
sks339IN_kahrp <- "sks339_KAHRP.gff3"
sks344IN_kahrp <- "sks344_KAHRP.gff3"
##
sks047IN_DBP <- "sks047_DBP.gff3"
sks048IN_DBP <- "sks048_DBP.gff3"
pknhIN_DBP <- "PKNH_DBP.gff3"
CultIN_DBP <- "StAPkA1H1_DBP.gff3"
sks050IN_DBP <- "sks050_DBP.gff3"
sks058IN_DBP <- "sks058_DBP.gff3"
sks070IN_DBP <- "sks070_DBP.gff3"
sks074IN_DBP <- "sks074_DBP.gff3"
sks078IN_DBP <- "sks078_DBP.gff3"
sks125IN_DBP <- "sks125_DBP.gff3"
sks325IN_DBP <- "sks325_DBP.gff3"
sks331IN_DBP <- "sks331_DBP.gff3"
sks333IN_DBP <- "sks333_DBP.gff3"
sks339IN_DBP <- "sks339_DBP.gff3"
sks344IN_DBP <- "sks344_DBP.gff3"
##
sks047IN_CyAP <- "sks047_CyAP.gff3"
sks048IN_CyAP <- "sks048_CyAP.gff3"
pknhIN_CyAP <- "PKNH_CyAP.gff3"
CultIN_CyAP <- "StAPkA1H1_CyAP.gff3"
sks050IN_CyAP <- "sks050_CyAP.gff3"
sks058IN_CyAP <- "sks058_CyAP.gff3"
sks070IN_CyAP <- "sks070_CyAP.gff3"
sks074IN_CyAP <- "sks074_CyAP.gff3"
sks078IN_CyAP <- "sks078_CyAP.gff3"
sks125IN_CyAP <- "sks125_CyAP.gff3"
sks325IN_CyAP <- "sks325_CyAP.gff3"
sks331IN_CyAP <- "sks331_CyAP.gff3"
sks333IN_CyAP <- "sks333_CyAP.gff3"
sks339IN_CyAP <- "sks339_CyAP.gff3"
sks344IN_CyAP <- "sks344_CyAP.gff3"
##
sks047IN_ETRAMP <- "sks047_ETRAMP.gff3"
sks048IN_ETRAMP <- "sks048_ETRAMP.gff3"
pknhIN_ETRAMP <- "PKNH_ETRAMP.gff3"
CultIN_ETRAMP <- "StAPkA1H1_ETRAMP.gff3"
sks050IN_ETRAMP <- "sks050_ETRAMP.gff3"
sks058IN_ETRAMP <- "sks058_ETRAMP.gff3"
sks070IN_ETRAMP <- "sks070_ETRAMP.gff3"
sks074IN_ETRAMP <- "sks074_ETRAMP.gff3"
sks078IN_ETRAMP <- "sks078_ETRAMP.gff3"
sks125IN_ETRAMP <- "sks125_ETRAMP.gff3"
sks325IN_ETRAMP <- "sks325_ETRAMP.gff3"
sks331IN_ETRAMP <- "sks331_ETRAMP.gff3"
sks333IN_ETRAMP <- "sks333_ETRAMP.gff3"
sks339IN_ETRAMP <- "sks339_ETRAMP.gff3"
sks344IN_ETRAMP <- "sks344_ETRAMP.gff3"
##
sks047IN_CIRC <- "sks047_circum.gff3"
sks048IN_CIRC <- "sks048_circum.gff3"
CultIN_CIRC <- "StAPkA1H1_circum.gff3"
pknhIN_CIRC <- "PKNH_circum.gff3"
sks050IN_CIRC <- "sks050_circum.gff3"
sks058IN_CIRC <- "sks058_circum.gff3"
sks070IN_CIRC <- "sks070_circum.gff3"
sks074IN_CIRC <- "sks074_circum.gff3"
sks078IN_CIRC <- "sks078_circum.gff3"
sks125IN_CIRC <- "sks125_circum.gff3"
sks325IN_CIRC <- "sks325_circum.gff3"
sks331IN_CIRC <- "sks331_circum.gff3"
sks333IN_CIRC <- "sks333_circum.gff3"
sks339IN_CIRC <- "sks339_circum.gff3"
sks344IN_CIRC <- "sks344_circum.gff3"
##
pknhIN_ABC <- "PKNH_ABC.gff3"
CultIN_ABC <- "StAPkA1H1_ABC.gff3"
sks047IN_ABC <- "sks047_ABC.gff3"
sks048IN_ABC <- "sks048_ABC.gff3"
sks050IN_ABC <- "sks050_ABC.gff3"
sks058IN_ABC <- "sks058_ABC.gff3"
sks070IN_ABC <- "sks070_ABC.gff3"
sks074IN_ABC <- "sks074_ABC.gff3"
sks078IN_ABC <- "sks078_ABC.gff3"
sks125IN_ABC <- "sks125_ABC.gff3"
sks325IN_ABC <- "sks325_ABC.gff3"
sks331IN_ABC <- "sks331_ABC.gff3"
sks333IN_ABC <- "sks333_ABC.gff3"
sks339IN_ABC <- "sks339_ABC.gff3"
sks344IN_ABC <- "sks344_ABC.gff3"
######## read in lines ##########
# read the first 'n' lines to get chromosome length information
pknhHEAD <- readLines(pknhGFF, n=36)
sks047HEAD <- readLines(sks047GFF, n=16)
sks048HEAD <- readLines(sks048GFF, n=16)
cultHEAD <- readLines(cultGFF, n=16)
sks050HEAD <- readLines(sks050GFF, n=16)
sks058HEAD <- readLines(sks058GFF, n=16)
sks070HEAD <- readLines(sks070GFF, n=16)
sks074HEAD <- readLines(sks074GFF, n=16)
sks078HEAD <- readLines(sks078GFF, n=16)
sks125HEAD <- readLines(sks125GFF, n=16)
sks325HEAD <- readLines(sks325GFF, n=16)
sks331HEAD <- readLines(sks331GFF, n=16)
sks333HEAD <- readLines(sks333GFF, n=16)
sks339HEAD <- readLines(sks339GFF, n=16)
sks344HEAD <- readLines(sks344GFF, n=16)
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
chroInfo050 <- sks050HEAD[grepl(sks050HEAD, pattern = "sequence-region   PKCLINC050")]
chroInfo058 <- sks058HEAD[grepl(sks058HEAD, pattern = "sequence-region   PKCLINC058")]
chroInfo070 <- sks070HEAD[grepl(sks070HEAD, pattern = "sequence-region   PKCLINC070")]
chroInfo074 <- sks074HEAD[grepl(sks074HEAD, pattern = "sequence-region   PKCLINC074")]
chroInfo078 <- sks078HEAD[grepl(sks078HEAD, pattern = "sequence-region   PKCLINC078")]
chroInfo125 <- sks125HEAD[grepl(sks125HEAD, pattern = "sequence-region   PKCLINC125")]
chroInfo325 <- sks325HEAD[grepl(sks325HEAD, pattern = "sequence-region   PKCLINC325")]
chroInfo331 <- sks331HEAD[grepl(sks331HEAD, pattern = "sequence-region   PKCLINC331")]
chroInfo333 <- sks333HEAD[grepl(sks333HEAD, pattern = "sequence-region   PKCLINC333")]
chroInfo339 <- sks339HEAD[grepl(sks339HEAD, pattern = "sequence-region   PKCLINC339")]
chroInfo344 <- sks344HEAD[grepl(sks344HEAD, pattern = "sequence-region   PKCLINC344")]
# split this information and place them into a dataframe
ff047 <- data.frame(do.call(rbind, strsplit(chroInfo047, split = " ")))
ff048 <- data.frame(do.call(rbind, strsplit(chroInfo048, split = " ")))
ffPKNH <- data.frame(do.call(rbind, strsplit(chroInfoPKNH, split = " ")))
ffCULT <- data.frame(do.call(rbind, strsplit(chroInfoCult, split = " ")))
ff050 <- data.frame(do.call(rbind, strsplit(chroInfo050, split = " ")))
ff058 <- data.frame(do.call(rbind, strsplit(chroInfo058, split = " ")))
ff070 <- data.frame(do.call(rbind, strsplit(chroInfo070, split = " ")))
ff074 <- data.frame(do.call(rbind, strsplit(chroInfo074, split = " ")))
ff078 <- data.frame(do.call(rbind, strsplit(chroInfo078, split = " ")))
ff125 <- data.frame(do.call(rbind, strsplit(chroInfo125, split = " ")))
ff325 <- data.frame(do.call(rbind, strsplit(chroInfo325, split = " ")))
ff331 <- data.frame(do.call(rbind, strsplit(chroInfo331, split = " ")))
ff333 <- data.frame(do.call(rbind, strsplit(chroInfo333, split = " ")))
ff339 <- data.frame(do.call(rbind, strsplit(chroInfo339, split = " ")))
ff344 <- data.frame(do.call(rbind, strsplit(chroInfo344, split = " ")))
# Only place the chromosome name, start and stop positions into the dataframe
ff047[,5] <- as.numeric(as.character(ff047[,5]))
ff047[,6] <- as.numeric(as.character(ff047[,6]))
sks047 <- toGRanges(ff047[,c(4,5,6)])
##
ff048[,5] <- as.numeric(as.character(ff048[,5]))
ff048[,6] <- as.numeric(as.character(ff048[,6]))
sks048 <- toGRanges(ff048[,c(4,5,6)])
##
ffPKNH[,5] <- as.numeric(as.character(ffPKNH[,5]))
ffPKNH[,6] <- as.numeric(as.character(ffPKNH[,6]))
PKNHy <- toGRanges(ffPKNH[,c(4,5,6)])
##
ffCULT[,5] <- as.numeric(as.character(ffCULT[,5]))
ffCULT[,6] <- as.numeric(as.character(ffCULT[,6]))
CULT <- toGRanges(ffCULT[,c(4,5,6)])
##
ff050[,5] <- as.numeric(as.character(ff050[,5]))
ff050[,6] <- as.numeric(as.character(ff050[,6]))
sks050 <- toGRanges(ff050[,c(4,5,6)])
##
ff058[,5] <- as.numeric(as.character(ff058[,5]))
ff058[,6] <- as.numeric(as.character(ff058[,6]))
sks058 <- toGRanges(ff058[,c(4,5,6)])
##
ff070[,5] <- as.numeric(as.character(ff070[,5]))
ff070[,6] <- as.numeric(as.character(ff070[,6]))
sks070 <- toGRanges(ff070[,c(4,5,6)])
##
ff074[,5] <- as.numeric(as.character(ff074[,5]))
ff074[,6] <- as.numeric(as.character(ff074[,6]))
sks074 <- toGRanges(ff074[,c(4,5,6)])
##
ff078[,5] <- as.numeric(as.character(ff078[,5]))
ff078[,6] <- as.numeric(as.character(ff078[,6]))
sks078 <- toGRanges(ff078[,c(4,5,6)])
##
ff125[,5] <- as.numeric(as.character(ff125[,5]))
ff125[,6] <- as.numeric(as.character(ff125[,6]))
sks125 <- toGRanges(ff125[,c(4,5,6)])
##
ff325[,5] <- as.numeric(as.character(ff325[,5]))
ff325[,6] <- as.numeric(as.character(ff325[,6]))
sks325 <- toGRanges(ff325[,c(4,5,6)])
##
ff331[,5] <- as.numeric(as.character(ff331[,5]))
ff331[,6] <- as.numeric(as.character(ff331[,6]))
sks331 <- toGRanges(ff331[,c(4,5,6)])
##
ff333[,5] <- as.numeric(as.character(ff333[,5]))
ff333[,6] <- as.numeric(as.character(ff333[,6]))
sks333 <- toGRanges(ff333[,c(4,5,6)])
##
ff339[,5] <- as.numeric(as.character(ff339[,5]))
ff339[,6] <- as.numeric(as.character(ff339[,6]))
sks339 <- toGRanges(ff339[,c(4,5,6)])
##
ff344[,5] <- as.numeric(as.character(ff344[,5]))
ff344[,6] <- as.numeric(as.character(ff344[,6]))
sks344 <- toGRanges(ff344[,c(4,5,6)])
##
########## Create features ##########
# create features just for the gene  information
Feat047 <- import(sks047GFF)
Feat048 <- import(sks048GFF)
FeatPKNH <- import(pknhGFF)
FeatCULT <- import(cultGFF)
Feat050 <- import(sks050GFF)
Feat058 <- import(sks058GFF)
Feat070 <- import(sks070GFF)
Feat074 <- import(sks074GFF)
Feat078 <- import(sks078GFF)
Feat125 <- import(sks125GFF)
Feat325 <- import(sks325GFF)
Feat331 <- import(sks331GFF)
Feat333 <- import(sks333GFF)
Feat339 <- import(sks339GFF)
Feat344 <- import(sks344GFF)
# show a table of the feature types present in the data
table(Feat047$type)
table(Feat048$type)
table(FeatPKNH$type)
table(FeatCULT$type)
table(Feat050$type)
table(Feat058$type)
table(Feat070$type)
table(Feat074$type)
table(Feat078$type)
table(Feat125$type)
table(Feat325$type)
table(Feat331$type)
table(Feat333$type)
table(Feat339$type)
table(Feat344$type)
########## Extract Genes ##########
# extract just the genes
Genes047 <- Feat047[Feat047$type=="gene"]
Genes048 <- Feat048[Feat048$type=="gene"]
GenesPKNH <- FeatPKNH[FeatPKNH$type=="gene"]
Genes0CULT <- FeatCULT[FeatCULT$type=="gene"]
Genes050 <- Feat050[Feat050$type=="gene"]
Genes058 <- Feat058[Feat058$type=="gene"]
Genes070 <- Feat070[Feat070$type=="gene"]
Genes074 <- Feat074[Feat074$type=="gene"]
Genes078 <- Feat078[Feat078$type=="gene"]
Genes125 <- Feat125[Feat125$type=="gene"]
Genes325 <- Feat325[Feat325$type=="gene"]
Genes331 <- Feat331[Feat331$type=="gene"]
Genes333 <- Feat333[Feat333$type=="gene"]
Genes339 <- Feat339[Feat339$type=="gene"]
Genes344 <- Feat344[Feat344$type=="gene"]
Genes047
Genes048
GenesPKNH
Genes0CULT
########################################################################################################
###### add in the MSPgenes #######
msp047IN <- read.delim(sks047IN_MSPs)
msp047 <- toGRanges(data.frame(chromosome = c(msp047IN$chromosome), start=c(msp047IN$start),
                               end=c(msp047IN$end), gene = c(msp047IN$short), strand=c(msp047IN$strand)))
msp047 <- sort(msp047)
##
msp048IN <- read.delim(sks048IN_MSPs)
msp048 <- toGRanges(data.frame(chromosome = c(msp048IN$chromosome), start=c(msp048IN$start),
                               end=c(msp048IN$end), gene = c(msp048IN$short), strand=c(msp048IN$strand)))
msp048 <- sort(msp048)
##
mspCultIN <- read.delim(CultIN_MSPs)
mspCult <- toGRanges(data.frame(chromosome = c(mspCultIN$chromosome), start=c(mspCultIN$start),
                                end=c(mspCultIN$end), gene = c(mspCultIN$short), strand=c(mspCultIN$strand)))
mspCult <- sort(mspCult)
##
mspPKnhIN <- read.delim(pknhIN_MSPs)
mspPKNH <- toGRanges(data.frame(chromosome = c(mspPKnhIN$chromosome), start=c(mspPKnhIN$start),
                                end=c(mspPKnhIN$end), gene = c(mspPKnhIN$short), strand=c(mspPKnhIN$strand)))
mspPKNH <- sort(mspPKNH)
##
msp050IN <- read.delim(sks050IN_MSPs)
msp050 <- toGRanges(data.frame(chromosome = c(msp050IN$chromosome), start=c(msp050IN$start),
                               end=c(msp050IN$end), gene = c(msp050IN$short), strand=c(msp050IN$strand)))
msp050 <- sort(msp050)
##
##
msp058IN <- read.delim(sks058IN_MSPs)
msp058 <- toGRanges(data.frame(chromosome = c(msp058IN$chromosome), start=c(msp058IN$start),
                               end=c(msp058IN$end), gene = c(msp058IN$short), strand=c(msp058IN$strand)))
msp058 <- sort(msp058)
##
##
msp070IN <- read.delim(sks070IN_MSPs)
msp070 <- toGRanges(data.frame(chromosome = c(msp070IN$chromosome), start=c(msp070IN$start),
                               end=c(msp070IN$end), gene = c(msp070IN$short), strand=c(msp070IN$strand)))
msp070 <- sort(msp070)
##
##
msp074IN <- read.delim(sks074IN_MSPs)
msp074 <- toGRanges(data.frame(chromosome = c(msp074IN$chromosome), start=c(msp074IN$start),
                               end=c(msp074IN$end), gene = c(msp074IN$short), strand=c(msp074IN$strand)))
msp074 <- sort(msp074)
##
##
msp078IN <- read.delim(sks078IN_MSPs)
msp078 <- toGRanges(data.frame(chromosome = c(msp078IN$chromosome), start=c(msp078IN$start),
                               end=c(msp078IN$end), gene = c(msp078IN$short), strand=c(msp078IN$strand)))
msp078 <- sort(msp078)
##
##
msp125IN <- read.delim(sks125IN_MSPs)
msp125 <- toGRanges(data.frame(chromosome = c(msp125IN$chromosome), start=c(msp125IN$start),
                               end=c(msp125IN$end), gene = c(msp125IN$short), strand=c(msp125IN$strand)))
msp125 <- sort(msp125)
##
##
msp325IN <- read.delim(sks325IN_MSPs)
msp325 <- toGRanges(data.frame(chromosome = c(msp325IN$chromosome), start=c(msp325IN$start),
                               end=c(msp325IN$end), gene = c(msp325IN$short), strand=c(msp325IN$strand)))
msp325 <- sort(msp325)
##
##
msp331IN <- read.delim(sks331IN_MSPs)
msp331 <- toGRanges(data.frame(chromosome = c(msp331IN$chromosome), start=c(msp331IN$start),
                               end=c(msp331IN$end), gene = c(msp331IN$short), strand=c(msp331IN$strand)))
msp331 <- sort(msp331)
##
##
msp333IN <- read.delim(sks333IN_MSPs)
msp333 <- toGRanges(data.frame(chromosome = c(msp333IN$chromosome), start=c(msp333IN$start),
                               end=c(msp333IN$end), gene = c(msp333IN$short), strand=c(msp333IN$strand)))
msp333 <- sort(msp333)
##
##
msp339IN <- read.delim(sks339IN_MSPs)
msp339 <- toGRanges(data.frame(chromosome = c(msp339IN$chromosome), start=c(msp339IN$start),
                               end=c(msp339IN$end), gene = c(msp339IN$short), strand=c(msp339IN$strand)))
msp339 <- sort(msp339)
##
##
msp344IN <- read.delim(sks344IN_MSPs)
msp344 <- toGRanges(data.frame(chromosome = c(msp344IN$chromosome), start=c(msp344IN$start),
                               end=c(msp344IN$end), gene = c(msp344IN$short), strand=c(msp344IN$strand)))
msp344 <- sort(msp344)
##

####### add in the  circumsporozoite genes ########
circ047IN <- read.delim(sks047IN_CIRC)
circ047 <- toGRanges(data.frame(chromosome = c(circ047IN$chromosome), start=c(circ047IN$start),
                                end=c(circ047IN$end), gene = c(circ047IN$short), strand=c(circ047IN$strand)))
circ047 <- sort(circ047)
##
circ048IN <- read.delim(sks048IN_CIRC)
circ048 <- toGRanges(data.frame(chromosome = c(circ048IN$chromosome), start=c(circ048IN$start),
                                end=c(circ048IN$end), gene = c(circ048IN$short), strand=c(circ048IN$strand)))
circ048 <- sort(circ048)
##
circCultIN <- read.delim(CultIN_CIRC)
circCult <- toGRanges(data.frame(chromosome = c(circCultIN$chromosome), start=c(circCultIN$start),
                                 end=c(circCultIN$end), gene = c(circCultIN$short), strand=c(circCultIN$strand)))
circCult <- sort(circCult)
##
circPKNHIN <- read.delim(pknhIN_CIRC)
circPKNH <- toGRanges(data.frame(chromosome = c(circPKNHIN$chromosome), start=c(circPKNHIN$start),
                                 end=c(circPKNHIN$end), gene = c(circPKNHIN$short), strand=c(circPKNHIN$strand)))
circPKNH <- sort(circPKNH)
##
circ050IN <- read.delim(sks050IN_CIRC)
circ050 <- toGRanges(data.frame(chromosome = c(circ050IN$chromosome), start=c(circ050IN$start),
                               end=c(circ050IN$end), gene = c(circ050IN$short), strand=c(circ050IN$strand)))
circ050 <- sort(circ050)
##
##
circ058IN <- read.delim(sks058IN_CIRC)
circ058 <- toGRanges(data.frame(chromosome = c(circ058IN$chromosome), start=c(circ058IN$start),
                               end=c(circ058IN$end), gene = c(circ058IN$short), strand=c(circ058IN$strand)))
circ058 <- sort(circ058)
##
##
circ070IN <- read.delim(sks070IN_CIRC)
circ070 <- toGRanges(data.frame(chromosome = c(circ070IN$chromosome), start=c(circ070IN$start),
                               end=c(circ070IN$end), gene = c(circ070IN$short), strand=c(circ070IN$strand)))
circ070 <- sort(circ070)
##
##
circ074IN <- read.delim(sks074IN_CIRC)
circ074 <- toGRanges(data.frame(chromosome = c(circ074IN$chromosome), start=c(circ074IN$start),
                               end=c(circ074IN$end), gene = c(circ074IN$short), strand=c(circ074IN$strand)))
circ074 <- sort(circ074)
##
##
circ078IN <- read.delim(sks078IN_CIRC)
circ078 <- toGRanges(data.frame(chromosome = c(circ078IN$chromosome), start=c(circ078IN$start),
                               end=c(circ078IN$end), gene = c(circ078IN$short), strand=c(circ078IN$strand)))
circ078 <- sort(circ078)
##
##
circ125IN <- read.delim(sks125IN_CIRC)
circ125 <- toGRanges(data.frame(chromosome = c(circ125IN$chromosome), start=c(circ125IN$start),
                               end=c(circ125IN$end), gene = c(circ125IN$short), strand=c(circ125IN$strand)))
circ125 <- sort(circ125)
##
##
circ325IN <- read.delim(sks325IN_CIRC)
circ325 <- toGRanges(data.frame(chromosome = c(circ325IN$chromosome), start=c(circ325IN$start),
                               end=c(circ325IN$end), gene = c(circ325IN$short), strand=c(circ325IN$strand)))
circ325 <- sort(circ325)
##
##
circ331IN <- read.delim(sks331IN_CIRC)
circ331 <- toGRanges(data.frame(chromosome = c(circ331IN$chromosome), start=c(circ331IN$start),
                               end=c(circ331IN$end), gene = c(circ331IN$short), strand=c(circ331IN$strand)))
circ331 <- sort(circ331)
##
##
circ333IN <- read.delim(sks333IN_CIRC)
circ333 <- toGRanges(data.frame(chromosome = c(circ333IN$chromosome), start=c(circ333IN$start),
                               end=c(circ333IN$end), gene = c(circ333IN$short), strand=c(circ333IN$strand)))
circ333 <- sort(circ333)
##
##
circ339IN <- read.delim(sks339IN_CIRC)
circ339 <- toGRanges(data.frame(chromosome = c(circ339IN$chromosome), start=c(circ339IN$start),
                               end=c(circ339IN$end), gene = c(circ339IN$short), strand=c(circ339IN$strand)))
circ339 <- sort(circ339)
##
##
circ344IN <- read.delim(sks344IN_CIRC)
circ344 <- toGRanges(data.frame(chromosome = c(circ344IN$chromosome), start=c(circ344IN$start),
                               end=c(circ344IN$end), gene = c(circ344IN$short), strand=c(circ344IN$strand)))
circ344 <- sort(circ344)
##
####### add in the  sporozoite invasion-associated protein genes ########
spor047IN <- read.delim(sks047IN_SPOR)
spor047 <- toGRanges(data.frame(chromosome = c(spor047IN$chromosome), start=c(spor047IN$start),
                                end=c(spor047IN$end), gene = c(spor047IN$short), strand=c(spor047IN$strand)))
spor047 <- sort(spor047)
##
spor048IN <- read.delim(sks048IN_SPOR)
spor048 <- toGRanges(data.frame(chromosome = c(spor048IN$chromosome), start=c(spor048IN$start),
                                end=c(spor048IN$end), gene = c(spor048IN$short), strand=c(spor048IN$strand)))
spor048 <- sort(spor048)
##
sporCultIN <- read.delim(CultIN_SPOR)
sporCult <- toGRanges(data.frame(chromosome = c(sporCultIN$chromosome), start=c(sporCultIN$start),
                                 end=c(sporCultIN$end), gene = c(sporCultIN$short), strand=c(sporCultIN$strand)))
sporCult <- sort(sporCult)
##
sporPKNHIN <- read.delim(pknhIN_SPOR)
sporPKNH <- toGRanges(data.frame(chromosome = c(sporPKNHIN$chromosome), start=c(sporPKNHIN$start),
                                 end=c(sporPKNHIN$end), gene = c(sporPKNHIN$short), strand=c(sporPKNHIN$strand)))
sporPKNH <- sort(sporPKNH)
##
spor050IN <- read.delim(sks050IN_SPOR)
spor050 <- toGRanges(data.frame(chromosome = c(spor050IN$chromosome), start=c(spor050IN$start),
                               end=c(spor050IN$end), gene = c(spor050IN$short), strand=c(spor050IN$strand)))
spor050 <- sort(spor050)
##
##
spor058IN <- read.delim(sks058IN_SPOR)
spor058 <- toGRanges(data.frame(chromosome = c(spor058IN$chromosome), start=c(spor058IN$start),
                               end=c(spor058IN$end), gene = c(spor058IN$short), strand=c(spor058IN$strand)))
spor058 <- sort(spor058)
##
##
spor070IN <- read.delim(sks070IN_SPOR)
spor070 <- toGRanges(data.frame(chromosome = c(spor070IN$chromosome), start=c(spor070IN$start),
                               end=c(spor070IN$end), gene = c(spor070IN$short), strand=c(spor070IN$strand)))
spor070 <- sort(spor070)
##
##
spor074IN <- read.delim(sks074IN_SPOR)
spor074 <- toGRanges(data.frame(chromosome = c(spor074IN$chromosome), start=c(spor074IN$start),
                               end=c(spor074IN$end), gene = c(spor074IN$short), strand=c(spor074IN$strand)))
spor074 <- sort(spor074)
##
##
spor078IN <- read.delim(sks078IN_SPOR)
spor078 <- toGRanges(data.frame(chromosome = c(spor078IN$chromosome), start=c(spor078IN$start),
                               end=c(spor078IN$end), gene = c(spor078IN$short), strand=c(spor078IN$strand)))
spor078 <- sort(spor078)
##
##
spor125IN <- read.delim(sks125IN_SPOR)
spor125 <- toGRanges(data.frame(chromosome = c(spor125IN$chromosome), start=c(spor125IN$start),
                               end=c(spor125IN$end), gene = c(spor125IN$short), strand=c(spor125IN$strand)))
spor125 <- sort(spor125)
##
##
spor325IN <- read.delim(sks325IN_SPOR)
spor325 <- toGRanges(data.frame(chromosome = c(spor325IN$chromosome), start=c(spor325IN$start),
                               end=c(spor325IN$end), gene = c(spor325IN$short), strand=c(spor325IN$strand)))
spor325 <- sort(spor325)
##
##
spor331IN <- read.delim(sks331IN_SPOR)
spor331 <- toGRanges(data.frame(chromosome = c(spor331IN$chromosome), start=c(spor331IN$start),
                               end=c(spor331IN$end), gene = c(spor331IN$short), strand=c(spor331IN$strand)))
spor331 <- sort(spor331)
##
##
spor333IN <- read.delim(sks333IN_SPOR)
spor333 <- toGRanges(data.frame(chromosome = c(spor333IN$chromosome), start=c(spor333IN$start),
                               end=c(spor333IN$end), gene = c(spor333IN$short), strand=c(spor333IN$strand)))
spor333 <- sort(spor333)
##
##
spor339IN <- read.delim(sks339IN_SPOR)
spor339 <- toGRanges(data.frame(chromosome = c(spor339IN$chromosome), start=c(spor339IN$start),
                               end=c(spor339IN$end), gene = c(spor339IN$short), strand=c(spor339IN$strand)))
spor339 <- sort(spor339)
##
##
spor344IN <- read.delim(sks344IN_SPOR)
spor344 <- toGRanges(data.frame(chromosome = c(spor344IN$chromosome), start=c(spor344IN$start),
                               end=c(spor344IN$end), gene = c(spor344IN$short), strand=c(spor344IN$strand)))
spor344 <- sort(spor344)
##
####### add in the  pknbp genes ########
pknb047IN <- read.delim(sks047IN_pknbp)
pknb047 <- toGRanges(data.frame(chromosome = c(pknb047IN$chromosome), start=c(pknb047IN$start),
                                end=c(pknb047IN$end), gene = c(pknb047IN$short), strand=c(pknb047IN$strand)))
pknb047 <- sort(pknb047)
##
pknb048IN <- read.delim(sks048IN_pknbp)
pknb048 <- toGRanges(data.frame(chromosome = c(pknb048IN$chromosome), start=c(pknb048IN$start),
                                end=c(pknb048IN$end), gene = c(pknb048IN$short), strand=c(pknb048IN$strand)))
pknb048 <- sort(pknb048)

##
pknbCultIN <- read.delim(CultIN_pknbp)
pknbCult <- toGRanges(data.frame(chromosome = c(pknbCultIN$chromosome), start=c(pknbCultIN$start),
                                 end=c(pknbCultIN$end), gene = c(pknbCultIN$short), strand=c(pknbCultIN$strand)))
pknbCult <- sort(pknbCult)
##
pknbPKNHIN <- read.delim(pknhIN_pknbp)
pknbPKNH <- toGRanges(data.frame(chromosome = c(pknbPKNHIN$chromosome), start=c(pknbPKNHIN$start),
                                 end=c(pknbPKNHIN$end), gene = c(pknbPKNHIN$short), strand=c(pknbPKNHIN$strand)))
pknbPKNH <- sort(pknbPKNH)
##
pknb050IN <- read.delim(sks050IN_pknbp)
pknb050 <- toGRanges(data.frame(chromosome = c(pknb050IN$chromosome), start=c(pknb050IN$start),
                               end=c(pknb050IN$end), gene = c(pknb050IN$short), strand=c(pknb050IN$strand)))
pknb050 <- sort(pknb050)
##
##
pknb058IN <- read.delim(sks058IN_pknbp)
pknb058 <- toGRanges(data.frame(chromosome = c(pknb058IN$chromosome), start=c(pknb058IN$start),
                               end=c(pknb058IN$end), gene = c(pknb058IN$short), strand=c(pknb058IN$strand)))
pknb058 <- sort(pknb058)
##
##
pknb070IN <- read.delim(sks070IN_pknbp)
pknb070 <- toGRanges(data.frame(chromosome = c(pknb070IN$chromosome), start=c(pknb070IN$start),
                               end=c(pknb070IN$end), gene = c(pknb070IN$short), strand=c(pknb070IN$strand)))
pknb070 <- sort(pknb070)
##
##
pknb074IN <- read.delim(sks074IN_pknbp)
pknb074 <- toGRanges(data.frame(chromosome = c(pknb074IN$chromosome), start=c(pknb074IN$start),
                               end=c(pknb074IN$end), gene = c(pknb074IN$short), strand=c(pknb074IN$strand)))
pknb074 <- sort(pknb074)
##
##
pknb078IN <- read.delim(sks078IN_pknbp)
pknb078 <- toGRanges(data.frame(chromosome = c(pknb078IN$chromosome), start=c(pknb078IN$start),
                               end=c(pknb078IN$end), gene = c(pknb078IN$short), strand=c(pknb078IN$strand)))
pknb078 <- sort(pknb078)
##
##
pknb125IN <- read.delim(sks125IN_pknbp)
pknb125 <- toGRanges(data.frame(chromosome = c(pknb125IN$chromosome), start=c(pknb125IN$start),
                               end=c(pknb125IN$end), gene = c(pknb125IN$short), strand=c(pknb125IN$strand)))
pknb125 <- sort(pknb125)
##
##
pknb325IN <- read.delim(sks325IN_pknbp)
pknb325 <- toGRanges(data.frame(chromosome = c(pknb325IN$chromosome), start=c(pknb325IN$start),
                               end=c(pknb325IN$end), gene = c(pknb325IN$short), strand=c(pknb325IN$strand)))
pknb325 <- sort(pknb325)
##
##
pknb331IN <- read.delim(sks331IN_pknbp)
pknb331 <- toGRanges(data.frame(chromosome = c(pknb331IN$chromosome), start=c(pknb331IN$start),
                               end=c(pknb331IN$end), gene = c(pknb331IN$short), strand=c(pknb331IN$strand)))
pknb331 <- sort(pknb331)
##
##
pknb333IN <- read.delim(sks333IN_pknbp)
pknb333 <- toGRanges(data.frame(chromosome = c(pknb333IN$chromosome), start=c(pknb333IN$start),
                               end=c(pknb333IN$end), gene = c(pknb333IN$short), strand=c(pknb333IN$strand)))
pknb333 <- sort(pknb333)
##
##
pknb339IN <- read.delim(sks339IN_pknbp)
pknb339 <- toGRanges(data.frame(chromosome = c(pknb339IN$chromosome), start=c(pknb339IN$start),
                               end=c(pknb339IN$end), gene = c(pknb339IN$short), strand=c(pknb339IN$strand)))
pknb339 <- sort(pknb339)
##
##
pknb344IN <- read.delim(sks344IN_pknbp)
pknb344 <- toGRanges(data.frame(chromosome = c(pknb344IN$chromosome), start=c(pknb344IN$start),
                               end=c(pknb344IN$end), gene = c(pknb344IN$short), strand=c(pknb344IN$strand)))
pknb344 <- sort(pknb344)
##
####### add in the knob-associated histidine rich protein genes ########
kahrp047IN <- read.delim(sks047IN_kahrp)
kahrp047 <- toGRanges(data.frame(chromosome = c(kahrp047IN$chromosome), start=c(kahrp047IN$start),
                                 end=c(kahrp047IN$end), gene = c(kahrp047IN$short), strand=c(kahrp047IN$strand)))
kahrp047 <- sort(kahrp047)
##
kahrp048IN <- read.delim(sks048IN_kahrp)
kahrp048 <- toGRanges(data.frame(chromosome = c(kahrp048IN$chromosome), start=c(kahrp048IN$start),
                                 end=c(kahrp048IN$end), gene = c(kahrp048IN$short), strand=c(kahrp048IN$strand)))
kahrp048 <- sort(kahrp048)
##
kahrpCultIN <- read.delim(CultIN_kahrp)
kahrpCult <- toGRanges(data.frame(chromosome = c(kahrpCultIN$chromosome), start=c(kahrpCultIN$start),
                                  end=c(kahrpCultIN$end), gene = c(kahrpCultIN$short), strand=c(kahrpCultIN$strand)))
kahrpCult <- sort(kahrpCult)
##
kahrpPKNHIN <- read.delim(pknhIN_kahrp)
kahrpPKNH <- toGRanges(data.frame(chromosome = c(kahrpPKNHIN$chromosome), start=c(kahrpPKNHIN$start),
                                  end=c(kahrpPKNHIN$end), gene = c(kahrpPKNHIN$short), strand=c(kahrpPKNHIN$strand)))
kahrpPKNH <- sort(kahrpPKNH)
##
kahrp050IN <- read.delim(sks050IN_kahrp)
kahrp050 <- toGRanges(data.frame(chromosome = c(kahrp050IN$chromosome), start=c(kahrp050IN$start),
                               end=c(kahrp050IN$end), gene = c(kahrp050IN$short), strand=c(kahrp050IN$strand)))
kahrp050 <- sort(kahrp050)
##
##
kahrp058IN <- read.delim(sks058IN_kahrp)
kahrp058 <- toGRanges(data.frame(chromosome = c(kahrp058IN$chromosome), start=c(kahrp058IN$start),
                               end=c(kahrp058IN$end), gene = c(kahrp058IN$short), strand=c(kahrp058IN$strand)))
kahrp058 <- sort(kahrp058)
##
##
kahrp070IN <- read.delim(sks070IN_kahrp)
kahrp070 <- toGRanges(data.frame(chromosome = c(kahrp070IN$chromosome), start=c(kahrp070IN$start),
                               end=c(kahrp070IN$end), gene = c(kahrp070IN$short), strand=c(kahrp070IN$strand)))
kahrp070 <- sort(kahrp070)
##
##
kahrp074IN <- read.delim(sks074IN_kahrp)
kahrp074 <- toGRanges(data.frame(chromosome = c(kahrp074IN$chromosome), start=c(kahrp074IN$start),
                               end=c(kahrp074IN$end), gene = c(kahrp074IN$short), strand=c(kahrp074IN$strand)))
kahrp074 <- sort(kahrp074)
##
##
kahrp078IN <- read.delim(sks078IN_kahrp)
kahrp078 <- toGRanges(data.frame(chromosome = c(kahrp078IN$chromosome), start=c(kahrp078IN$start),
                               end=c(kahrp078IN$end), gene = c(kahrp078IN$short), strand=c(kahrp078IN$strand)))
kahrp078 <- sort(kahrp078)
##
##
kahrp125IN <- read.delim(sks125IN_kahrp)
kahrp125 <- toGRanges(data.frame(chromosome = c(kahrp125IN$chromosome), start=c(kahrp125IN$start),
                               end=c(kahrp125IN$end), gene = c(kahrp125IN$short), strand=c(kahrp125IN$strand)))
kahrp125 <- sort(kahrp125)
##
##
kahrp325IN <- read.delim(sks325IN_kahrp)
kahrp325 <- toGRanges(data.frame(chromosome = c(kahrp325IN$chromosome), start=c(kahrp325IN$start),
                               end=c(kahrp325IN$end), gene = c(kahrp325IN$short), strand=c(kahrp325IN$strand)))
kahrp325 <- sort(kahrp325)
##
##
kahrp331IN <- read.delim(sks331IN_kahrp)
kahrp331 <- toGRanges(data.frame(chromosome = c(kahrp331IN$chromosome), start=c(kahrp331IN$start),
                               end=c(kahrp331IN$end), gene = c(kahrp331IN$short), strand=c(kahrp331IN$strand)))
kahrp331 <- sort(kahrp331)
##
##
kahrp333IN <- read.delim(sks333IN_kahrp)
kahrp333 <- toGRanges(data.frame(chromosome = c(kahrp333IN$chromosome), start=c(kahrp333IN$start),
                               end=c(kahrp333IN$end), gene = c(kahrp333IN$short), strand=c(kahrp333IN$strand)))
kahrp333 <- sort(kahrp333)
##
##
kahrp339IN <- read.delim(sks339IN_kahrp)
kahrp339 <- toGRanges(data.frame(chromosome = c(kahrp339IN$chromosome), start=c(kahrp339IN$start),
                               end=c(kahrp339IN$end), gene = c(kahrp339IN$short), strand=c(kahrp339IN$strand)))
kahrp339 <- sort(kahrp339)
##
##
kahrp344IN <- read.delim(sks344IN_kahrp)
kahrp344 <- toGRanges(data.frame(chromosome = c(kahrp344IN$chromosome), start=c(kahrp344IN$start),
                               end=c(kahrp344IN$end), gene = c(kahrp344IN$short), strand=c(kahrp344IN$strand)))
kahrp344 <- sort(kahrp344)
##
####### add in the CLAG genes ########
clag047IN <- read.delim(sks047IN_CyAP)
clag047 <- toGRanges(data.frame(chromosome = c(clag047IN$chromosome), start=c(clag047IN$start),
                                end=c(clag047IN$end), gene = c(clag047IN$short), strand=c(clag047IN$strand)))
clag047 <- sort(clag047)
##
clag048IN <- read.delim(sks048IN_CyAP)
clag048 <- toGRanges(data.frame(chromosome = c(clag048IN$chromosome), start=c(clag048IN$start),
                                end=c(clag048IN$end), gene = c(clag048IN$short), strand=c(clag048IN$strand)))
clag048 <- sort(clag048)
##
clagCultIN <- read.delim(CultIN_CyAP)
clagCult <- toGRanges(data.frame(chromosome = c(clagCultIN$chromosome), start=c(clagCultIN$start),
                                 end=c(clagCultIN$end), gene = c(clagCultIN$short), strand=c(clagCultIN$strand)))
clagCult <- sort(clagCult)
##
clagPKNHIN <- read.delim(pknhIN_CyAP)
clagPKNH <- toGRanges(data.frame(chromosome = c(clagPKNHIN$chromosome), start=c(clagPKNHIN$start),
                                 end=c(clagPKNHIN$end), gene = c(clagPKNHIN$short), strand=c(clagPKNHIN$strand)))
clagPKNH <- sort(clagPKNH)
##
clag050IN <- read.delim(sks050IN_CyAP)
clag050 <- toGRanges(data.frame(chromosome = c(clag050IN$chromosome), start=c(clag050IN$start),
                               end=c(clag050IN$end), gene = c(clag050IN$short), strand=c(clag050IN$strand)))
clag050 <- sort(clag050)
##
##
clag058IN <- read.delim(sks058IN_CyAP)
clag058 <- toGRanges(data.frame(chromosome = c(clag058IN$chromosome), start=c(clag058IN$start),
                               end=c(clag058IN$end), gene = c(clag058IN$short), strand=c(clag058IN$strand)))
clag058 <- sort(clag058)
##
##
clag070IN <- read.delim(sks070IN_CyAP)
clag070 <- toGRanges(data.frame(chromosome = c(clag070IN$chromosome), start=c(clag070IN$start),
                               end=c(clag070IN$end), gene = c(clag070IN$short), strand=c(clag070IN$strand)))
clag070 <- sort(clag070)
##
##
clag074IN <- read.delim(sks074IN_CyAP)
clag074 <- toGRanges(data.frame(chromosome = c(clag074IN$chromosome), start=c(clag074IN$start),
                               end=c(clag074IN$end), gene = c(clag074IN$short), strand=c(clag074IN$strand)))
clag074 <- sort(clag074)
##
##
clag078IN <- read.delim(sks078IN_CyAP)
clag078 <- toGRanges(data.frame(chromosome = c(clag078IN$chromosome), start=c(clag078IN$start),
                               end=c(clag078IN$end), gene = c(clag078IN$short), strand=c(clag078IN$strand)))
clag078 <- sort(clag078)
##
##
clag125IN <- read.delim(sks125IN_CyAP)
clag125 <- toGRanges(data.frame(chromosome = c(clag125IN$chromosome), start=c(clag125IN$start),
                               end=c(clag125IN$end), gene = c(clag125IN$short), strand=c(clag125IN$strand)))
clag125 <- sort(clag125)
##
##
clag325IN <- read.delim(sks325IN_CyAP)
clag325 <- toGRanges(data.frame(chromosome = c(clag325IN$chromosome), start=c(clag325IN$start),
                               end=c(clag325IN$end), gene = c(clag325IN$short), strand=c(clag325IN$strand)))
clag325 <- sort(clag325)
##
##
clag331IN <- read.delim(sks331IN_CyAP)
clag331 <- toGRanges(data.frame(chromosome = c(clag331IN$chromosome), start=c(clag331IN$start),
                               end=c(clag331IN$end), gene = c(clag331IN$short), strand=c(clag331IN$strand)))
clag331 <- sort(clag331)
##
##
clag333IN <- read.delim(sks333IN_CyAP)
clag333 <- toGRanges(data.frame(chromosome = c(clag333IN$chromosome), start=c(clag333IN$start),
                               end=c(clag333IN$end), gene = c(clag333IN$short), strand=c(clag333IN$strand)))
clag333 <- sort(clag333)
##
##
clag339IN <- read.delim(sks339IN_CyAP)
clag339 <- toGRanges(data.frame(chromosome = c(clag339IN$chromosome), start=c(clag339IN$start),
                               end=c(clag339IN$end), gene = c(clag339IN$short), strand=c(clag339IN$strand)))
clag339 <- sort(clag339)
##
##
clag344IN <- read.delim(sks344IN_CyAP)
clag344 <- toGRanges(data.frame(chromosome = c(clag344IN$chromosome), start=c(clag344IN$start),
                               end=c(clag344IN$end), gene = c(clag344IN$short), strand=c(clag344IN$strand)))
clag344 <- sort(clag344)
##
####### add in the Duffy binding proteins genes ########
dbp047IN <- read.delim(sks047IN_DBP)
dbp047 <- toGRanges(data.frame(chromosome = c(dbp047IN$chromosome), start=c(dbp047IN$start),
                               end=c(dbp047IN$end), gene = c(dbp047IN$short), strand=c(dbp047IN$strand)))
dbp047 <- sort(dbp047)
##
dbp048IN <- read.delim(sks048IN_DBP)
dbp048 <- toGRanges(data.frame(chromosome = c(dbp048IN$chromosome), start=c(dbp048IN$start),
                               end=c(dbp048IN$end), gene = c(dbp048IN$short), strand=c(dbp048IN$strand)))
dbp048 <- sort(dbp048)
##
dbpCultIN <- read.delim(CultIN_DBP)
dbpCult <- toGRanges(data.frame(chromosome = c(dbpCultIN$chromosome), start=c(dbpCultIN$start),
                                end=c(dbpCultIN$end), gene = c(dbpCultIN$short), strand=c(dbpCultIN$strand)))
dbpCult <- sort(dbpCult)
##
dbpPKNHIN <- read.delim(pknhIN_DBP)
dbpPKNH <- toGRanges(data.frame(chromosome = c(dbpPKNHIN$chromosome), start=c(dbpPKNHIN$start),
                                end=c(dbpPKNHIN$end), gene = c(dbpPKNHIN$short), strand=c(dbpPKNHIN$strand)))
dbpPKNH <- sort(dbpPKNH)
##
dbp050IN <- read.delim(sks050IN_DBP)
dbp050 <- toGRanges(data.frame(chromosome = c(dbp050IN$chromosome), start=c(dbp050IN$start),
                               end=c(dbp050IN$end), gene = c(dbp050IN$short), strand=c(dbp050IN$strand)))
dbp050 <- sort(dbp050)
##
##
dbp058IN <- read.delim(sks058IN_DBP)
dbp058 <- toGRanges(data.frame(chromosome = c(dbp058IN$chromosome), start=c(dbp058IN$start),
                               end=c(dbp058IN$end), gene = c(dbp058IN$short), strand=c(dbp058IN$strand)))
dbp058 <- sort(dbp058)
##
##
dbp070IN <- read.delim(sks070IN_DBP)
dbp070 <- toGRanges(data.frame(chromosome = c(dbp070IN$chromosome), start=c(dbp070IN$start),
                               end=c(dbp070IN$end), gene = c(dbp070IN$short), strand=c(dbp070IN$strand)))
dbp070 <- sort(dbp070)
##
##
dbp074IN <- read.delim(sks074IN_DBP)
dbp074 <- toGRanges(data.frame(chromosome = c(dbp074IN$chromosome), start=c(dbp074IN$start),
                               end=c(dbp074IN$end), gene = c(dbp074IN$short), strand=c(dbp074IN$strand)))
dbp074 <- sort(dbp074)
##
##
dbp078IN <- read.delim(sks078IN_DBP)
dbp078 <- toGRanges(data.frame(chromosome = c(dbp078IN$chromosome), start=c(dbp078IN$start),
                               end=c(dbp078IN$end), gene = c(dbp078IN$short), strand=c(dbp078IN$strand)))
dbp078 <- sort(dbp078)
##
##
dbp125IN <- read.delim(sks125IN_DBP)
dbp125 <- toGRanges(data.frame(chromosome = c(dbp125IN$chromosome), start=c(dbp125IN$start),
                               end=c(dbp125IN$end), gene = c(dbp125IN$short), strand=c(dbp125IN$strand)))
dbp125 <- sort(dbp125)
##
##
dbp325IN <- read.delim(sks325IN_DBP)
dbp325 <- toGRanges(data.frame(chromosome = c(dbp325IN$chromosome), start=c(dbp325IN$start),
                               end=c(dbp325IN$end), gene = c(dbp325IN$short), strand=c(dbp325IN$strand)))
dbp325 <- sort(dbp325)
##
##
dbp331IN <- read.delim(sks331IN_DBP)
dbp331 <- toGRanges(data.frame(chromosome = c(dbp331IN$chromosome), start=c(dbp331IN$start),
                               end=c(dbp331IN$end), gene = c(dbp331IN$short), strand=c(dbp331IN$strand)))
dbp331 <- sort(dbp331)
##
##
dbp333IN <- read.delim(sks333IN_DBP)
dbp333 <- toGRanges(data.frame(chromosome = c(dbp333IN$chromosome), start=c(dbp333IN$start),
                               end=c(dbp333IN$end), gene = c(dbp333IN$short), strand=c(dbp333IN$strand)))
dbp333 <- sort(dbp333)
##
##
dbp339IN <- read.delim(sks339IN_DBP)
dbp339 <- toGRanges(data.frame(chromosome = c(dbp339IN$chromosome), start=c(dbp339IN$start),
                               end=c(dbp339IN$end), gene = c(dbp339IN$short), strand=c(dbp339IN$strand)))
dbp339 <- sort(dbp339)
##
##
dbp344IN <- read.delim(sks344IN_DBP)
dbp344 <- toGRanges(data.frame(chromosome = c(dbp344IN$chromosome), start=c(dbp344IN$start),
                               end=c(dbp344IN$end), gene = c(dbp344IN$short), strand=c(dbp344IN$strand)))
dbp344 <- sort(dbp344)
##
####### add in the ETRAMP genes ########
etramp047IN <- read.delim(sks047IN_ETRAMP)
etramp047 <- toGRanges(data.frame(chromosome = c(etramp047IN$chromosome), start=c(etramp047IN$start),
                                  end=c(etramp047IN$end), gene = c(etramp047IN$short), strand=c(etramp047IN$strand)))
etramp047 <- sort(etramp047)
##
etramp048IN <- read.delim(sks048IN_ETRAMP)
etramp048 <- toGRanges(data.frame(chromosome = c(etramp048IN$chromosome), start=c(etramp048IN$start),
                                  end=c(etramp048IN$end), gene = c(etramp048IN$short), strand=c(etramp048IN$strand)))
etramp048 <- sort(etramp048)
##
etrampCultIN <- read.delim(CultIN_ETRAMP)
etrampCult <- toGRanges(data.frame(chromosome = c(etrampCultIN$chromosome), start=c(etrampCultIN$start),
                                   end=c(etrampCultIN$end), gene = c(etrampCultIN$short), strand=c(etrampCultIN$strand)))
etrampCult <- sort(etrampCult)
##
etrampPKNHIN <- read.delim(pknhIN_ETRAMP)
etrampPKNH <- toGRanges(data.frame(chromosome = c(etrampPKNHIN$chromosome), start=c(etrampPKNHIN$start),
                                   end=c(etrampPKNHIN$end), gene = c(etrampPKNHIN$short), strand=c(etrampPKNHIN$strand)))
etrampPKNH <- sort(etrampPKNH)
##
etramp050IN <- read.delim(sks050IN_ETRAMP)
etramp050 <- toGRanges(data.frame(chromosome = c(etramp050IN$chromosome), start=c(etramp050IN$start),
                               end=c(etramp050IN$end), gene = c(etramp050IN$short), strand=c(etramp050IN$strand)))
etramp050 <- sort(etramp050)
##
##
etramp058IN <- read.delim(sks058IN_ETRAMP)
etramp058 <- toGRanges(data.frame(chromosome = c(etramp058IN$chromosome), start=c(etramp058IN$start),
                               end=c(etramp058IN$end), gene = c(etramp058IN$short), strand=c(etramp058IN$strand)))
etramp058 <- sort(etramp058)
##
##
etramp070IN <- read.delim(sks070IN_ETRAMP)
etramp070 <- toGRanges(data.frame(chromosome = c(etramp070IN$chromosome), start=c(etramp070IN$start),
                               end=c(etramp070IN$end), gene = c(etramp070IN$short), strand=c(etramp070IN$strand)))
etramp070 <- sort(etramp070)
##
##
etramp074IN <- read.delim(sks074IN_ETRAMP)
etramp074 <- toGRanges(data.frame(chromosome = c(etramp074IN$chromosome), start=c(etramp074IN$start),
                               end=c(etramp074IN$end), gene = c(etramp074IN$short), strand=c(etramp074IN$strand)))
etramp074 <- sort(etramp074)
##
##
etramp078IN <- read.delim(sks078IN_ETRAMP)
etramp078 <- toGRanges(data.frame(chromosome = c(etramp078IN$chromosome), start=c(etramp078IN$start),
                               end=c(etramp078IN$end), gene = c(etramp078IN$short), strand=c(etramp078IN$strand)))
etramp078 <- sort(etramp078)
##
##
etramp125IN <- read.delim(sks125IN_ETRAMP)
etramp125 <- toGRanges(data.frame(chromosome = c(etramp125IN$chromosome), start=c(etramp125IN$start),
                               end=c(etramp125IN$end), gene = c(etramp125IN$short), strand=c(etramp125IN$strand)))
etramp125 <- sort(etramp125)
##
##
etramp325IN <- read.delim(sks325IN_ETRAMP)
etramp325 <- toGRanges(data.frame(chromosome = c(etramp325IN$chromosome), start=c(etramp325IN$start),
                               end=c(etramp325IN$end), gene = c(etramp325IN$short), strand=c(etramp325IN$strand)))
etramp325 <- sort(etramp325)
##
##
etramp331IN <- read.delim(sks331IN_ETRAMP)
etramp331 <- toGRanges(data.frame(chromosome = c(etramp331IN$chromosome), start=c(etramp331IN$start),
                               end=c(etramp331IN$end), gene = c(etramp331IN$short), strand=c(etramp331IN$strand)))
etramp331 <- sort(etramp331)
##
##
etramp333IN <- read.delim(sks333IN_ETRAMP)
etramp333 <- toGRanges(data.frame(chromosome = c(etramp333IN$chromosome), start=c(etramp333IN$start),
                               end=c(etramp333IN$end), gene = c(etramp333IN$short), strand=c(etramp333IN$strand)))
etramp333 <- sort(etramp333)
##
##
etramp339IN <- read.delim(sks339IN_ETRAMP)
etramp339 <- toGRanges(data.frame(chromosome = c(etramp339IN$chromosome), start=c(etramp339IN$start),
                               end=c(etramp339IN$end), gene = c(etramp339IN$short), strand=c(etramp339IN$strand)))
etramp339 <- sort(etramp339)
##
##
etramp344IN <- read.delim(sks344IN_ETRAMP)
etramp344 <- toGRanges(data.frame(chromosome = c(etramp344IN$chromosome), start=c(etramp344IN$start),
                               end=c(etramp344IN$end), gene = c(etramp344IN$short), strand=c(etramp344IN$strand)))
etramp344 <- sort(etramp344)
##
####### add in the ABC transporter genes ########
ABC047IN <- read.delim(sks047IN_ABC)
ABC047 <- toGRanges(data.frame(chromosome = c(ABC047IN$chromosome), start=c(ABC047IN$start),
                                 end=c(ABC047IN$end), gene = c(ABC047IN$short), strand=c(ABC047IN$strand)))
ABC047 <- sort(ABC047)
##
ABC048IN <- read.delim(sks048IN_ABC)
ABC048 <- toGRanges(data.frame(chromosome = c(ABC048IN$chromosome), start=c(ABC048IN$start),
                                 end=c(ABC048IN$end), gene = c(ABC048IN$short), strand=c(ABC048IN$strand)))
ABC048 <- sort(ABC048)
##
ABCCultIN <- read.delim(CultIN_ABC)
ABCCult <- toGRanges(data.frame(chromosome = c(ABCCultIN$chromosome), start=c(ABCCultIN$start),
                                  end=c(ABCCultIN$end), gene = c(ABCCultIN$short), strand=c(ABCCultIN$strand)))
ABCCult <- sort(ABCCult)
##
ABCPKNHIN <- read.delim(pknhIN_ABC)
ABCPKNH <- toGRanges(data.frame(chromosome = c(ABCPKNHIN$chromosome), start=c(ABCPKNHIN$start),
                                  end=c(ABCPKNHIN$end), gene = c(ABCPKNHIN$short), strand=c(ABCPKNHIN$strand)))
ABCPKNH <- sort(ABCPKNH)
##
ABC050IN <- read.delim(sks050IN_ABC)
ABC050 <- toGRanges(data.frame(chromosome = c(ABC050IN$chromosome), start=c(ABC050IN$start),
                               end=c(ABC050IN$end), gene = c(ABC050IN$short), strand=c(ABC050IN$strand)))
ABC050 <- sort(ABC050)
##
##
ABC058IN <- read.delim(sks058IN_ABC)
ABC058 <- toGRanges(data.frame(chromosome = c(ABC058IN$chromosome), start=c(ABC058IN$start),
                               end=c(ABC058IN$end), gene = c(ABC058IN$short), strand=c(ABC058IN$strand)))
ABC058 <- sort(ABC058)
##
##
ABC070IN <- read.delim(sks070IN_ABC)
ABC070 <- toGRanges(data.frame(chromosome = c(ABC070IN$chromosome), start=c(ABC070IN$start),
                               end=c(ABC070IN$end), gene = c(ABC070IN$short), strand=c(ABC070IN$strand)))
ABC070 <- sort(ABC070)
##
##
ABC074IN <- read.delim(sks074IN_ABC)
ABC074 <- toGRanges(data.frame(chromosome = c(ABC074IN$chromosome), start=c(ABC074IN$start),
                               end=c(ABC074IN$end), gene = c(ABC074IN$short), strand=c(ABC074IN$strand)))
ABC074 <- sort(ABC074)
##
##
ABC078IN <- read.delim(sks078IN_ABC)
ABC078 <- toGRanges(data.frame(chromosome = c(ABC078IN$chromosome), start=c(ABC078IN$start),
                               end=c(ABC078IN$end), gene = c(ABC078IN$short), strand=c(ABC078IN$strand)))
ABC078 <- sort(ABC078)
##
##
ABC125IN <- read.delim(sks125IN_ABC)
ABC125 <- toGRanges(data.frame(chromosome = c(ABC125IN$chromosome), start=c(ABC125IN$start),
                               end=c(ABC125IN$end), gene = c(ABC125IN$short), strand=c(ABC125IN$strand)))
ABC125 <- sort(ABC125)
##
##
ABC325IN <- read.delim(sks325IN_ABC)
ABC325 <- toGRanges(data.frame(chromosome = c(ABC325IN$chromosome), start=c(ABC325IN$start),
                               end=c(ABC325IN$end), gene = c(ABC325IN$short), strand=c(ABC325IN$strand)))
ABC325 <- sort(ABC325)
##
##
ABC331IN <- read.delim(sks331IN_ABC)
ABC331 <- toGRanges(data.frame(chromosome = c(ABC331IN$chromosome), start=c(ABC331IN$start),
                               end=c(ABC331IN$end), gene = c(ABC331IN$short), strand=c(ABC331IN$strand)))
ABC331 <- sort(ABC331)
##
##
ABC333IN <- read.delim(sks333IN_ABC)
ABC333 <- toGRanges(data.frame(chromosome = c(ABC333IN$chromosome), start=c(ABC333IN$start),
                               end=c(ABC333IN$end), gene = c(ABC333IN$short), strand=c(ABC333IN$strand)))
ABC333 <- sort(ABC333)
##
##
ABC339IN <- read.delim(sks339IN_ABC)
ABC339 <- toGRanges(data.frame(chromosome = c(ABC339IN$chromosome), start=c(ABC339IN$start),
                               end=c(ABC339IN$end), gene = c(ABC339IN$short), strand=c(ABC339IN$strand)))
ABC339 <- sort(ABC339)
##
##
ABC344IN <- read.delim(sks344IN_ABC)
ABC344 <- toGRanges(data.frame(chromosome = c(ABC344IN$chromosome), start=c(ABC344IN$start),
                               end=c(ABC344IN$end), gene = c(ABC344IN$short), strand=c(ABC344IN$strand)))
ABC344 <- sort(ABC344)
##
############################################################################################################
####### sks047 #######
# plot all the genes on the chromosomes and save to file
tiff("sks047.tif", units="px", width=6400, height=3200, res=300)
pp <- getDefaultPlotParams(plot.type = 2)
# set boundaries for the graph
pp$leftmargin <- 0.12
pp$data1outmargin <- 250
pp$data1height <- 300
pp$data2height <- 300
pp$data2outmargin <- 1500
pp$topmargin <- 450
pp$bottommargin <- 400

Graph047 <- plotKaryotype(genome = sks047, ideogram.plotter = NULL, 
                          plot.type = 2, plot.params = pp, cex=1.0)
kpAddCytobandsAsLine(Graph047)
kpAddBaseNumbers(Graph047, tick.dist = 500000, add.units = TRUE,
                 minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
kpAddMainTitle(Graph047, "sks047", cex=1.2)
kpPlotRegions(Graph047, data=Genes047[strand(Genes047)=="+"], avoid.overlapping = FALSE)
kpPlotRegions(Graph047, data=Genes047[strand(Genes047)=="-"], avoid.overlapping = FALSE, data.panel = 2)
kpAddLabels(Graph047, "strand +", cex=0.75, col="#888888", r1=3)
kpAddLabels(Graph047, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3.2)
# MSP genes
kpPlotMarkers(Graph047, data = msp047[strand(msp047$strand)=="+"], labels = msp047$gene, label.color = "purple",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
kpPlotMarkers(Graph047, data = msp047[strand(msp047$strand)=="-"], labels = msp047$gene, label.color = "brown",
              text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown", offset = 0.3)
# #circumsporozoite
kpPlotMarkers(Graph047, data = circ047[strand(circ047$strand)=="+"], labels = circ047$gene, label.color = "dodgerblue2",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
kpPlotMarkers(Graph047, data = circ047[strand(circ047$strand)=="-"], labels = circ047$gene, label.color = "forestgreen",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
# sporozoite invasion genes
kpPlotMarkers(Graph047, data = spor047[strand(spor047$strand)=="+"], labels = spor047$gene, label.color = "indianred",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
kpPlotMarkers(Graph047, data = spor047[strand(spor047$strand)=="-"], labels = spor047$gene, label.color = "darkorange1",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.4, line.color = "darkorange1")
# pknbp
kpPlotMarkers(Graph047, data = pknb047[strand(pknb047$strand)=="+"], labels = pknb047$gene, label.color = "steelblue3",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3", offset = 0.5)
# kpPlotMarkers(Graph047, data = pknb047[strand(pknb047$strand)=="-"], labels = pknb047$gene, label.color = "yellowgreen",
#               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
#               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
#kahrp
kpPlotMarkers(Graph047, data = kahrp047[strand(kahrp047$strand)=="+"], labels = kahrp047$gene, label.color = "red2",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
# kpPlotMarkers(Graph047, data = kahrp047[strand(kahrp047$strand)=="-"], labels = kahrp047$gene, label.color = "royalblue1",
#               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
#               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
#etramp
kpPlotMarkers(Graph047, data = etramp047[strand(etramp047$strand)=="+"], labels = etramp047$gene, label.color = "darkgoldenrod2",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
kpPlotMarkers(Graph047, data = etramp047[strand(etramp047$strand)=="-"], labels = etramp047$gene, label.color = "navyblue",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
#dbp
kpPlotMarkers(Graph047, data = dbp047[strand(dbp047$strand)=="+"], labels = dbp047$gene, label.color = "darkred",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
kpPlotMarkers(Graph047, data = dbp047[strand(dbp047$strand)=="-"], labels = dbp047$gene, label.color = "mediumorchid2",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
#clag
kpPlotMarkers(Graph047, data = clag047[strand(clag047$strand)=="-"], labels = clag047$gene, label.color = "seagreen",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
#abc
kpPlotMarkers(Graph047, data = ABC047[strand(ABC047$strand)=="+"], labels = ABC047$gene, label.color = "red",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
kpPlotMarkers(Graph047, data = ABC047[strand(ABC047$strand)=="-"], labels = ABC047$gene, label.color = "olivedrab4",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
dev.off()
####### sks048 #######
  tiff("sks048.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph048 <- plotKaryotype(genome = sks048, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph048)
  kpAddMainTitle(Graph048, "sks048", cex=1.2)
  kpAddBaseNumbers(Graph048, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(Graph048, data=Genes048[strand(Genes048)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph048, data=Genes048[strand(Genes048)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph048, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph048, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph048, data = msp048[strand(msp048$strand)=="+"], labels = msp048$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph048, data = msp048[strand(msp048$strand)=="-"], labels = msp048$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(Graph048, data = circ048[strand(circ048$strand)=="+"], labels = circ048$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(Graph048, data = circ048[strand(circ048$strand)=="-"], labels = circ048$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(Graph048, data = spor048[strand(spor048$strand)=="+"], labels = spor048$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  kpPlotMarkers(Graph048, data = spor048[strand(spor048$strand)=="-"], labels = spor048$gene, label.color = "darkorange1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  kpPlotMarkers(Graph048, data = pknb048[strand(pknb048$strand)=="+"], labels = pknb048$gene, label.color = "steelblue3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(Graph048, data = pknb048[strand(pknb048$strand)=="-"], labels = pknb048$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  kpPlotMarkers(Graph048, data = kahrp048[strand(kahrp048$strand)=="+"], labels = kahrp048$gene, label.color = "red2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # kpPlotMarkers(Graph048, data = kahrp048[strand(kahrp048$strand)=="-"], labels = kahrp048$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(Graph048, data = etramp048[strand(etramp048$strand)=="+"], labels = etramp048$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(Graph048, data = etramp048[strand(etramp048$strand)=="-"], labels = etramp048$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  kpPlotMarkers(Graph048, data = dbp048[strand(dbp048$strand)=="+"], labels = dbp048$gene, label.color = "darkred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph048, data = dbp048[strand(dbp048$strand)=="-"], labels = dbp048$gene, label.color = "mediumorchid2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(Graph048, data = clag048[strand(clag048$strand)=="-"], labels = clag048$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(Graph048, data = ABC048[strand(ABC048$strand)=="+"], labels = ABC048$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph048, data = ABC048[strand(ABC048$strand)=="-"], labels = ABC048$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()
####### PKNH #######
  tiff("PKNH.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.1
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  GraphPKNH <- plotKaryotype(genome = PKNHy, ideogram.plotter = NULL, 
                             plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(GraphPKNH)
  kpAddMainTitle(GraphPKNH, "PKNH", cex=1.2)
  kpAddBaseNumbers(GraphPKNH, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(GraphPKNH, data=GenesPKNH[strand(GenesPKNH)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(GraphPKNH, data=GenesPKNH[strand(GenesPKNH)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(GraphPKNH, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(GraphPKNH, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(GraphPKNH, data = mspPKNH[strand(mspPKNH$strand)=="+"], labels = mspPKNH$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(GraphPKNH, data = mspPKNH[strand(mspPKNH$strand)=="-"], labels = mspPKNH$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # circumsporozoite
  kpPlotMarkers(GraphPKNH, data = circPKNH[strand(circPKNH$strand)=="+"], labels = circPKNH$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(GraphPKNH, data = circPKNH[strand(circPKNH$strand)=="-"], labels = circPKNH$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(GraphPKNH, data = sporPKNH[strand(sporPKNH$strand)=="+"], labels = sporPKNH$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  kpPlotMarkers(GraphPKNH, data = sporPKNH[strand(sporPKNH$strand)=="-"], labels = sporPKNH$gene, label.color = "darkorange1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  kpPlotMarkers(GraphPKNH, data = pknbPKNH[strand(pknbPKNH$strand)=="+"], labels = pknbPKNH$gene, label.color = "steelblue3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(GraphPKNH, data = pknbPKNH[strand(pknbPKNH$strand)=="-"], labels = pknbPKNH$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  kpPlotMarkers(GraphPKNH, data = kahrpPKNH[strand(kahrpPKNH$strand)=="+"], labels = kahrpPKNH$gene, label.color = "red2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # kpPlotMarkers(GraphPKNH, data = kahrpPKNH[strand(kahrpPKNH$strand)=="-"], labels = kahrpPKNH$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(GraphPKNH, data = etrampPKNH[strand(etrampPKNH$strand)=="+"], labels = etrampPKNH$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(GraphPKNH, data = etrampPKNH[strand(etrampPKNH$strand)=="-"], labels = etrampPKNH$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  kpPlotMarkers(GraphPKNH, data = dbpPKNH[strand(dbpPKNH$strand)=="+"], labels = dbpPKNH$gene, label.color = "darkred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(GraphPKNH, data = dbpPKNH[strand(dbpPKNH$strand)=="-"], labels = dbpPKNH$gene, label.color = "mediumorchid2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(GraphPKNH, data = clagPKNH[strand(clagPKNH$strand)=="-"], labels = clagPKNH$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(GraphPKNH, data = ABCPKNH[strand(ABCPKNH$strand)=="+"], labels = ABCPKNH$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(GraphPKNH, data = ABCPKNH[strand(ABCPKNH$strand)=="-"], labels = ABCPKNH$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()
####### StAPkA1H1 #######
  tiff("Cultured_Pk.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.1
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  GraphCULT <- plotKaryotype(genome = CULT, ideogram.plotter = NULL, 
                             plot.type = 2, plot.params = pp, cex=0.9)
  kpAddCytobandsAsLine(GraphCULT)
  kpAddMainTitle(GraphCULT, "StAPkA1H1", cex=1.2)
  kpAddBaseNumbers(GraphCULT, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(GraphCULT, data=Genes0CULT[strand(Genes0CULT)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(GraphCULT, data=Genes0CULT[strand(Genes0CULT)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(GraphCULT, "strand +", cex=0.8, col="#888888", r1=3)
  kpAddLabels(GraphCULT, "strand -", cex=0.8, data.panel = 2, col="#888888", r1=3.2)
  # MSP genes
  kpPlotMarkers(GraphCULT, data = mspCult[strand(mspCult$strand)=="+"], labels = mspCult$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(GraphCULT, data = mspCult[strand(mspCult$strand)=="-"], labels = mspCult$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(GraphCULT, data = circCult[strand(circCult$strand)=="+"], labels = circCult$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(GraphCULT, data = circCult[strand(circCult$strand)=="-"], labels = circCult$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(GraphCULT, data = sporCult[strand(sporCult$strand)=="+"], labels = sporCult$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  kpPlotMarkers(GraphCULT, data = sporCult[strand(sporCult$strand)=="-"], labels = sporCult$gene, label.color = "darkorange1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  kpPlotMarkers(GraphCULT, data = pknbCult[strand(pknbCult$strand)=="+"], labels = pknbCult$gene, label.color = "steelblue3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(GraphCULT, data = pknbCult[strand(pknbCult$strand)=="-"], labels = pknbCult$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  kpPlotMarkers(GraphCULT, data = kahrpCult[strand(kahrpCult$strand)=="+"], labels = kahrpCult$gene, label.color = "red2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # kpPlotMarkers(GraphCULT, data = kahrpCult[strand(kahrpCult$strand)=="-"], labels = kahrpCult$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(GraphCULT, data = etrampCult[strand(etrampCult$strand)=="+"], labels = etrampCult$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(GraphCULT, data = etrampCult[strand(etrampCult$strand)=="-"], labels = etrampCult$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  kpPlotMarkers(GraphCULT, data = dbpCult[strand(dbpCult$strand)=="+"], labels = dbpCult$gene, label.color = "darkred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(GraphCULT, data = dbpCult[strand(dbpCult$strand)=="-"], labels = dbpCult$gene, label.color = "mediumorchid2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(GraphCULT, data = clagCult[strand(clagCult$strand)=="-"], labels = clagCult$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(GraphCULT, data = ABCCult[strand(ABCCult$strand)=="+"], labels = ABCCult$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(GraphCULT, data = ABCCult[strand(ABCCult$strand)=="-"], labels = ABCCult$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()

####### sks058 #######
  tiff("sks058.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph058 <- plotKaryotype(genome = sks058, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph058)
  kpAddMainTitle(Graph058, "sks058", cex=1.2)
  kpAddBaseNumbers(Graph058, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(Graph058, data=Genes058[strand(Genes058)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph058, data=Genes058[strand(Genes058)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph058, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph058, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph058, data = msp058[strand(msp058$strand)=="+"], labels = msp058$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph058, data = msp058[strand(msp058$strand)=="-"], labels = msp058$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(Graph058, data = circ058[strand(circ058$strand)=="+"], labels = circ058$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(Graph058, data = circ058[strand(circ058$strand)=="-"], labels = circ058$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(Graph058, data = spor058[strand(spor058$strand)=="+"], labels = spor058$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  kpPlotMarkers(Graph058, data = spor058[strand(spor058$strand)=="-"], labels = spor058$gene, label.color = "darkorange1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  kpPlotMarkers(Graph058, data = pknb058[strand(pknb058$strand)=="+"], labels = pknb058$gene, label.color = "steelblue3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(Graph058, data = pknb058[strand(pknb058$strand)=="-"], labels = pknb058$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  kpPlotMarkers(Graph058, data = kahrp058[strand(kahrp058$strand)=="+"], labels = kahrp058$gene, label.color = "red2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # kpPlotMarkers(Graph058, data = kahrp058[strand(kahrp058$strand)=="-"], labels = kahrp058$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(Graph058, data = etramp058[strand(etramp058$strand)=="+"], labels = etramp058$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(Graph058, data = etramp058[strand(etramp058$strand)=="-"], labels = etramp058$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  kpPlotMarkers(Graph058, data = dbp058[strand(dbp058$strand)=="+"], labels = dbp058$gene, label.color = "darkred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph058, data = dbp058[strand(dbp058$strand)=="-"], labels = dbp058$gene, label.color = "mediumorchid2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(Graph058, data = clag058[strand(clag058$strand)=="-"], labels = clag058$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(Graph058, data = ABC058[strand(ABC058$strand)=="+"], labels = ABC058$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph058, data = ABC058[strand(ABC058$strand)=="-"], labels = ABC058$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()
####### sks070 #######
  tiff("sks070.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph070 <- plotKaryotype(genome = sks070, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph070)
  kpAddMainTitle(Graph070, "sks070", cex=1.2)
  kpAddBaseNumbers(Graph070, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(Graph070, data=Genes070[strand(Genes070)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph070, data=Genes070[strand(Genes070)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph070, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph070, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph070, data = msp070[strand(msp070$strand)=="+"], labels = msp070$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph070, data = msp070[strand(msp070$strand)=="-"], labels = msp070$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(Graph070, data = circ070[strand(circ070$strand)=="+"], labels = circ070$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(Graph070, data = circ070[strand(circ070$strand)=="-"], labels = circ070$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(Graph070, data = spor070[strand(spor070$strand)=="+"], labels = spor070$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  kpPlotMarkers(Graph070, data = spor070[strand(spor070$strand)=="-"], labels = spor070$gene, label.color = "darkorange1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  kpPlotMarkers(Graph070, data = pknb070[strand(pknb070$strand)=="+"], labels = pknb070$gene, label.color = "steelblue3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(Graph070, data = pknb070[strand(pknb070$strand)=="-"], labels = pknb070$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  kpPlotMarkers(Graph070, data = kahrp070[strand(kahrp070$strand)=="+"], labels = kahrp070$gene, label.color = "red2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # kpPlotMarkers(Graph070, data = kahrp070[strand(kahrp070$strand)=="-"], labels = kahrp070$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(Graph070, data = etramp070[strand(etramp070$strand)=="+"], labels = etramp070$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(Graph070, data = etramp070[strand(etramp070$strand)=="-"], labels = etramp070$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  kpPlotMarkers(Graph070, data = dbp070[strand(dbp070$strand)=="+"], labels = dbp070$gene, label.color = "darkred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph070, data = dbp070[strand(dbp070$strand)=="-"], labels = dbp070$gene, label.color = "mediumorchid2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(Graph070, data = clag070[strand(clag070$strand)=="-"], labels = clag070$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(Graph070, data = ABC070[strand(ABC070$strand)=="+"], labels = ABC070$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph070, data = ABC070[strand(ABC070$strand)=="-"], labels = ABC070$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()
####### sks074 #######
  tiff("sks074.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph074 <- plotKaryotype(genome = sks074, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph074)
  kpAddMainTitle(Graph074, "sks074", cex=1.2)
  kpAddBaseNumbers(Graph074, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(Graph074, data=Genes074[strand(Genes074)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph074, data=Genes074[strand(Genes074)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph074, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph074, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph074, data = msp074[strand(msp074$strand)=="+"], labels = msp074$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph074, data = msp074[strand(msp074$strand)=="-"], labels = msp074$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(Graph074, data = circ074[strand(circ074$strand)=="+"], labels = circ074$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(Graph074, data = circ074[strand(circ074$strand)=="-"], labels = circ074$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(Graph074, data = spor074[strand(spor074$strand)=="+"], labels = spor074$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  # kpPlotMarkers(Graph074, data = spor074[strand(spor074$strand)=="-"], labels = spor074$gene, label.color = "darkorange1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  # kpPlotMarkers(Graph074, data = pknb074[strand(pknb074$strand)=="+"], labels = pknb074$gene, label.color = "steelblue3",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(Graph074, data = pknb074[strand(pknb074$strand)=="-"], labels = pknb074$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  kpPlotMarkers(Graph074, data = kahrp074[strand(kahrp074$strand)=="+"], labels = kahrp074$gene, label.color = "red2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # kpPlotMarkers(Graph074, data = kahrp074[strand(kahrp074$strand)=="-"], labels = kahrp074$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(Graph074, data = etramp074[strand(etramp074$strand)=="+"], labels = etramp074$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(Graph074, data = etramp074[strand(etramp074$strand)=="-"], labels = etramp074$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  kpPlotMarkers(Graph074, data = dbp074[strand(dbp074$strand)=="+"], labels = dbp074$gene, label.color = "darkred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  # kpPlotMarkers(Graph074, data = dbp074[strand(dbp074$strand)=="-"], labels = dbp074$gene, label.color = "mediumorchid2",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(Graph074, data = clag074[strand(clag074$strand)=="-"], labels = clag074$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(Graph074, data = ABC074[strand(ABC074$strand)=="+"], labels = ABC074$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph074, data = ABC074[strand(ABC074$strand)=="-"], labels = ABC074$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()
  
####### sks078 #######
  tiff("sks078.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph078 <- plotKaryotype(genome = sks078, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph078)
  kpAddMainTitle(Graph078, "sks078", cex=1.2)
  kpAddBaseNumbers(Graph078, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(Graph078, data=Genes078[strand(Genes078)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph078, data=Genes078[strand(Genes078)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph078, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph078, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph078, data = msp078[strand(msp078$strand)=="+"], labels = msp078$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph078, data = msp078[strand(msp078$strand)=="-"], labels = msp078$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(Graph078, data = circ078[strand(circ078$strand)=="+"], labels = circ078$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(Graph078, data = circ078[strand(circ078$strand)=="-"], labels = circ078$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(Graph078, data = spor078[strand(spor078$strand)=="+"], labels = spor078$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  # kpPlotMarkers(Graph078, data = spor078[strand(spor078$strand)=="-"], labels = spor078$gene, label.color = "darkorange1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  # kpPlotMarkers(Graph078, data = pknb078[strand(pknb078$strand)=="+"], labels = pknb078$gene, label.color = "steelblue3",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # # kpPlotMarkers(Graph078, data = pknb078[strand(pknb078$strand)=="-"], labels = pknb078$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  # kpPlotMarkers(Graph078, data = kahrp078[strand(kahrp078$strand)=="+"], labels = kahrp078$gene, label.color = "red2",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # # kpPlotMarkers(Graph078, data = kahrp078[strand(kahrp078$strand)=="-"], labels = kahrp078$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(Graph078, data = etramp078[strand(etramp078$strand)=="+"], labels = etramp078$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(Graph078, data = etramp078[strand(etramp078$strand)=="-"], labels = etramp078$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  # kpPlotMarkers(Graph078, data = dbp078[strand(dbp078$strand)=="+"], labels = dbp078$gene, label.color = "darkred",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  # kpPlotMarkers(Graph078, data = dbp078[strand(dbp078$strand)=="-"], labels = dbp078$gene, label.color = "mediumorchid2",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(Graph078, data = clag078[strand(clag078$strand)=="-"], labels = clag078$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(Graph078, data = ABC078[strand(ABC078$strand)=="+"], labels = ABC078$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph078, data = ABC078[strand(ABC078$strand)=="-"], labels = ABC078$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()
####### sks125 #######
  tiff("sks125.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph125 <- plotKaryotype(genome = sks125, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph125)
  kpAddMainTitle(Graph125, "sks125", cex=1.2)
  kpAddBaseNumbers(Graph125, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(Graph125, data=Genes125[strand(Genes125)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph125, data=Genes125[strand(Genes125)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph125, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph125, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph125, data = msp125[strand(msp125$strand)=="+"], labels = msp125$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph125, data = msp125[strand(msp125$strand)=="-"], labels = msp125$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(Graph125, data = circ125[strand(circ125$strand)=="+"], labels = circ125$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(Graph125, data = circ125[strand(circ125$strand)=="-"], labels = circ125$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(Graph125, data = spor125[strand(spor125$strand)=="+"], labels = spor125$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  # kpPlotMarkers(Graph125, data = spor125[strand(spor125$strand)=="-"], labels = spor125$gene, label.color = "darkorange1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # # pknbp
  # kpPlotMarkers(Graph125, data = pknb125[strand(pknb125$strand)=="+"], labels = pknb125$gene, label.color = "steelblue3",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(Graph125, data = pknb125[strand(pknb125$strand)=="-"], labels = pknb125$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  # kpPlotMarkers(Graph125, data = kahrp125[strand(kahrp125$strand)=="+"], labels = kahrp125$gene, label.color = "red2",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # # kpPlotMarkers(Graph125, data = kahrp125[strand(kahrp125$strand)=="-"], labels = kahrp125$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(Graph125, data = etramp125[strand(etramp125$strand)=="+"], labels = etramp125$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(Graph125, data = etramp125[strand(etramp125$strand)=="-"], labels = etramp125$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  kpPlotMarkers(Graph125, data = dbp125[strand(dbp125$strand)=="+"], labels = dbp125$gene, label.color = "darkred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  # kpPlotMarkers(Graph125, data = dbp125[strand(dbp125$strand)=="-"], labels = dbp125$gene, label.color = "mediumorchid2",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(Graph125, data = clag125[strand(clag125$strand)=="-"], labels = clag125$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(Graph125, data = ABC125[strand(ABC125$strand)=="+"], labels = ABC125$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph125, data = ABC125[strand(ABC125$strand)=="-"], labels = ABC125$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()
####### sks325 #######
  tiff("sks325.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph325 <- plotKaryotype(genome = sks325, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph325)
  kpAddMainTitle(Graph325, "sks325", cex=1.2)
  kpAddBaseNumbers(Graph325, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(Graph325, data=Genes325[strand(Genes325)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph325, data=Genes325[strand(Genes325)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph325, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph325, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph325, data = msp325[strand(msp325$strand)=="+"], labels = msp325$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph325, data = msp325[strand(msp325$strand)=="-"], labels = msp325$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(Graph325, data = circ325[strand(circ325$strand)=="+"], labels = circ325$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(Graph325, data = circ325[strand(circ325$strand)=="-"], labels = circ325$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(Graph325, data = spor325[strand(spor325$strand)=="+"], labels = spor325$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  # kpPlotMarkers(Graph325, data = spor325[strand(spor325$strand)=="-"], labels = spor325$gene, label.color = "darkorange1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  kpPlotMarkers(Graph325, data = pknb325[strand(pknb325$strand)=="+"], labels = pknb325$gene, label.color = "steelblue3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(Graph325, data = pknb325[strand(pknb325$strand)=="-"], labels = pknb325$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  # kpPlotMarkers(Graph325, data = kahrp325[strand(kahrp325$strand)=="+"], labels = kahrp325$gene, label.color = "red2",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # kpPlotMarkers(Graph325, data = kahrp325[strand(kahrp325$strand)=="-"], labels = kahrp325$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(Graph325, data = etramp325[strand(etramp325$strand)=="+"], labels = etramp325$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(Graph325, data = etramp325[strand(etramp325$strand)=="-"], labels = etramp325$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  # kpPlotMarkers(Graph325, data = dbp325[strand(dbp325$strand)=="+"], labels = dbp325$gene, label.color = "darkred",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  # kpPlotMarkers(Graph325, data = dbp325[strand(dbp325$strand)=="-"], labels = dbp325$gene, label.color = "mediumorchid2",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(Graph325, data = clag325[strand(clag325$strand)=="-"], labels = clag325$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(Graph325, data = ABC325[strand(ABC325$strand)=="+"], labels = ABC325$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph325, data = ABC325[strand(ABC325$strand)=="-"], labels = ABC325$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()
  
####### sks331 #######
  tiff("sks331.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph331 <- plotKaryotype(genome = sks331, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph331)
  kpAddMainTitle(Graph331, "sks331", cex=1.2)
  kpAddBaseNumbers(Graph331, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(Graph331, data=Genes331[strand(Genes331)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph331, data=Genes331[strand(Genes331)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph331, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph331, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph331, data = msp331[strand(msp331$strand)=="+"], labels = msp331$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph331, data = msp331[strand(msp331$strand)=="-"], labels = msp331$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(Graph331, data = circ331[strand(circ331$strand)=="+"], labels = circ331$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(Graph331, data = circ331[strand(circ331$strand)=="-"], labels = circ331$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(Graph331, data = spor331[strand(spor331$strand)=="+"], labels = spor331$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  # kpPlotMarkers(Graph331, data = spor331[strand(spor331$strand)=="-"], labels = spor331$gene, label.color = "darkorange1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  # kpPlotMarkers(Graph331, data = pknb331[strand(pknb331$strand)=="+"], labels = pknb331$gene, label.color = "steelblue3",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(Graph331, data = pknb331[strand(pknb331$strand)=="-"], labels = pknb331$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  # kpPlotMarkers(Graph331, data = kahrp331[strand(kahrp331$strand)=="+"], labels = kahrp331$gene, label.color = "red2",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # kpPlotMarkers(Graph331, data = kahrp331[strand(kahrp331$strand)=="-"], labels = kahrp331$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(Graph331, data = etramp331[strand(etramp331$strand)=="+"], labels = etramp331$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(Graph331, data = etramp331[strand(etramp331$strand)=="-"], labels = etramp331$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  # kpPlotMarkers(Graph331, data = dbp331[strand(dbp331$strand)=="+"], labels = dbp331$gene, label.color = "darkred",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  # kpPlotMarkers(Graph331, data = dbp331[strand(dbp331$strand)=="-"], labels = dbp331$gene, label.color = "mediumorchid2",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(Graph331, data = clag331[strand(clag331$strand)=="-"], labels = clag331$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(Graph331, data = ABC331[strand(ABC331$strand)=="+"], labels = ABC331$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph331, data = ABC331[strand(ABC331$strand)=="-"], labels = ABC331$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()

####### sks333 #######
  tiff("sks333.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph333 <- plotKaryotype(genome = sks333, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph333)
  kpAddMainTitle(Graph333, "sks333", cex=1.2)
  kpAddBaseNumbers(Graph333, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(Graph333, data=Genes333[strand(Genes333)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph333, data=Genes333[strand(Genes333)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph333, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph333, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph333, data = msp333[strand(msp333$strand)=="+"], labels = msp333$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph333, data = msp333[strand(msp333$strand)=="-"], labels = msp333$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(Graph333, data = circ333[strand(circ333$strand)=="+"], labels = circ333$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(Graph333, data = circ333[strand(circ333$strand)=="-"], labels = circ333$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(Graph333, data = spor333[strand(spor333$strand)=="+"], labels = spor333$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  kpPlotMarkers(Graph333, data = spor333[strand(spor333$strand)=="-"], labels = spor333$gene, label.color = "darkorange1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  kpPlotMarkers(Graph333, data = pknb333[strand(pknb333$strand)=="+"], labels = pknb333$gene, label.color = "steelblue3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(Graph333, data = pknb333[strand(pknb333$strand)=="-"], labels = pknb333$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  kpPlotMarkers(Graph333, data = kahrp333[strand(kahrp333$strand)=="+"], labels = kahrp333$gene, label.color = "red2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # kpPlotMarkers(Graph333, data = kahrp333[strand(kahrp333$strand)=="-"], labels = kahrp333$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(Graph333, data = etramp333[strand(etramp333$strand)=="+"], labels = etramp333$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(Graph333, data = etramp333[strand(etramp333$strand)=="-"], labels = etramp333$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  kpPlotMarkers(Graph333, data = dbp333[strand(dbp333$strand)=="+"], labels = dbp333$gene, label.color = "darkred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  # kpPlotMarkers(Graph333, data = dbp333[strand(dbp333$strand)=="-"], labels = dbp333$gene, label.color = "mediumorchid2",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(Graph333, data = clag333[strand(clag333$strand)=="-"], labels = clag333$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(Graph333, data = ABC333[strand(ABC333$strand)=="+"], labels = ABC333$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph333, data = ABC333[strand(ABC333$strand)=="-"], labels = ABC333$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()
  
####### sks339 #######
  tiff("sks339.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph339 <- plotKaryotype(genome = sks339, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph339)
  kpAddMainTitle(Graph339, "sks339", cex=1.2)
  kpAddBaseNumbers(Graph339, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(Graph339, data=Genes339[strand(Genes339)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph339, data=Genes339[strand(Genes339)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph339, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph339, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph339, data = msp339[strand(msp339$strand)=="+"], labels = msp339$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph339, data = msp339[strand(msp339$strand)=="-"], labels = msp339$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(Graph339, data = circ339[strand(circ339$strand)=="+"], labels = circ339$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(Graph339, data = circ339[strand(circ339$strand)=="-"], labels = circ339$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(Graph339, data = spor339[strand(spor339$strand)=="+"], labels = spor339$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  kpPlotMarkers(Graph339, data = spor339[strand(spor339$strand)=="-"], labels = spor339$gene, label.color = "darkorange1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  kpPlotMarkers(Graph339, data = pknb339[strand(pknb339$strand)=="+"], labels = pknb339$gene, label.color = "steelblue3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(Graph339, data = pknb339[strand(pknb339$strand)=="-"], labels = pknb339$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  kpPlotMarkers(Graph339, data = kahrp339[strand(kahrp339$strand)=="+"], labels = kahrp339$gene, label.color = "red2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # kpPlotMarkers(Graph339, data = kahrp339[strand(kahrp339$strand)=="-"], labels = kahrp339$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(Graph339, data = etramp339[strand(etramp339$strand)=="+"], labels = etramp339$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(Graph339, data = etramp339[strand(etramp339$strand)=="-"], labels = etramp339$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  kpPlotMarkers(Graph339, data = dbp339[strand(dbp339$strand)=="+"], labels = dbp339$gene, label.color = "darkred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph339, data = dbp339[strand(dbp339$strand)=="-"], labels = dbp339$gene, label.color = "mediumorchid2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(Graph339, data = clag339[strand(clag339$strand)=="-"], labels = clag339$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(Graph339, data = ABC339[strand(ABC339$strand)=="+"], labels = ABC339$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph339, data = ABC339[strand(ABC339$strand)=="-"], labels = ABC339$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()
  
####### sks344 #######
  tiff("sks344.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph344 <- plotKaryotype(genome = sks344, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph344)
  kpAddMainTitle(Graph344, "sks344", cex=1.2)
  kpAddBaseNumbers(Graph344, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
  kpPlotRegions(Graph344, data=Genes344[strand(Genes344)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph344, data=Genes344[strand(Genes344)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph344, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph344, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph344, data = msp344[strand(msp344$strand)=="+"], labels = msp344$gene, label.color = "purple",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph344, data = msp344[strand(msp344$strand)=="-"], labels = msp344$gene, label.color = "brown",
                text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  # #circumsporozoite
  kpPlotMarkers(Graph344, data = circ344[strand(circ344$strand)=="+"], labels = circ344$gene, label.color = "dodgerblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "dodgerblue2")
  kpPlotMarkers(Graph344, data = circ344[strand(circ344$strand)=="-"], labels = circ344$gene, label.color = "forestgreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.1, line.color = "forestgreen")
  # sporozoite invasion genes
  kpPlotMarkers(Graph344, data = spor344[strand(spor344$strand)=="+"], labels = spor344$gene, label.color = "indianred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "indianred")
  kpPlotMarkers(Graph344, data = spor344[strand(spor344$strand)=="-"], labels = spor344$gene, label.color = "darkorange1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # pknbp
  kpPlotMarkers(Graph344, data = pknb344[strand(pknb344$strand)=="+"], labels = pknb344$gene, label.color = "steelblue3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
  # kpPlotMarkers(Graph344, data = pknb344[strand(pknb344$strand)=="-"], labels = pknb344$gene, label.color = "yellowgreen",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #kahrp
  kpPlotMarkers(Graph344, data = kahrp344[strand(kahrp344$strand)=="+"], labels = kahrp344$gene, label.color = "red2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "red2")
  # kpPlotMarkers(Graph344, data = kahrp344[strand(kahrp344$strand)=="-"], labels = kahrp344$gene, label.color = "royalblue1",
  #               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
  #etramp
  kpPlotMarkers(Graph344, data = etramp344[strand(etramp344$strand)=="+"], labels = etramp344$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2")
  kpPlotMarkers(Graph344, data = etramp344[strand(etramp344$strand)=="-"], labels = etramp344$gene, label.color = "navyblue",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "navyblue")
  #dbp
  kpPlotMarkers(Graph344, data = dbp344[strand(dbp344$strand)=="+"], labels = dbp344$gene, label.color = "darkred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph344, data = dbp344[strand(dbp344$strand)=="-"], labels = dbp344$gene, label.color = "mediumorchid2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.6, line.color = "mediumorchid2", ignore.chromosome.ends = TRUE)
  #clag
  kpPlotMarkers(Graph344, data = clag344[strand(clag344$strand)=="-"], labels = clag344$gene, label.color = "seagreen",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  #abc
  kpPlotMarkers(Graph344, data = ABC344[strand(ABC344$strand)=="+"], labels = ABC344$gene, label.color = "red",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "red", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph344, data = ABC344[strand(ABC344$strand)=="-"], labels = ABC344$gene, label.color = "olivedrab4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "olivedrab4")
  dev.off()
  
###### Combined Output ######
# png("Multiple gene families in isolates2.png", width = 2400, height = 1200)
# tiff("Cultured_Pk.tif", units="px", width=6400, height=3200, res=300)
# p1 <- as.ggplot(expression(plot047()))
# p2 <- as.ggplot(expression(plot048()))
# p3 <- as.ggplot(expression(plotCULT()))
# p4 <- as.ggplot(expression(plotPKNH()))
# plot_grid(p4, p3, p1, p2, ncol=2)
# dev.off()
#### make a table #####
all.table <- data.frame(
  Genes = c("Circumsporozoite protein", "Cytoadherence linked asexual protein/gene", "Duffy binding/Duffy-antigen protein", "Early transcribed membrane protein",
            "Knob-associated histidine-rich protein","Merozoite surface protein", "Reticulocyte binding protein", 
            "Sporozoite invasion-associated protein", "ABC Transporter"),
  abbr. = c("CSP/CS-TRAP", "CLAG", "DBP", "ETRAMP", "KAHRP", "MSP", "Pknbp", "SPIAP", "ABCtrp"),
  PKNH = c(nrow(circPKNHIN),nrow(clagPKNHIN), nrow(dbpPKNHIN), nrow(etrampPKNHIN), nrow(kahrpPKNHIN),
           nrow(mspPKnhIN), nrow(pknbPKNHIN), nrow(sporPKNHIN), nrow(ABCPKNHIN)),
  PkA1H1 = c(nrow(circCultIN), nrow(clagCultIN), nrow(dbpCultIN), nrow(etrampCultIN), nrow(kahrpCultIN),
             nrow(mspCultIN), nrow(pknbCultIN), nrow(sporCultIN), nrow(ABCCultIN)),
  sks047 = c(nrow(circ047IN), nrow(clag047IN), nrow(dbp047IN), nrow(etramp047IN), nrow(kahrp047IN),
             nrow(msp047IN), nrow(pknb047IN), nrow(spor047IN), nrow(ABC047IN)),
  sks048 = c(nrow(circ048IN), nrow(clag048IN), nrow(dbp048IN), nrow(etramp048IN), nrow(kahrp048IN),
             nrow(msp048IN), nrow(pknb048IN), nrow(spor048IN), nrow(ABC048IN)),
  sks050 = c(nrow(circ050IN), nrow(clag050IN), nrow(dbp050IN), nrow(etramp050IN), nrow(kahrp050IN),
             nrow(msp050IN), nrow(pknb050IN), nrow(spor050IN), nrow(ABC050IN)),
  sks058 = c(nrow(circ058IN), nrow(clag058IN), nrow(dbp058IN), nrow(etramp058IN), nrow(kahrp058IN),
             nrow(msp058IN), nrow(pknb058IN), nrow(spor058IN), nrow(ABC058IN)),
  sks070 = c(nrow(circ070IN), nrow(clag070IN), nrow(dbp070IN), nrow(etramp070IN), nrow(kahrp070IN),
             nrow(msp070IN), nrow(pknb070IN), nrow(spor070IN), nrow(ABC070IN)),
  sks074 = c(nrow(circ074IN), nrow(clag074IN), nrow(dbp074IN), nrow(etramp074IN), nrow(kahrp074IN),
             nrow(msp074IN), nrow(pknb074IN), nrow(spor074IN), nrow(ABC074IN)),
  sks078 = c(nrow(circ078IN), nrow(clag078IN), nrow(dbp078IN), nrow(etramp078IN), nrow(kahrp078IN),
             nrow(msp078IN), nrow(pknb078IN), nrow(spor078IN), nrow(ABC078IN)),
  sks125 = c(nrow(circ125IN), nrow(clag125IN), nrow(dbp125IN), nrow(etramp125IN), nrow(kahrp125IN),
             nrow(msp125IN), nrow(pknb125IN), nrow(spor125IN), nrow(ABC125IN)),
  sks325 = c(nrow(circ325IN), nrow(clag325IN), nrow(dbp325IN), nrow(etramp325IN), nrow(kahrp325IN),
             nrow(msp325IN), nrow(pknb325IN), nrow(spor325IN), nrow(ABC325IN)),
  sks331 = c(nrow(circ331IN), nrow(clag331IN), nrow(dbp331IN), nrow(etramp331IN), nrow(kahrp331IN),
             nrow(msp331IN), nrow(pknb331IN), nrow(spor331IN), nrow(ABC331IN)),
  sks333 = c(nrow(circ333IN), nrow(clag333IN), nrow(dbp333IN), nrow(etramp333IN), nrow(kahrp333IN),
             nrow(msp333IN), nrow(pknb333IN), nrow(spor333IN), nrow(ABC333IN)),
  sks339 = c(nrow(circ339IN), nrow(clag339IN), nrow(dbp339IN), nrow(etramp339IN), nrow(kahrp339IN),
             nrow(msp339IN), nrow(pknb339IN), nrow(spor339IN), nrow(ABC339IN)),
  sks344 = c(nrow(circ344IN), nrow(clag344IN), nrow(dbp344IN), nrow(etramp344IN), nrow(kahrp344IN),
             nrow(msp344IN), nrow(pknb344IN), nrow(spor344IN), nrow(ABC344IN))
)
# formattable(all.table, align = c("l", "c", "r", "r", "r", "r","r"), list(
#   PKNH = color_tile("white", "forestgreen"),
#   PkA1H1 = color_tile("white", "forestgreen"),
#   sks047 = color_tile("white", "forestgreen"),
#   sks048 = color_tile("white", "forestgreen")
# ))

formattable(all.table, align = c("l", "c", "r", "r", "r", "r", "r", "r", "r", "r", "r", "r", "r", "r", "r", "r", "r"))

# graphhy <- plotKaryotype(genome=CULT,ideogram.plotter = NULL, 
#                          plot.type = 2, cex=0.9, chromosomes = "PKA1H1_STAND_01")
# kpAddBaseNumbers(graphhy, tick.dist = 50000, add.units = TRUE,minor.tick.dist = 10000)
# CULT

#####################################
# example of singular plotting for one isolate
###########################################
# pp <- getDefaultPlotParams(plot.type = 2)
# # set boundaries for the graph
# pp$leftmargin <- 0.1
# pp$data1outmargin <- 350
# pp$data2outmargin <- 450
# pp$topmargin <- 450
# pp$bottommargin <- 400
# png("sks047.png", width = 1800, height = 900)
# Graph047 <- plotKaryotype(genome = sks047, ideogram.plotter = NULL,
#                           plot.type = 2, plot.params = pp, cex=0.9)
# #kpAddBaseNumbers(Graph047, add.units = TRUE)
# kpAddCytobandsAsLine(Graph047)
# kpAddMainTitle(Graph047, "sks047 chromosome layout", cex=2)
# kpPlotRegions(Graph047, data=Genes047[strand(Genes047)=="+"], avoid.overlapping = FALSE, col="royal blue")
# kpPlotRegions(Graph047, data=Genes047[strand(Genes047)=="-"], avoid.overlapping = FALSE, data.panel = 2, col="maroon")
# kpAddLabels(Graph047, "strand +", cex=0.8, col="#888888", r1=2.5)
# kpAddLabels(Graph047, "strand -", cex=0.8, data.panel = 2, col="#888888", r1=2.5)
# kpPlotMarkers(myGraph, data = vva[strand(vva$strand)=="+"], labels = vva$gene, label.color = "purple",
#               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
#               data.panel = 1, label.dist = 0.001, max.iter = 1000)
# kpPlotMarkers(myGraph, data = vva[strand(vva$strand)=="-"], labels = vva$gene, label.color = "brown",
#               text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
#               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9)
# dev.off()
