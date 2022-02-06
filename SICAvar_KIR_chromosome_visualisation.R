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
sks047IN_KIR <- "sks047_KIR.gff3"
sks048IN_KIR <- "sks048_KIR.gff3"
CultIN_KIR <- "StAPkA1H1_KIR.gff3"
pknhIN_KIR <- "PKNH_KIR.gff3"
sks050IN_KIR <- "sks050_KIR.gff3"
sks058IN_KIR <- "sks058_KIR.gff3"
sks070IN_KIR <- "sks070_KIR.gff3"
sks074IN_KIR <- "sks074_KIR.gff3"
sks078IN_KIR <- "sks078_KIR.gff3"
sks125IN_KIR <- "sks125_KIR.gff3"
sks325IN_KIR <- "sks325_KIR.gff3"
sks331IN_KIR <- "sks331_KIR.gff3"
sks333IN_KIR <- "sks333_KIR.gff3"
sks339IN_KIR <- "sks339_KIR.gff3"
sks344IN_KIR <- "sks344_KIR.gff3"
#
sks047IN_SICAvars <- "sks047_SICAvars.gff3"
sks048IN_SICAvars <- "sks048_SICAvars.gff3"
CultIN_SICAvars <- "StAPkA1H1_SICAvars.gff3"
pknhIN_SICAvars <- "Pknowlesi_SICAvars.gff3"
sks050IN_SICAvars <- "sks050_SICAvars.gff3"
sks058IN_SICAvars <- "sks058_SICAvars.gff3"
sks070IN_SICAvars <- "sks070_SICAvars.gff3"
sks074IN_SICAvars <- "sks074_SICAvars.gff3"
sks078IN_SICAvars <- "sks078_SICAvars.gff3"
sks125IN_SICAvars <- "sks125_SICAvars.gff3"
sks325IN_SICAvars <- "sks325_SICAvars.gff3"
sks331IN_SICAvars <- "sks331_SICAvars.gff3"
sks333IN_SICAvars <- "sks333_SICAvars.gff3"
sks339IN_SICAvars <- "sks339_SICAvars.gff3"
sks344IN_SICAvars <- "sks344_SICAvars.gff3"
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
########################################################################################################
####### add in the KIRs genes ########
kir047IN <- read.delim(sks047IN_KIR)
kir047 <- toGRanges(data.frame(chromosome = c(kir047IN$chromosome), start=c(kir047IN$start),
                               end=c(kir047IN$end), gene = c(kir047IN$short), strand=c(kir047IN$strand)))
kir047 <- sort(kir047)
##
kir048IN <- read.delim(sks048IN_KIR)
kir048 <- toGRanges(data.frame(chromosome = c(kir048IN$chromosome), start=c(kir048IN$start),
                               end=c(kir048IN$end), gene = c(kir048IN$short), strand=c(kir048IN$strand)))
kir048 <- sort(kir048)
##
kirCultIN <- read.delim(CultIN_KIR)
kirCult <- toGRanges(data.frame(chromosome = c(kirCultIN$chromosome), start=c(kirCultIN$start),
                                end=c(kirCultIN$end), gene = c(kirCultIN$short), strand=c(kirCultIN$strand)))
kirCult <- sort(kirCult)
##
kirPKnhIN <- read.delim(pknhIN_KIR)
kirPKNH <- toGRanges(data.frame(chromosome = c(kirPKnhIN$chromosome), start=c(kirPKnhIN$start),
                                end=c(kirPKnhIN$end), gene = c(kirPKnhIN$short), strand=c(kirPKnhIN$strand)))
kirPKNH <- sort(kirPKNH)
##
kir050IN <- read.delim(sks050IN_KIR)
kir050 <- toGRanges(data.frame(chromosome = c(kir050IN$chromosome), start=c(kir050IN$start),
                               end=c(kir050IN$end), gene = c(kir050IN$short), strand=c(kir050IN$strand)))
kir050 <- sort(kir050)
##
##
kir058IN <- read.delim(sks058IN_KIR)
kir058 <- toGRanges(data.frame(chromosome = c(kir058IN$chromosome), start=c(kir058IN$start),
                               end=c(kir058IN$end), gene = c(kir058IN$short), strand=c(kir058IN$strand)))
kir058 <- sort(kir058)
##
##
kir070IN <- read.delim(sks070IN_KIR)
kir070 <- toGRanges(data.frame(chromosome = c(kir070IN$chromosome), start=c(kir070IN$start),
                               end=c(kir070IN$end), gene = c(kir070IN$short), strand=c(kir070IN$strand)))
kir070 <- sort(kir070)
##
##
kir074IN <- read.delim(sks074IN_KIR)
kir074 <- toGRanges(data.frame(chromosome = c(kir074IN$chromosome), start=c(kir074IN$start),
                               end=c(kir074IN$end), gene = c(kir074IN$short), strand=c(kir074IN$strand)))
kir074 <- sort(kir074)
##
##
kir078IN <- read.delim(sks078IN_KIR)
kir078 <- toGRanges(data.frame(chromosome = c(kir078IN$chromosome), start=c(kir078IN$start),
                               end=c(kir078IN$end), gene = c(kir078IN$short), strand=c(kir078IN$strand)))
kir078 <- sort(kir078)
##
##
kir125IN <- read.delim(sks125IN_KIR)
kir125 <- toGRanges(data.frame(chromosome = c(kir125IN$chromosome), start=c(kir125IN$start),
                               end=c(kir125IN$end), gene = c(kir125IN$short), strand=c(kir125IN$strand)))
kir125 <- sort(kir125)
##
##
kir325IN <- read.delim(sks325IN_KIR)
kir325 <- toGRanges(data.frame(chromosome = c(kir325IN$chromosome), start=c(kir325IN$start),
                               end=c(kir325IN$end), gene = c(kir325IN$short), strand=c(kir325IN$strand)))
kir325 <- sort(kir325)
##
##
kir331IN <- read.delim(sks331IN_KIR)
kir331 <- toGRanges(data.frame(chromosome = c(kir331IN$chromosome), start=c(kir331IN$start),
                               end=c(kir331IN$end), gene = c(kir331IN$short), strand=c(kir331IN$strand)))
kir331 <- sort(kir331)
##
##
kir333IN <- read.delim(sks333IN_KIR)
kir333 <- toGRanges(data.frame(chromosome = c(kir333IN$chromosome), start=c(kir333IN$start),
                               end=c(kir333IN$end), gene = c(kir333IN$short), strand=c(kir333IN$strand)))
kir333 <- sort(kir333)
##
##
kir339IN <- read.delim(sks339IN_KIR)
kir339 <- toGRanges(data.frame(chromosome = c(kir339IN$chromosome), start=c(kir339IN$start),
                               end=c(kir339IN$end), gene = c(kir339IN$short), strand=c(kir339IN$strand)))
kir339 <- sort(kir339)
##
##
kir344IN <- read.delim(sks344IN_KIR)
kir344 <- toGRanges(data.frame(chromosome = c(kir344IN$chromosome), start=c(kir344IN$start),
                               end=c(kir344IN$end), gene = c(kir344IN$short), strand=c(kir344IN$strand)))
kir344 <- sort(kir344)
####### add in the SICAvar genes ########
SICA047IN <- read.delim(sks047IN_SICAvars)
SICA047 <- toGRanges(data.frame(chromosome = c(SICA047IN$chromosome), start=c(SICA047IN$start),
                               end=c(SICA047IN$end), gene = c(SICA047IN$short), strand=c(SICA047IN$strand)))
SICA047 <- sort(SICA047)
##
SICA048IN <- read.delim(sks048IN_SICAvars)
SICA048 <- toGRanges(data.frame(chromosome = c(SICA048IN$chromosome), start=c(SICA048IN$start),
                               end=c(SICA048IN$end), gene = c(SICA048IN$short), strand=c(SICA048IN$strand)))
SICA048 <- sort(SICA048)
##
SICACultIN <- read.delim(CultIN_SICAvars)
SICACult <- toGRanges(data.frame(chromosome = c(SICACultIN$chromosome), start=c(SICACultIN$start),
                                end=c(SICACultIN$end), gene = c(SICACultIN$short), strand=c(SICACultIN$strand)))
SICACult <- sort(SICACult)
##
SICAPKnhIN <- read.delim(pknhIN_SICAvars)
SICAPKNH <- toGRanges(data.frame(chromosome = c(SICAPKnhIN$chromosome), start=c(SICAPKnhIN$start),
                                end=c(SICAPKnhIN$end), gene = c(SICAPKnhIN$short), strand=c(SICAPKnhIN$strand)))
SICAPKNH <- sort(SICAPKNH)
##
SICA050IN <- read.delim(sks050IN_SICAvars)
SICA050 <- toGRanges(data.frame(chromosome = c(SICA050IN$chromosome), start=c(SICA050IN$start),
                               end=c(SICA050IN$end), gene = c(SICA050IN$short), strand=c(SICA050IN$strand)))
SICA050 <- sort(SICA050)
##
##
SICA058IN <- read.delim(sks058IN_SICAvars)
SICA058 <- toGRanges(data.frame(chromosome = c(SICA058IN$chromosome), start=c(SICA058IN$start),
                               end=c(SICA058IN$end), gene = c(SICA058IN$short), strand=c(SICA058IN$strand)))
SICA058 <- sort(SICA058)
##
##
SICA070IN <- read.delim(sks070IN_SICAvars)
SICA070 <- toGRanges(data.frame(chromosome = c(SICA070IN$chromosome), start=c(SICA070IN$start),
                               end=c(SICA070IN$end), gene = c(SICA070IN$short), strand=c(SICA070IN$strand)))
SICA070 <- sort(SICA070)
##
##
SICA074IN <- read.delim(sks074IN_SICAvars)
SICA074 <- toGRanges(data.frame(chromosome = c(SICA074IN$chromosome), start=c(SICA074IN$start),
                               end=c(SICA074IN$end), gene = c(SICA074IN$short), strand=c(SICA074IN$strand)))
SICA074 <- sort(SICA074)
##
##
SICA078IN <- read.delim(sks078IN_SICAvars)
SICA078 <- toGRanges(data.frame(chromosome = c(SICA078IN$chromosome), start=c(SICA078IN$start),
                               end=c(SICA078IN$end), gene = c(SICA078IN$short), strand=c(SICA078IN$strand)))
SICA078 <- sort(SICA078)
##
##
SICA125IN <- read.delim(sks125IN_SICAvars)
SICA125 <- toGRanges(data.frame(chromosome = c(SICA125IN$chromosome), start=c(SICA125IN$start),
                               end=c(SICA125IN$end), gene = c(SICA125IN$short), strand=c(SICA125IN$strand)))
SICA125 <- sort(SICA125)
##
##
SICA325IN <- read.delim(sks325IN_SICAvars)
SICA325 <- toGRanges(data.frame(chromosome = c(SICA325IN$chromosome), start=c(SICA325IN$start),
                               end=c(SICA325IN$end), gene = c(SICA325IN$short), strand=c(SICA325IN$strand)))
SICA325 <- sort(SICA325)
##
##
SICA331IN <- read.delim(sks331IN_SICAvars)
SICA331 <- toGRanges(data.frame(chromosome = c(SICA331IN$chromosome), start=c(SICA331IN$start),
                               end=c(SICA331IN$end), gene = c(SICA331IN$short), strand=c(SICA331IN$strand)))
SICA331 <- sort(SICA331)
##
##
SICA333IN <- read.delim(sks333IN_SICAvars)
SICA333 <- toGRanges(data.frame(chromosome = c(SICA333IN$chromosome), start=c(SICA333IN$start),
                               end=c(SICA333IN$end), gene = c(SICA333IN$short), strand=c(SICA333IN$strand)))
SICA333 <- sort(SICA333)
##
##
SICA339IN <- read.delim(sks339IN_SICAvars)
SICA339 <- toGRanges(data.frame(chromosome = c(SICA339IN$chromosome), start=c(SICA339IN$start),
                               end=c(SICA339IN$end), gene = c(SICA339IN$short), strand=c(SICA339IN$strand)))
SICA339 <- sort(SICA339)
##
##
SICA344IN <- read.delim(sks344IN_SICAvars)
SICA344 <- toGRanges(data.frame(chromosome = c(SICA344IN$chromosome), start=c(SICA344IN$start),
                               end=c(SICA344IN$end), gene = c(SICA344IN$short), strand=c(SICA344IN$strand)))
SICA344 <- sort(SICA344)

# show_col(viridis_pal()(20))
# q_colors =  4 # for no particular reason
# v_colors =  viridis(q_colors, option = "magma")
# v_colors

############################################################################################################
# plot all the genes on the chromosomes and save to file
#plot047 <- function() {
  #png("sks047SICAvars.png", width = 1600, height = 800)
  tiff("sks047SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
                   minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpAddMainTitle(Graph047, "sks047", cex=1.2)
  kpPlotRegions(Graph047, data=Genes047[strand(Genes047)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph047, data=Genes047[strand(Genes047)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph047, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph047, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3.2)
  #kirs #magenta, turqoise1
  kpPlotMarkers(Graph047, data = kir047[strand(kir047$strand)=="+"], labels = kir047$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph047, data = kir047[strand(kir047$strand)=="-"], labels = kir047$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 5)
  #sicavars darkred seagreen
  kpPlotMarkers(Graph047, data = SICA047[strand(SICA047$strand)=="+"], labels = SICA047$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph047, data = SICA047[strand(SICA047$strand)=="-"], labels = SICA047$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  dev.off()
#}

#plot048 <- function() {
  #png("sks048SICAvars.png", width = 1600, height = 800)
  tiff("sks048SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(Graph048, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph048, data=Genes048[strand(Genes048)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph048, data=Genes048[strand(Genes048)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph048, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph048, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
    #kirs
  kpPlotMarkers(Graph048, data = kir048[strand(kir048$strand)=="+"], labels = kir048$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph048, data = kir048[strand(kir048$strand)=="-"], labels = kir048$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph048, data = SICA048[strand(SICA048$strand)=="+"], labels = SICA048$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph048, data = SICA048[strand(SICA048$strand)=="-"], labels = SICA048$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
#}

#plotPKNH <- function() {
  # png("PKNHSICAvars.png", width = 1600, height = 800)
  tiff("PKNHSICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(GraphPKNH, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(GraphPKNH, data=GenesPKNH[strand(GenesPKNH)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(GraphPKNH, data=GenesPKNH[strand(GenesPKNH)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(GraphPKNH, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(GraphPKNH, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(GraphPKNH, data = kirPKNH[strand(kirPKNH$strand)=="+"], labels = kirPKNH$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(GraphPKNH, data = kirPKNH[strand(kirPKNH$strand)=="-"], labels = kirPKNH$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(GraphPKNH, data = SICAPKNH[strand(SICAPKNH$strand)=="+"], labels = SICAPKNH$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(GraphPKNH, data = SICAPKNH[strand(SICAPKNH$strand)=="-"], labels = SICAPKNH$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  dev.off()
#}

#plotCULT <- function() {
  # png(".png", width = 1600, height = 800)
  tiff("StAPkA1H1SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.1
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  GraphCULT <- plotKaryotype(genome = CULT, ideogram.plotter = NULL, 
                             plot.type = 2, plot.params = pp, cex=1)
  kpAddCytobandsAsLine(GraphCULT)
  kpAddMainTitle(GraphCULT, "StAPkA1H1", cex=1.2)
  kpAddBaseNumbers(GraphCULT, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(GraphCULT, data=Genes0CULT[strand(Genes0CULT)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(GraphCULT, data=Genes0CULT[strand(Genes0CULT)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(GraphCULT, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(GraphCULT, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3.2)
  #kirs
  kpPlotMarkers(GraphCULT, data = kirCult[strand(kirCult$strand)=="+"], labels = kirCult$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(GraphCULT, data = kirCult[strand(kirCult$strand)=="-"], labels = kirCult$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 5)
  #sicavars
  kpPlotMarkers(GraphCULT, data = SICACult[strand(SICACult$strand)=="+"], labels = SICACult$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(GraphCULT, data = SICACult[strand(SICACult$strand)=="-"], labels = SICACult$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  dev.off()
#}
  
  #plot050 <- function() {
  #png("sks050SICAvars.png", width = 1600, height = 800)
  tiff("sks050SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
  pp$topmargin <- 450
  pp$bottommargin <- 400
  Graph050 <- plotKaryotype(genome = sks050, ideogram.plotter = NULL, 
                            plot.type = 2, plot.params = pp, cex=1.0)
  kpAddCytobandsAsLine(Graph050)
  kpAddMainTitle(Graph050, "sks050", cex=1.2)
  kpAddBaseNumbers(Graph050, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph050, data=Genes050[strand(Genes050)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph050, data=Genes050[strand(Genes050)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph050, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph050, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(Graph050, data = kir050[strand(kir050$strand)=="+"], labels = kir050$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph050, data = kir050[strand(kir050$strand)=="-"], labels = kir050$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph050, data = SICA050[strand(SICA050$strand)=="+"], labels = SICA050$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph050, data = SICA050[strand(SICA050$strand)=="-"], labels = SICA050$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
  
  #plot058 <- function() {
  #png("sks058SICAvars.png", width = 1600, height = 800)
  tiff("sks058SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(Graph058, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph058, data=Genes058[strand(Genes058)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph058, data=Genes058[strand(Genes058)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph058, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph058, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(Graph058, data = kir058[strand(kir058$strand)=="+"], labels = kir058$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph058, data = kir058[strand(kir058$strand)=="-"], labels = kir058$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph058, data = SICA058[strand(SICA058$strand)=="+"], labels = SICA058$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph058, data = SICA058[strand(SICA058$strand)=="-"], labels = SICA058$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
  
  #plot070 <- function() {
  #png("sks070SICAvars.png", width = 1600, height = 800)
  tiff("sks070SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(Graph070, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph070, data=Genes070[strand(Genes070)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph070, data=Genes070[strand(Genes070)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph070, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph070, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(Graph070, data = kir070[strand(kir070$strand)=="+"], labels = kir070$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph070, data = kir070[strand(kir070$strand)=="-"], labels = kir070$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph070, data = SICA070[strand(SICA070$strand)=="+"], labels = SICA070$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph070, data = SICA070[strand(SICA070$strand)=="-"], labels = SICA070$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
  
  #plot074 <- function() {
  #png("sks074SICAvars.png", width = 1600, height = 800)
  tiff("sks074SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(Graph074, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph074, data=Genes074[strand(Genes074)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph074, data=Genes074[strand(Genes074)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph074, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph074, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(Graph074, data = kir074[strand(kir074$strand)=="+"], labels = kir074$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph074, data = kir074[strand(kir074$strand)=="-"], labels = kir074$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph074, data = SICA074[strand(SICA074$strand)=="+"], labels = SICA074$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph074, data = SICA074[strand(SICA074$strand)=="-"], labels = SICA074$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
  
  #plot078 <- function() {
  #png("sks078SICAvars.png", width = 1600, height = 800)
  tiff("sks078SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(Graph078, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph078, data=Genes078[strand(Genes078)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph078, data=Genes078[strand(Genes078)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph078, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph078, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(Graph078, data = kir078[strand(kir078$strand)=="+"], labels = kir078$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph078, data = kir078[strand(kir078$strand)=="-"], labels = kir078$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph078, data = SICA078[strand(SICA078$strand)=="+"], labels = SICA078$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph078, data = SICA078[strand(SICA078$strand)=="-"], labels = SICA078$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
  #plot125 <- function() {
  #png("sks125SICAvars.png", width = 1600, height = 800)
  tiff("sks125SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(Graph125, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph125, data=Genes125[strand(Genes125)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph125, data=Genes125[strand(Genes125)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph125, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph125, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(Graph125, data = kir125[strand(kir125$strand)=="+"], labels = kir125$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph125, data = kir125[strand(kir125$strand)=="-"], labels = kir125$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph125, data = SICA125[strand(SICA125$strand)=="+"], labels = SICA125$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph125, data = SICA125[strand(SICA125$strand)=="-"], labels = SICA125$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
  
  #plot325 <- function() {
  #png("sks325SICAvars.png", width = 1600, height = 800)
  tiff("sks325SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(Graph325, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph325, data=Genes325[strand(Genes325)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph325, data=Genes325[strand(Genes325)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph325, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph325, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(Graph325, data = kir325[strand(kir325$strand)=="+"], labels = kir325$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph325, data = kir325[strand(kir325$strand)=="-"], labels = kir325$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph325, data = SICA325[strand(SICA325$strand)=="+"], labels = SICA325$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph325, data = SICA325[strand(SICA325$strand)=="-"], labels = SICA325$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
  
  #plot331 <- function() {
  #png("sks331SICAvars.png", width = 1600, height = 800)
  tiff("sks331SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(Graph331, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph331, data=Genes331[strand(Genes331)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph331, data=Genes331[strand(Genes331)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph331, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph331, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(Graph331, data = kir331[strand(kir331$strand)=="+"], labels = kir331$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph331, data = kir331[strand(kir331$strand)=="-"], labels = kir331$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph331, data = SICA331[strand(SICA331$strand)=="+"], labels = SICA331$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph331, data = SICA331[strand(SICA331$strand)=="-"], labels = SICA331$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
  
  
  #plot333 <- function() {
  #png("sks333SICAvars.png", width = 1600, height = 800)
  tiff("sks333SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(Graph333, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph333, data=Genes333[strand(Genes333)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph333, data=Genes333[strand(Genes333)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph333, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph333, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(Graph333, data = kir333[strand(kir333$strand)=="+"], labels = kir333$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph333, data = kir333[strand(kir333$strand)=="-"], labels = kir333$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph333, data = SICA333[strand(SICA333$strand)=="+"], labels = SICA333$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph333, data = SICA333[strand(SICA333$strand)=="-"], labels = SICA333$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
  
  #plot339 <- function() {
  #png("sks339SICAvars.png", width = 1600, height = 800)
  tiff("sks339SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(Graph339, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph339, data=Genes339[strand(Genes339)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph339, data=Genes339[strand(Genes339)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph339, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph339, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(Graph339, data = kir339[strand(kir339$strand)=="+"], labels = kir339$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph339, data = kir339[strand(kir339$strand)=="-"], labels = kir339$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph339, data = SICA339[strand(SICA339$strand)=="+"], labels = SICA339$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph339, data = SICA339[strand(SICA339$strand)=="-"], labels = SICA339$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
  
  #plot344 <- function() {
  #png("sks344SICAvars.png", width = 1600, height = 800)
  tiff("sks344SICA_KIRs.tif", units="px", width=6400, height=3200, res=300)
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
  kpAddBaseNumbers(Graph344, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350, cex=1)
  kpPlotRegions(Graph344, data=Genes344[strand(Genes344)=="+"], avoid.overlapping = FALSE)
  kpPlotRegions(Graph344, data=Genes344[strand(Genes344)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  kpAddLabels(Graph344, "strand +", cex=0.75, col="#888888", r1=3)
  kpAddLabels(Graph344, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  #kirs
  kpPlotMarkers(Graph344, data = kir344[strand(kir344$strand)=="+"], labels = kir344$gene, label.color = "#a90288",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "#a90288", adjust.label.position = 5)
  kpPlotMarkers(Graph344, data = kir344[strand(kir344$strand)=="-"], labels = kir344$gene, label.color = "#471dcf",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#471dcf", adjust.label.position = 1)
  #sicavars
  kpPlotMarkers(Graph344, data = SICA344[strand(SICA344$strand)=="+"], labels = SICA344$gene, label.color = "#2fc058",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "#2fc058", ignore.chromosome.ends = TRUE)
  kpPlotMarkers(Graph344, data = SICA344[strand(SICA344$strand)=="-"], labels = SICA344$gene, label.color = "#db7409",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "#db7409")
  
  dev.off()
  
  rm(list = ls())
####################################################################################################################################################
# png("SICAvar_and_KIRs.png", width = 2400, height = 1200)
# p1 <- as.ggplot(expression(plot047()))
# p2 <- as.ggplot(expression(plot048()))
# p3 <- as.ggplot(expression(plotCULT()))
# p4 <- as.ggplot(expression(plotPKNH()))
# plot_grid(p4, p3, p1, p2, ncol=2)
# dev.off()

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
  
  # sks048 for the biomalpar poster
  # png("CultSICAvars.png", width = 3200, height = 1600)
  # pp <- getDefaultPlotParams(plot.type = 2)
  # # set boundaries for the graph
  # pp$leftmargin <- 0.1
  # pp$data1outmargin <- 350
  # pp$data2outmargin <- 1000
  # pp$topmargin <- 450
  # pp$bottommargin <- 400
  # GraphPKNH <- plotKaryotype(genome = PKNHy, ideogram.plotter = NULL, 
  #                           plot.type = 2, plot.params = pp, cex=2)
  # kpAddCytobandsAsLine(GraphPKNH)
  # kpAddMainTitle(GraphPKNH, "PKNH", cex=2)
  # kpAddBaseNumbers(GraphPKNH, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350,cex = 1.5)
  # kpPlotRegions(GraphPKNH, data=GenesPKNH[strand(GenesPKNH)=="+"], avoid.overlapping = FALSE)
  # kpPlotRegions(GraphPKNH, data=GenesPKNH[strand(GenesPKNH)=="-"], avoid.overlapping = FALSE, data.panel = 2)
  # kpAddLabels(GraphPKNH, "strand +", cex=1.5, col="#888888", r1=3)
  # kpAddLabels(GraphPKNH, "strand -", cex=1.5, data.panel = 2, col="#888888", r1=3)
  # #kirs
  # # kpPlotMarkers(Graph048, data = kir048[strand(kir048$strand)=="+"], labels = kir048$gene, label.color = "magenta",
  # #               text.orientation = "horizontal", r1=1.5, cex=1.8, marker.parts = c(0.2, 0.7, 0.1),
  # #               data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "magenta", adjust.label.position = 5)
  # # kpPlotMarkers(Graph048, data = kir048[strand(kir048$strand)=="-"], labels = kir048$gene, label.color = "turquoise1",
  # #               text.orientation = "horizontal", r1=1.5, cex=1.8, marker.parts = c(0.2, 0.7, 0.1),
  # #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "turquoise1", adjust.label.position = 1)
  # #sicavars
  # kpPlotMarkers(GraphPKNH, data = SICAPKNH[strand(SICAPKNH$strand)=="+"], labels = SICAPKNH$gene, label.color = "darkred",
  #               text.orientation = "horizontal", r1=2.5, cex=1.5, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 1, label.dist = 0.00001, max.iter = 1000, line.color = "darkred", ignore.chromosome.ends = TRUE)
  # kpPlotMarkers(GraphPKNH, data = SICAPKNH[strand(SICAPKNH$strand)=="-"], labels = SICAPKNH$gene, label.color = "seagreen",
  #               text.orientation = "horizontal", r1=2.5, cex=1.5, marker.parts = c(0.2, 0.7, 0.1),
  #               data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "seagreen")
  # 
  # dev.off()