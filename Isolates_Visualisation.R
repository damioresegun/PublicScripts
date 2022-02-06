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
#BiocManager::install("cowplot")
library(karyoploteR)
library(rtracklayer)
library(ggplotify)
library(cowplot)
# load in your gff file
pknhGFF <- "Pknowlesi_pseudo.out.gff3"
sks047GFF <- "sks047_pseudo.out.gff3"
sks048GFF <- "sks048_pseudo.out.gff3"
cultGFF <- "Cultured_pseudo.out.gff3"
# load in the gene files
sks047IN_MSPs <- "sks047_MSPs.gff3"
sks048IN_MSPs <- "sks048_MSPs.gff3"
CultIN_MSPs <- "Cultured_MSPs.gff3"
pknhIN_MSPs <- "PKNH_MSPs.gff3"
sks047IN_CIRC <- "sks047_circum.gff3"
sks048IN_CIRC <- "sks048_circum.gff3"
CultIN_CIRC <- "Cultured_circum.gff3"
pknhIN_CIRC <- "PKNH_circum.gff3"
sks047IN_ERY <- "sks047_erythro.gff3"
sks048IN_ERY <- "sks048_erythro.gff3"
CultIN_ERY <- "Cultured_erythro.gff3"
pknhIN_ERY <- "PKNH_erythro.gff3"
sks047IN_SPOR <- "sks047_sporo.gff3"
sks048IN_SPOR <- "sks048_sporo.gff3"
CultIN_SPOR <- "Cultured_sporo.gff3"
pknhIN_SPOR <- "PKNH_sporo.gff3"
sks047IN_pknbp <- "sks047_pknbp.gff3"
sks048IN_pknbp <- "sks048_pknbp.gff3"
CultIN_pknbp <- "Cultured_pknbp.gff3"
pknhIN_pknbp <- "PKNH_pknbp.gff3"
sks047IN_kahrp <- "sks047_KAHRP.gff3"
sks048IN_kahrp <- "sks048_KAHRP.gff3"
CultIN_kahrp <- "Cultured_KAHRP.gff3"
pknhIN_kahrp <- "PKNH_KAHRP.gff3"
sks047IN_KIR <- "sks047_KIRs.gff3"
sks048IN_KIR <- "sks048_KIRs.gff3"
CultIN_KIR <- "Cultured_KIRs.gff3"
pknhIN_KIR <- "PKNH_KIRs.gff3"
sks047IN_MDR <- "sks047_MDR.gff3"
sks048IN_MDR <- "sks048_MDR.gff3"
pknhIN_MDR <- "PKNH_MDR.gff3"
CultIN_MDR <- "Cultured_MDR.gff3"
sks047IN_DBP <- "sks047_DBP.gff3"
sks048IN_DBP <- "sks048_DBP.gff3"
pknhIN_DBP <- "PKNH_DBP.gff3"
CultIN_DBP <- "Cultured_DBP.gff3"
sks047IN_CyAP <- "sks047_CyAP.gff3"
sks048IN_CyAP <- "sks048_CyAP.gff3"
pknhIN_CyAP <- "PKNH_CyAP.gff3"
CultIN_CyAP <- "Cultured_CyAP.gff3"
sks047IN_TrpRA <- "sks047_TrpRA.gff3"
sks048IN_TrpRA <- "sks048_TrpRA.gff3"
pknhIN_TrpRA <- "PKNH_TrpRA.gff3"
CultIN_TrpRA <- "Cultured_TrpRA.gff3"
sks047IN_ETRAMP <- "sks047_ETRAMP.gff3"
sks048IN_ETRAMP <- "sks048_ETRAMP.gff3"
pknhIN_ETRAMP <- "PKNH_ETRAMP.gff3"
CultIN_ETRAMP <- "Cultured_ETRAMP.gff3"
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
PKNH <- toGRanges(ffPKNH[,c(4,5,6)])
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
####### add in the  erythrocyte binding protein genes ########
ery047IN <- read.delim(sks047IN_ERY)
ery047 <- toGRanges(data.frame(chromosome = c(ery047IN$chromosome), start=c(ery047IN$start),
                                end=c(ery047IN$end), gene = c(ery047IN$short), strand=c(ery047IN$strand)))
ery047 <- sort(ery047)
##
ery048IN <- read.delim(sks048IN_ERY)
ery048 <- toGRanges(data.frame(chromosome = c(ery048IN$chromosome), start=c(ery048IN$start),
                                end=c(ery048IN$end), gene = c(ery048IN$short), strand=c(ery048IN$strand)))
ery048 <- sort(ery048)
##
eryCultIN <- read.delim(CultIN_ERY)
eryCult <- toGRanges(data.frame(chromosome = c(eryCultIN$chromosome), start=c(eryCultIN$start),
                                 end=c(eryCultIN$end), gene = c(eryCultIN$short), strand=c(eryCultIN$strand)))
eryCult <- sort(eryCult)
##
eryPKNHIN <- read.delim(pknhIN_ERY)
eryPKNH <- toGRanges(data.frame(chromosome = c(eryPKNHIN$chromosome), start=c(eryPKNHIN$start),
                                 end=c(eryPKNHIN$end), gene = c(eryPKNHIN$short), strand=c(eryPKNHIN$strand)))
eryPKNH <- sort(eryPKNH)

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

####### add in the multidrug resistance genes ########
mdr047IN <- read.delim(sks047IN_MDR)
mdr047 <- toGRanges(data.frame(chromosome = c(mdr047IN$chromosome), start=c(mdr047IN$start),
                                 end=c(mdr047IN$end), gene = c(mdr047IN$short), strand=c(mdr047IN$strand)))
mdr047 <- sort(mdr047)
##
mdr048IN <- read.delim(sks048IN_MDR)
mdr048 <- toGRanges(data.frame(chromosome = c(mdr048IN$chromosome), start=c(mdr048IN$start),
                                 end=c(mdr048IN$end), gene = c(mdr048IN$short), strand=c(mdr048IN$strand)))
mdr048 <- sort(mdr048)
##
mdrCultIN <- read.delim(CultIN_MDR)
mdrCult <- toGRanges(data.frame(chromosome = c(mdrCultIN$chromosome), start=c(mdrCultIN$start),
                                  end=c(mdrCultIN$end), gene = c(mdrCultIN$short), strand=c(mdrCultIN$strand)))
mdrCult <- sort(mdrCult)
##
mdrPKNHIN <- read.delim(pknhIN_MDR)
mdrPKNH <- toGRanges(data.frame(chromosome = c(mdrPKNHIN$chromosome), start=c(mdrPKNHIN$start),
                                  end=c(mdrPKNHIN$end), gene = c(mdrPKNHIN$short), strand=c(mdrPKNHIN$strand)))
mdrPKNH <- sort(mdrPKNH)

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
kirPKNHIN <- read.delim(pknhIN_KIR)
kirPKNH <- toGRanges(data.frame(chromosome = c(kirPKNHIN$chromosome), start=c(kirPKNHIN$start),
                                end=c(kirPKNHIN$end), gene = c(kirPKNHIN$short), strand=c(kirPKNHIN$strand)))
kirPKNH <- sort(kirPKNH)
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

####### add in the tryptophan genes ########
trpRA047IN <- read.delim(sks047IN_TrpRA)
trpRA047 <- toGRanges(data.frame(chromosome = c(trpRA047IN$chromosome), start=c(trpRA047IN$start),
                                  end=c(trpRA047IN$end), gene = c(trpRA047IN$short), strand=c(trpRA047IN$strand)))
trpRA047 <- sort(trpRA047)
##
trpRA048IN <- read.delim(sks048IN_TrpRA)
trpRA048 <- toGRanges(data.frame(chromosome = c(trpRA048IN$chromosome), start=c(trpRA048IN$start),
                                  end=c(trpRA048IN$end), gene = c(trpRA048IN$short), strand=c(trpRA048IN$strand)))
trpRA048 <- sort(trpRA048)
##
trpRACultIN <- read.delim(CultIN_TrpRA)
trpRACult <- toGRanges(data.frame(chromosome = c(trpRACultIN$chromosome), start=c(trpRACultIN$start),
                                   end=c(trpRACultIN$end), gene = c(trpRACultIN$short), strand=c(trpRACultIN$strand)))
trpRACult <- sort(trpRACult)
##
trpRAPKNHIN <- read.delim(pknhIN_TrpRA)
trpRAPKNH <- toGRanges(data.frame(chromosome = c(trpRAPKNHIN$chromosome), start=c(trpRAPKNHIN$start),
                                   end=c(trpRAPKNHIN$end), gene = c(trpRAPKNHIN$short), strand=c(trpRAPKNHIN$strand)))
trpRAPKNH <- sort(trpRAPKNH)
############################################################################################################
# plot all the genes on the chromosomes and save to file
#plot047 <- function() {
  png("sks047.png", width = 1600, height = 800)
  pp <- getDefaultPlotParams(plot.type = 2)
  # set boundaries for the graph
  pp$leftmargin <- 0.12
  pp$data1outmargin <- 350
  pp$data2outmargin <- 1000
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
  kpAddLabels(Graph047, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3)
  # MSP genes
  kpPlotMarkers(Graph047, data = msp047[strand(msp047$strand)=="+"], labels = msp047$gene, label.color = "purple",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
  kpPlotMarkers(Graph047, data = msp047[strand(msp047$strand)=="-"], labels = msp047$gene, label.color = "brown",
              text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
  #circumsporozoite
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
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.2, line.color = "darkorange1")
  # erythrocyte
  kpPlotMarkers(Graph047, data = ery047[strand(ery047$strand)=="+"], labels = ery047$gene, label.color = "firebrick3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "firebrick3", offset = 0.5)
  kpPlotMarkers(Graph047, data = ery047[strand(ery047$strand)=="-"], labels = ery047$gene, label.color = "gold4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.9, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 1, line.color = "gold4")
  # pknbp
  kpPlotMarkers(Graph047, data = pknb047[strand(pknb047$strand)=="+"], labels = pknb047$gene, label.color = "steelblue3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "steelblue3")
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
  #mdr
  kpPlotMarkers(Graph047, data = mdr047[strand(mdr047$strand)=="+"], labels = mdr047$gene, label.color = "royalblue1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "royalblue1")
  kpPlotMarkers(Graph047, data = mdr047[strand(mdr047$strand)=="-"], labels = mdr047$gene, label.color = "sienna1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.2, offset = 0.3, line.color = "sienna1")
  #kirs
  kpPlotMarkers(Graph047, data = kir047[strand(kir047$strand)=="+"], labels = kir047$gene, label.color = "magenta",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "magenta", adjust.label.position = 5)
  kpPlotMarkers(Graph047, data = kir047[strand(kir047$strand)=="-"], labels = kir047$gene, label.color = "turquoise1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "turquoise1", adjust.label.position = 5)
  #trpra
  kpPlotMarkers(Graph047, data = trpRA047[strand(trpRA047$strand)=="+"], labels = trpRA047$gene, label.color = "aquamarine3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "aquamarine3", adjust.label.position = 5)
  kpPlotMarkers(Graph047, data = trpRA047[strand(trpRA047$strand)=="-"], labels = trpRA047$gene, label.color = "coral2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "coral2", adjust.label.position = 5)
  #etramp
  kpPlotMarkers(Graph047, data = etramp047[strand(etramp047$strand)=="+"], labels = etramp047$gene, label.color = "darkgoldenrod2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkgoldenrod2", adjust.label.position = 5)
  kpPlotMarkers(Graph047, data = etramp047[strand(etramp047$strand)=="-"], labels = etramp047$gene, label.color = "darkolivegreen2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "darkolivegreen2", adjust.label.position = 5)
  #dbp
  kpPlotMarkers(Graph047, data = dbp047[strand(dbp047$strand)=="+"], labels = dbp047$gene, label.color = "darkred",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "darkred", adjust.label.position = 5)
  kpPlotMarkers(Graph047, data = dbp047[strand(dbp047$strand)=="-"], labels = dbp047$gene, label.color = "mediumorchid2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "mediumorchid2", adjust.label.position = 5)
  #clag
  kpPlotMarkers(Graph047, data = clag047[strand(clag047$strand)=="-"], labels = clag047$gene, label.color = "lightsteelblue2",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "lightsteelblue2", adjust.label.position = 5)
  
  dev.off()
#  }

plot048 <- function() {
 # png("sks048.png", width = 1600, height = 800)
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
  #circumsporozoite
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
  # erythrocyte
  kpPlotMarkers(Graph048, data = ery048[strand(ery048$strand)=="+"], labels = ery048$gene, label.color = "firebrick3",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "firebrick3", offset = 0.5)
  kpPlotMarkers(Graph048, data = ery048[strand(ery048$strand)=="-"], labels = ery048$gene, label.color = "gold4",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.9, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 1, line.color = "gold4")
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
  #mdr
  kpPlotMarkers(Graph048, data = mdr048[strand(mdr048$strand)=="+"], labels = mdr048$gene, label.color = "royalblue1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "royalblue1")
  kpPlotMarkers(Graph048, data = mdr048[strand(mdr048$strand)=="-"], labels = mdr048$gene, label.color = "sienna1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.2, offset = 0.3, line.color = "sienna1")
  #kirs
  kpPlotMarkers(Graph048, data = kir048[strand(kir048$strand)=="+"], labels = kir048$gene, label.color = "magenta",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "magenta", adjust.label.position = 5)
  kpPlotMarkers(Graph048, data = kir048[strand(kir048$strand)=="-"], labels = kir048$gene, label.color = "turquoise1",
                text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
                data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "turquoise1", adjust.label.position = 1)
#dev.off()
  }

plotPKNH <- function() {
 # png("PKNH.png", width = 1600, height = 800)
pp <- getDefaultPlotParams(plot.type = 2)
# set boundaries for the graph
pp$leftmargin <- 0.1
pp$data1outmargin <- 350
pp$data2outmargin <- 1000
pp$topmargin <- 450
pp$bottommargin <- 400
GraphPKNH <- plotKaryotype(genome = PKNH, ideogram.plotter = NULL, 
                          plot.type = 2, plot.params = pp, cex=1.0)
kpAddCytobandsAsLine(GraphPKNH)
kpAddMainTitle(GraphPKNH, "PKNH", cex=1.2)
kpAddBaseNumbers(GraphPKNH, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
kpPlotRegions(GraphPKNH, data=GenesPKNH[strand(GenesPKNH)=="+"], avoid.overlapping = FALSE)
kpPlotRegions(GraphPKNH, data=GenesPKNH[strand(GenesPKNH)=="-"], avoid.overlapping = FALSE, data.panel = 2)
kpAddLabels(GraphPKNH, "strand +", cex=0.75, col="#888888", r1=3.2)
kpAddLabels(GraphPKNH, "strand -", cex=0.75, data.panel = 2, col="#888888", r1=3.2)
# MSP genes
kpPlotMarkers(GraphPKNH, data = mspPKNH[strand(mspPKNH$strand)=="+"], labels = mspPKNH$gene, label.color = "purple",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
kpPlotMarkers(GraphPKNH, data = mspPKNH[strand(mspPKNH$strand)=="-"], labels = mspPKNH$gene, label.color = "brown",
              text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
#circumsporozoite
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
# erythrocyte
kpPlotMarkers(GraphPKNH, data = eryPKNH[strand(eryPKNH$strand)=="+"], labels = eryPKNH$gene, label.color = "firebrick3",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "firebrick3", offset = 0.5)
kpPlotMarkers(GraphPKNH, data = eryPKNH[strand(eryPKNH$strand)=="-"], labels = eryPKNH$gene, label.color = "gold4",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.9, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 1, line.color = "gold4")
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
#mdr
kpPlotMarkers(GraphPKNH, data = mdrPKNH[strand(mdrPKNH$strand)=="+"], labels = mdrPKNH$gene, label.color = "royalblue1",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "royalblue1")
kpPlotMarkers(GraphPKNH, data = mdrPKNH[strand(mdrPKNH$strand)=="-"], labels = mdrPKNH$gene, label.color = "sienna1",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.2, offset = 0.3, line.color = "sienna1")
#kirs
kpPlotMarkers(GraphPKNH, data = kirPKNH[strand(kirPKNH$strand)=="+"], labels = kirPKNH$gene, label.color = "magenta",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "magenta", adjust.label.position = 5)
kpPlotMarkers(GraphPKNH, data = kirPKNH[strand(kirPKNH$strand)=="-"], labels = kirPKNH$gene, label.color = "turquoise1",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "turquoise1", adjust.label.position = 1)
#dev.off()
}

plotCULT <- function() {
 # png("Cultured_Pk.png", width = 1600, height = 800)
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
kpAddMainTitle(GraphCULT, "PkA1H1 Cultured", cex=1.2)
kpAddBaseNumbers(GraphCULT, tick.dist = 500000, add.units = TRUE,minor.tick.dist = 100000, tick.len = 600, minor.tick.len = 350)
kpPlotRegions(GraphCULT, data=Genes0CULT[strand(Genes0CULT)=="+"], avoid.overlapping = FALSE)
kpPlotRegions(GraphCULT, data=Genes0CULT[strand(Genes0CULT)=="-"], avoid.overlapping = FALSE, data.panel = 2)
kpAddLabels(GraphCULT, "strand +", cex=0.8, col="#888888", r1=3)
kpAddLabels(GraphCULT, "strand -", cex=0.8, data.panel = 2, col="#888888", r1=3)
# MSP genes
kpPlotMarkers(GraphCULT, data = mspCult[strand(mspCult$strand)=="+"], labels = mspCult$gene, label.color = "purple",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.0001, max.iter = 1000, line.color = "purple")
kpPlotMarkers(GraphCULT, data = mspCult[strand(mspCult$strand)=="-"], labels = mspCult$gene, label.color = "brown",
              text.orientation = "horizontal", r1=5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, line.color = "brown")
#circumsporozoite
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
# erythrocyte
kpPlotMarkers(GraphCULT, data = eryCult[strand(eryCult$strand)=="+"], labels = eryCult$gene, label.color = "firebrick3",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "firebrick3", offset = 0.5)
kpPlotMarkers(GraphCULT, data = eryCult[strand(eryCult$strand)=="-"], labels = eryCult$gene, label.color = "gold4",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.9, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 1, line.color = "gold4")
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
#mdr
kpPlotMarkers(GraphCULT, data = mdrCult[strand(mdrCult$strand)=="+"], labels = mdrCult$gene, label.color = "royalblue1",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.001, max.iter = 1000, line.color = "royalblue1")
kpPlotMarkers(GraphCULT, data = mdrCult[strand(mdrCult$strand)=="-"], labels = mdrCult$gene, label.color = "sienna1",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.2, offset = 0.3, line.color = "sienna1")
#kirs
kpPlotMarkers(GraphCULT, data = kirCult[strand(kirCult$strand)=="+"], labels = kirCult$gene, label.color = "magenta",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 1, label.dist = 0.00001, max.iter = 1000, offset = 0.2, line.color = "magenta", adjust.label.position = 5)
kpPlotMarkers(GraphCULT, data = kirCult[strand(kirCult$strand)=="-"], labels = kirCult$gene, label.color = "turquoise1",
              text.orientation = "horizontal", r1=2.5, cex=0.9, marker.parts = c(0.2, 0.7, 0.1),
              data.panel = 2, label.dist = 0.00001, max.iter = 1000, label.margin = 0.9, offset = 0.3, line.color = "turquoise1", adjust.label.position = 5)
#dev.off()
}

png("Comparison of genes in isolates.png", width = 2600, height = 1300)
p1 <- as.ggplot(expression(plot047()))
p2 <- as.ggplot(expression(plot048()))
p3 <- as.ggplot(expression(plotCULT()))
p4 <- as.ggplot(expression(plotPKNH()))
plot_grid(p4, p3, p1, p2, ncol=2)
dev.off()



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
