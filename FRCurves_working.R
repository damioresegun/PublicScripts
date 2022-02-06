# clear the environment
rm(list = ls()
# use tidyverse package
library("tidyverse")
library("optparse")
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Text file containing FRCbam output", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out_graph.png", 
              help="Output file name in PNG [default= %default]", metavar="character"),
  make_option(c("-p", "--preset"), type="character", default="isolate", 
              help="Preset name for each of the isolates being entered", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt <- "/storage/home/users/dro/FRCbam/sks047/frcs/sks047_FRC.txt,/storage/home/users/dro/FRCbam/CANUsks047/sks047_FRC.txt"
# if (is.null(opt$file)){
#   print_help(opt_parser)
#   stop("At least the full path to the input and output files are needed", call.=FALSE)
# }
filelist <- strsplit(opt,"[,:]")
filelisty <- as.data.frame(filelist)
frt <- as.matrix(filelisty)
temp.file <- paste(filelist, filelisty, sep=",")
#temp <- read.csv(temp.file)
print(filelisty[1])
print(frt[2])
print(temp.file)

# load in the FRC text files from FRC bam otuput
df <- read.table("sks047_FRC.txt")
dfp <- read.table("Canusks047_FRC.txt")
# change the column names of the variables to be threshold and coverage
colnames(df) <- c("Feature_Threshold", "Approximate_Coverage_%")
colnames(dfp) <- c("Feature_Threshold", "Approximate_Coverage_%")
# check the structure of the variable
str(dfp)
#unlist the values of the dataframe to make them available for use by ggplot2
df$Feature_Threshold <- unlist(df$Feature_Threshold)
dfp$Feature_Threshold <- unlist(dfp$Feature_Threshold)
df$`Approximate_Coverage_%` <- unlist(df$`Approximate_Coverage_%`)
dfp$`Approximate_Coverage_%` <- unlist(dfp$`Approximate_Coverage_%`)
# check the structure again to ensure nothing is affected
str(dfp)
#set the isolate IDs
isolate <- c("sks047", "Canusks047")
# add a column of the isolate name
df <- mutate(df, isolate[1])
colnames(df)[3] <- c("Isolate_ID")
dfp <- mutate(dfp, isolate[2])
colnames(dfp)[3] <- c("Isolate_ID")
#plot the graph using the features as the x axis and coverage as the y axis, colours will be based on isolate ID 
p <- ggplot(df, aes(x = df$Feature_Threshold, y = df$`Approximate_Coverage_%`, color = df$Isolate_ID)) +
  geom_step() +
  geom_step(data = dfp, aes(x = dfp$Feature_Threshold, y = dfp$`Approximate_Coverage_%`, colour = dfp$Isolate_ID)) +
  labs(title = 'FRCurve of features'
       ,y = 'Approximate Coverage %'
       ,x = 'Feature Threshold'
       ,subtitle = str_c("showing fragmentation of the genomes")
       ,caption = "amazing"
       ,color = "Isolates:"
  ) +
  scale_color_manual(labels = isolate, values = c("blue", "red")) +
  #theme_bw() +
  theme(text = element_text(color = "#444444", family = 'Arial')
        ,plot.title = element_text(size = 18, color = '#333333')
        ,plot.subtitle = element_text(size = 10)
        ,axis.title = element_text(size = 10, color = '#333333')
        ,axis.title.y = element_text(angle = 90, vjust = 0)
        ,legend.position = "bottom"
  )
#p + scale_linetype_discrete(name = "Dose (mg)")
print(p)
