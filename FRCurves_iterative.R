##INFO: Please set the working directory to the folder containing the output of the FRCbam. 
##Be aware that the only the text files should have _FRC.txt in their filenames.

# clear the environment
rm(list = ls())
# use tidyverse package
library("tidyverse")
# load in the FRC text files from FRC bam outputput
filelist = list.files(pattern = "_FRC*.txt")
# load the files into a list of dataframes
datalist = lapply(filelist, function(x)read.table(x, header=T))
# set the column names
colnames <- c("Feature_Threshold", "Approximate_Coverage_%")
# change the column names of all the dataframes in the list
clr <- lapply(datalist, setNames, colnames)

# for each dataframe in the list:
for (i in 1:length(clr)) {
  print(i)
  print(clr[[i]]$Feature_Threshold)
  #unlist the values of the dataframe to make them available for use by ggplot2
  clr[[i]]$Feature_Threshold <- unlist(clr[[i]]$Feature_Threshold)
  clr[[i]]$`Approximate_Coverage_%` <- unlist(clr[[i]]$`Approximate_Coverage_%`)
  #set the isolate IDs
  isolate <- c("CANUsks047", "CANUsks048", "CANUsks058", "sks047", "sks048", "sks058")
  # add a column of the isolate name
  clr[[i]] <- mutate(clr[[i]], isolate[i])
  colnames(clr[[i]])[3] <- c("Isolate_ID")
  print(clr[[i]])
  
}
# set the colours of each isolate
tcolr <- c("blue", "red" , "orange", "green", "brown", "purple")

#plot the graph using the features as the x axis and coverage as the y axis, colours will be based on isolate ID 
p <- ggplot(clr[[1]], aes(x = clr[[1]]$Feature_Threshold, y = clr[[1]]$`Approximate_Coverage_%`, color = clr[[1]]$Isolate_ID)) +
  geom_step() +
  geom_step(data = clr[[2]], aes(x = clr[[2]]$Feature_Threshold, y = clr[[2]]$`Approximate_Coverage_%`, colour = clr[[2]]$Isolate_ID)) +
  geom_step(data = clr[[3]], aes(x = clr[[3]]$Feature_Threshold, y = clr[[3]]$`Approximate_Coverage_%`, colour = clr[[3]]$Isolate_ID)) +
  #geom_step(data = clr[[4]], aes(x = clr[[4]]$Feature_Threshold, y = clr[[4]]$`Approximate_Coverage_%`, colour = clr[[4]]$Isolate_ID)) +
  #geom_step(data = clr[[5]], aes(x = clr[[5]]$Feature_Threshold, y = clr[[5]]$`Approximate_Coverage_%`, colour = clr[[5]]$Isolate_ID)) +
  #geom_step(data = clr[[6]], aes(x = clr[[6]]$Feature_Threshold, y = clr[[6]]$`Approximate_Coverage_%`, colour = clr[[6]]$Isolate_ID)) +
  
  # set the title of the pplot, the subtitle, caption and figure legend
  labs(title = 'FRCurve of features'
       ,y = 'Approximate Coverage %'
       ,x = 'Feature Threshold'
       ,subtitle = str_c("showing fragmentation of the genomes")
       ,caption = "FRCurve comparison of the Pilon corrected isolates within the dataset between the two available assemblers -- Canu and Flye. 
       The steeper the line of the plot, the better the assembler is performing. Hence, across all three isolates, 
       the Canu assembler appears to be performing better"
       ,color = "Isolates:"
  ) +
  #set the colour of each isolate based on the colours set previously
  scale_color_manual(labels = isolate, values = tcolr) +
  # set the theme parameters, font size, text etc
  theme_bw() +
  theme(text = element_text(color = "#444444", family = 'Arial')
        ,plot.title = element_text(size = 18, color = '#333333')
        ,plot.subtitle = element_text(size = 10)
        ,axis.title = element_text(size = 10, color = '#333333')
        ,axis.title.y = element_text(angle = 90, vjust = 0)
        ,legend.position = "bottom"
  )
# save to a PNG file
ggsave("FRCurve_Comparisons.png")
