# clear the environment
rm(list = ls())
# use tidyverse package
library("tidyverse")
# load in the FRC text files from FRC bam otuput
df <- read.table("sks047_FRC.txt")
df2 <- read.table("sks048_FRC.txt")
df3 <- read.table("sks058_FRC.txt")
dfp1 <- read.table("CANUsks047_FRC.txt")
dfp2 <- read.table("CANUsks048_FRC.txt")
dfp3 <- read.table("CANUsks058_FRC.txt")
# change the column names of the variables to be threshold and coverage
colnames(df) <- c("Feature_Threshold", "Approximate_Coverage_%")
colnames(df2) <- c("Feature_Threshold", "Approximate_Coverage_%")
colnames(df3) <- c("Feature_Threshold", "Approximate_Coverage_%")
colnames(dfp1) <- c("Feature_Threshold", "Approximate_Coverage_%")
colnames(dfp2) <- c("Feature_Threshold", "Approximate_Coverage_%")
colnames(dfp3) <- c("Feature_Threshold", "Approximate_Coverage_%")
# check the structure of the variable
str(dfp1)
#unlist the values of the dataframe to make them available for use by ggplot2
df$Feature_Threshold <- unlist(df$Feature_Threshold)
df2$Feature_Threshold <- unlist(df2$Feature_Threshold)
df3$Feature_Threshold <- unlist(df3$Feature_Threshold)
dfp1$Feature_Threshold <- unlist(dfp1$Feature_Threshold)
dfp2$Feature_Threshold <- unlist(dfp2$Feature_Threshold)
dfp3$Feature_Threshold <- unlist(dfp3$Feature_Threshold)
df$`Approximate_Coverage_%` <- unlist(df$`Approximate_Coverage_%`)
df2$`Approximate_Coverage_%` <- unlist(df2$`Approximate_Coverage_%`)
df3$`Approximate_Coverage_%` <- unlist(df3$`Approximate_Coverage_%`)
dfp1$`Approximate_Coverage_%` <- unlist(dfp1$`Approximate_Coverage_%`)
dfp2$`Approximate_Coverage_%` <- unlist(dfp2$`Approximate_Coverage_%`)
dfp3$`Approximate_Coverage_%` <- unlist(dfp3$`Approximate_Coverage_%`)

# check the structure again to ensure nothing is affected
str(dfp1)
#set the isolate IDs
isolate <- c("sks047", "sks048", "sks058", "Canusks047", "Canusks048", "Canusks058")
# add a column of the isolate name
df <- mutate(df, isolate[1])
colnames(df)[3] <- c("Isolate_ID")

df2 <- mutate(df2, isolate[2])
colnames(df2)[3] <- c("Isolate_ID")

df3 <- mutate(df3, isolate[3])
colnames(df3)[3] <- c("Isolate_ID")

dfp1 <- mutate(dfp1, isolate[4])
colnames(dfp1)[3] <- c("Isolate_ID")

dfp2 <- mutate(dfp2, isolate[5])
colnames(dfp2)[3] <- c("Isolate_ID")

dfp3 <- mutate(dfp3, isolate[6])
colnames(dfp3)[3] <- c("Isolate_ID")

#plot the graph using the features as the x axis and coverage as the y axis, colours will be based on isolate ID 
p <- ggplot(df, aes(x = df$Feature_Threshold, y = df$`Approximate_Coverage_%`, color = df$Isolate_ID)) +
  geom_step() +
  geom_step(data = df2, aes(x = df2$Feature_Threshold, y = df2$`Approximate_Coverage_%`, colour = df2$Isolate_ID)) +
  geom_step(data = df3, aes(x = df3$Feature_Threshold, y = df3$`Approximate_Coverage_%`, colour = df3$Isolate_ID)) +
  geom_step(data = dfp1, aes(x = dfp1$Feature_Threshold, y = dfp1$`Approximate_Coverage_%`, colour = dfp1$Isolate_ID)) +
  geom_step(data = dfp2, aes(x = dfp2$Feature_Threshold, y = dfp2$`Approximate_Coverage_%`, colour = dfp2$Isolate_ID)) +
  geom_step(data = dfp3, aes(x = dfp3$Feature_Threshold, y = dfp3$`Approximate_Coverage_%`, colour = dfp3$Isolate_ID)) +
  labs(title = 'FRCurve of features'
       ,y = 'Approximate Coverage %'
       ,x = 'Feature Threshold'
       ,subtitle = str_c("showing fragmentation of the genomes")
       ,caption = "amazing"
       ,color = "Isolates:"
  ) +
  #scale_color_manual(labels = isolate, values = c("blue", "red" , "orange", "green", "brown", "purple")) +
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
ggsave("testing.png")
length(dfp1$Feature_Threshold)
