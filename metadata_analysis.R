### packages ####
library(readxl)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(plotly)
library(htmlwidgets)
library(lubridate)
library(readr)

### data ####
setwd('C:\\Users\\wikto\\Desktop\\wirusy')
metadata_demo <- read_xlsx('Nextclade pipeline demo v2.xlsx', sheet=1)
metadata <- read_tsv('metadata_filtered.tsv', col_names = F)

names <- metadata_demo[3,]
data <- metadata_demo[4:dim(metadata_demo)[1],]

colnames(data) <- names
colnames(metadata) <- colnames(data)


### date coversion ####
table(metadata$Variant)

#delete 4 genomes with before 22.10.2019
metadata_v2 <- metadata[order(metadata$`Collection date`),]
metadata_v3 <- metadata_v2[5:dim(metadata_v2)[1],]

### Converts character to date
test.date <- as.Date(metadata_v3$`Collection date`, format="%Y-%m-%d",tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))

#finding first date
start_date <- min(test.date, na.rm = T)

# Calculate the number of weeks since the start_date for each date
week_groups <- floor(as.numeric(difftime(test.date, start_date, units = "days")) / 7) + 1

# Combine dates with their respective week group
result <- data.frame(Date = test.date, WeekGroup = week_groups)

metadata_v3$week <- week_groups


### plot Total number of sequenced genomes per week ####

genomes_per_week <- data.frame(table(metadata_v3$week))
colnames(genomes_per_week) <- c("week", "freq")
genomes_per_week$week <- as.numeric(as.character(genomes_per_week$week))

png(paste("Total_number_of_sequenced_genomes_per_week_v5.png", sep = ""), width = 1000, height = 600)

p<-ggplot(genomes_per_week, aes(x=week, y=freq)) + 
  geom_bar(stat = "identity", fill = "steelblue", width = 1)+
  theme_classic()+
  scale_color_manual(values = "steelblue")+
  theme(plot.title = element_text(hjust=0.5, face="bold", size=17), 
        axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.text = element_text(size=14))+
  labs(title='', y="Total number of sequenced genomes", x= 'Week number')+
  scale_x_continuous(breaks = seq(0, max(as.numeric(genomes_per_week$week)), by = 10))

print(p)
dev.off()

### percent of VOC genomes ####
types <- c("Alpha", "Beta", "Delta", "Gamma", "Omicron")

metadata_VOC <- metadata_v3[grepl("VOC", metadata_v3$Variant),]
table(metadata_VOC$Variant)


data_to_plot <- as.data.frame(table(metadata_VOC$Variant, metadata_VOC$week))
colnames(data_to_plot) <- c("variant", "week", "freq")
data_to_plot$week <- as.numeric(as.character(data_to_plot$week))

library(RColorBrewer)
library(plotly)
colors <- c('#a6bddb','#74a9cf','#3690c0','#0570b0','#034e7b')

png(paste("Percent_of_VOC_genomes_per_week_v1.png", sep = ""), width = 1000, height = 600)

p <- ggplot(data_to_plot, aes(x = week, y=freq, fill=variant, color=variant)) + 
  geom_bar(position="fill", stat="identity")+
  theme_classic()+
  scale_color_manual(values=colors, guide = F)+
  #scale_color_blue(start = 0.8, end = 0.2, guide = F) +
  #scale_fill_grey(start = 0.8, end = 0.2, labels = types) +
  scale_fill_manual(values = colors,labels = types) + 
  labs(fill = "Variant")+
  theme(axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.text = element_text(size=14))+
  labs(title='', y="Percent of VOC genomes", x= 'Week number')
print(p)
dev.off()

### percent of VOI genomes ####
types <- c("Epsilon", "Eta", "Iota", "Kappa", "Lambda", "Mu", "Other", "Theta", "Zeta")

library("viridis") 

metadata_VOI <- metadata_v3[grepl("VOI", metadata_v3$Variant),]
table(metadata_VOI$Variant)

temp <- matrix(0, dim(metadata_VOI)[1])
for (i in 1:length(types)){
  temp[grepl(types[i],metadata_VOI$Variant)] <- types[i]
}

temp[temp==0] <- "Other"

metadata_VOI$types <- temp
table(metadata_VOI$types)

data_to_plot <- as.data.frame(table(metadata_VOI$types, metadata_VOI$week))
colnames(data_to_plot) <- c("variant", "week", "freq")
data_to_plot$week <- as.numeric(as.character(data_to_plot$week))

png(paste("Percent_of_VOI_genomes_per_week_v1.png", sep = ""), width = 1000, height = 600)

p <- ggplot(data_to_plot, aes(x = week, y=freq, fill=variant, color=variant)) + 
  geom_bar(position="fill", stat="identity")+
  theme_classic()+
  scale_color_brewer(palette = "RdBu", guide = F)+
  #scale_color_blue(start = 0.8, end = 0.2, guide = F) +
  #scale_fill_grey(start = 0.8, end = 0.2, labels = types) +
  scale_fill_brewer(palette = "RdBu",labels = types) + 
  labs(fill = "Variant")+
  theme(axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.text = element_text(size=14))+
  labs(title='', y="Percent of VOC genomes", x= 'Week number')
print(p)
dev.off()




### percent of all genomes ####
data_to_plot <- as.data.frame(table(metadata_v3$Variant, metadata_v3$week))
colnames(data_to_plot) <- c("variant", "week", "freq")
data_to_plot$week <- as.numeric(as.character(data_to_plot$week))

ggplot(data_to_plot, aes(x = week, y=freq, fill=variant)) + 
  geom_bar(position="fill", stat="identity")+
  theme_classic()+
  theme(legend.position = "none")

scale_fill_brewer(palette = "Set3") + 
  ylab("Percent of genomes") 


ggplotly(p)

### wykres voc,VOI, VUM, VOHC ####

types <- c("VOC", "VOI", "VUM", )

temp <- matrix(0, dim(metadata_v3)[1])
for (i in 1:length(types)){
  temp[grepl(types[i],metadata_v3$Variant)] <- types[i]
}
table(temp)

temp[temp==0] <- "Other"

metadata_v3$types <- temp
table(metadata_VOI$types)

data_to_plot <- as.data.frame(table(metadata_VOI$types, metadata_VOI$week))
colnames(data_to_plot) <- c("variant", "week", "freq")
data_to_plot$week <- as.numeric(as.character(data_to_plot$week))

png(paste("Percent_of_VOI_genomes_per_week_v1.png", sep = ""), width = 1000, height = 600)

p <- ggplot(data_to_plot, aes(x = week, y=freq, fill=variant, color=variant)) + 
  geom_bar(position="fill", stat="identity")+
  theme_classic()+
  scale_color_brewer(palette = "RdBu", guide = F)+
  #scale_color_blue(start = 0.8, end = 0.2, guide = F) +
  #scale_fill_grey(start = 0.8, end = 0.2, labels = types) +
  scale_fill_brewer(palette = "RdBu",labels = types) + 
  labs(fill = "Variant")+
  theme(axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),
        axis.text = element_text(size=14))+
  labs(title='', y="Percent of VOC genomes", x= 'Week number')
print(p)
dev.off()