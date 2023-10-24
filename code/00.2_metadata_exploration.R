#Author: Dzhansu Hasanova 10092023
#Metadata exploration

library(ggplot2)
library(gridExtra)

metadata <- read_excel("Documents/ETH/HS23/data/metadata.xlsx")

colnames(metadata) <- gsub(" ", "_", colnames(metadata))

# Select baseline samples
vc <- baseline$patient

metadata_baseline <- metadata[metadata$Sample_ID %in% vc,]


ggplot(metadata_baseline, aes(x=Sex)) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()


ggplot(metadata_baseline, aes(x=Histology)) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))

ggplot(metadata_baseline, aes(x=Molecular_subtype)) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))

ggplot(metadata_baseline, aes(x=Briding_therapy)) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))

ggplot(metadata_baseline, aes(x=Number_of_prior_therapies)) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()

ggplot(metadata_baseline, aes(x=Prior_autologous_transplant)) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()

ggplot(metadata_baseline, aes(x=Ann_Arbor_Stage)) + 
  geom_bar(stat = "count", fill = "deepskyblue4")+
  theme_bw()

ggplot(metadata_baseline, aes(x=age_at_diagnosis)) + 
  geom_density()+
  theme_bw()
