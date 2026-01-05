#R_code (v4.2.0) for Figure 3
#load all package
library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(tidyverse)
library(scales)
library(ggforce)
library(ggsci)
library(ape)
library(picante)

#work dir
setwd("Fig. 3/")

#number of host
data<-read_excel("Fig. 3.xlsx")
data$Virus_tax <- factor(data$Virus_tax,levels=c("Biggiephage","Jabbarphage","Judaphage","Whopperphage",
                                                 "Kabirphage","Mahaphage","Unclassified"))
data$Host_tax <- factor(data$Host_tax,levels=c('Planctomycetota', 'Nitrospirota','Chloroflexota',
                                               'Thermotogota','Patescibacteria','Actinobacteriota',
                                               'Cyanobacteria','Firmicutes','Bacteroidota','Proteobacteria'))
ggplot(data, aes(x=Host_tax, y=Number))+
  geom_bar(stat="identity", position="stack", aes(fill=Virus_tax),alpha=0.8,width=0.8)+
  scale_fill_manual(breaks=c("Biggiephage","Jabbarphage","Judaphage","Whopperphage","Kabirphage","Mahaphage"),
  values = c("#78be59","#f5c533","#CA8EFF", "#f6b26b","#9cceff","#2894ff"))+
  labs(x="",y="# of viral genomes")+
  theme(legend.background=element_rect(color="white"),)+
  theme(legend.key = element_blank())+
  theme(axis.line = element_line(color=))+
  theme_test()+coord_flip()+theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
                     axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
                     axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
                     axis.text.y=element_text(size=12,colour = 'black'),
                     panel.border = element_rect(size=0.8),
                     legend.text = element_text(size=12),
                     legend.title = element_text(size=12))+ylim(0,30)+theme(legend.position = "none")

#fig pie
sales <- c(79.5,20.5)
group <- c("Reference", "Current")
share <- sales/sum(sales)*100
data <- data.frame(sales,share,group)
pie <- ggplot()+
  geom_arc_bar(data=data,aes(x0 = 0, y0 = 0, r0 = 0, r = 1,amount=share,fill=group),stat="pie",alpha=0.9)+
  coord_fixed()+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("Reference","Current"),values=c("#7AC5CD","#f6b26b"))

#fig PD
tree<-read.tree(file="Terl210_iqtree.contree")

#calculate pd
edges <- tree$edge
edge_lengths <- tree$edge.length
tip_labels <- tree$tip.label

get_node_label <- function(node) {
  if (node <= length(tip_labels)) {
    return(tip_labels[node])
  } else {
    return(paste("Node", node, sep = "_"))
  }
}

edge_info <- data.frame(
  Start = sapply(edges[, 1], get_node_label),
  End = sapply(edges[, 2], get_node_label),
  Length = edge_lengths
)

#percentage of PD
sales <- c(38.3,36.6,25.1)
group <- c("Reference","Shared","Current")
share <- sales/sum(sales)*100
data <- data.frame(
  sales,share,group)
ggplot()+
  geom_arc_bar(data=data,aes(x0 = 0, y0 = 0, r0 = 0, r = 1,amount=share,fill=group),stat="pie",alpha=0.9)+
  coord_fixed()+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("Reference","Shared","Current"),values=c("#7AC5CD","#E0E0E0","#f6b26b"))

#PD for each lineage
data <- read_excel("Fig. 3.xlsx")
data$Host_tax <- factor(data$Host_tax,levels=c('Planctomycetota', 'Nitrospirota','Chloroflexota',
                                               'Thermotogota','Patescibacteria','Actinobacteriota',
                                               'Cyanobacteria','Firmicutes','Bacteroidota','Proteobacteria'))
ggplot(data, aes(x=Host_tax, y=PD))+
  geom_bar(stat="identity", fill="#d0d0d0",alpha=0.8,width=0.8)+
  labs(x="",y="Phylogenetic diversity (PD)")+
  theme_test()+theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.x = element_text(size=13),
                     axis.title.y = element_text(size=13),
                     axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
                     axis.text.y = element_text(size=12,colour = 'black'),
                     panel.border = element_rect(size=0.8),
                     legend.text = element_text(size=12),
                     legend.title = element_text(size=12),
                     legend.position = "none",
                     axis.ticks.y.left = element_blank())+
  coord_flip(ylim = c(0, 15))

