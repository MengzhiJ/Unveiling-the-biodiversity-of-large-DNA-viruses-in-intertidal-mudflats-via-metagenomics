#R_code (v4.2.0) for Figure 4
#loading package
library(Hmisc)
library(psych)
library(data.table)
library(dplyr)
library(vegan)
library(reshape2)
library(tidyverse)
library(scales)
library(ggforce)
library(ggsci)
library(ggpmisc)
library(ggraph)
library(tidygraph)
library(igraph)
library(ggplot2)
library(readxl)
library(ggpubr)

#work dir
setwd('Fig. 5/')

#Mantel_test (control geographic distance)
#extract_func
extract_columns <- function(df, pattern) {
  cols_to_extract <- grep(pattern, names(df), value = TRUE)
  if (length(cols_to_extract) == 0) {
    return(data.frame())
  }
  df_extracted <- df %>%
    select(all_of(cols_to_extract))
  return(df_extracted)
}

#correlation between euk_NCLDV lineages
votu<-t(read.delim('huge_mantel.txt',row.names = 1)) #species as row, sample as col
votu_bray<-vegdist(votu,method = "bray")

geo<-read.csv('distance.csv',row.names = 1)
geo<-vegdist(geo,method = 'euclidean')

euk<-read.delim('euk_mantel.txt',row.names = 1)#species as col, sample as row

#Imi_42_exist_species_remove_singleton
lineages <- c("Amoebozoa","Annelida","Apicomplexa", "Arthropoda","Ascomycota",
              "Bacillariophyta","Basidiomycota","Bigyra","Bolidophyceae",
              "Chlorophyta","Choanoflagellata","Chordata","Chytridiomycota",
              "Ciliophora","Cnidaria","Cryptophyceae",
              "Dictyochophyceae","Dinophyceae", "Discoba","Echinodermata",
              "Eustigmatophyceae","Haptophyta", "Metamonada","Microsporidia","Mollusca",
              "Mucoromycota","Nematoda","Oomycota", "Streptophyta",
              "Pelagophyceae","Perkinsozoa","Phaeophyceae","Pinguiophyceae",
              "Placozoa","Platyhelminthes","Porifera",
              "Raphidophyceae","Rhizaria","Rhodophyta","Rotifera","Tardigrada","Zoopagomycota"
)

euk_bray <- lapply(lineages, function(taxa) {
  dat <- extract_columns(euk, taxa)
  if (ncol(dat) == 0) return(NULL)
  vegdist(dat, method = "bray")
})

#partial_control_geo
r<-c()
p<-c()
for (i in 1:42){
  r[i]<-mantel.partial(votu_bray,euk_bray[[i]],geo,permutations = 999,method="pearson",na.rm=T)$statistic
  p[i]<-mantel.partial(votu_bray,euk_bray[[i]],geo,permutations = 999,method="pearson",na.rm=T)$signif
}

p<-as.data.frame(p)
r<-as.data.frame(r)
rownames(p)<-lineages
rownames(r)<-lineages

#Pim_42_exist_species_remove_singletons
lineages <- c("Amoebozoa","Annelida","Apicomplexa", "Arthropoda","Ascomycota",
              "Bacillariophyta","Basidiomycota","Bigyra","Bolidophyceae",
              "Chlorophyta","Choanoflagellata","Chordata","Chytridiomycota",
              "Ciliophora","Cnidaria","Cryptophyceae",
              "Dictyochophyceae","Dinophyceae", "Discoba","Echinodermata",
              "Eustigmatophyceae","Haptophyta","Oomycota",
              "Metamonada","Microsporidia","Mollusca","Mucoromycota","Nematoda",
              "Pelagophyceae","Perkinsozoa","Phaeophyceae","Pinguiophyceae",
              "Placozoa","Platyhelminthes","Porifera",
              "Raphidophyceae","Rhizaria","Rhodophyta","Rotifera",
              "Streptophyta","Tardigrada","Zoopagomycota"
)

euk_bray <- lapply(lineages, function(taxa) {
  dat <- extract_columns(euk, taxa)
  if (ncol(dat) == 0) return(NULL)
  vegdist(dat, method = "bray")
})

#partial_control_geo
r<-c()
p<-c()
for (i in 1:42){
  r[i]<-mantel.partial(votu_bray,euk_bray[[i]],geo,permutations = 999,method="pearson",na.rm=T)$statistic
  p[i]<-mantel.partial(votu_bray,euk_bray[[i]],geo,permutations = 999,method="pearson",na.rm=T)$signif
}

p<-as.data.frame(p)
r<-as.data.frame(r)
rownames(p)<-lineages
rownames(r)<-lineages

#Pand_41_exist_species_remove_singletons
lineages <- c("Amoebozoa","Annelida","Apicomplexa", "Arthropoda","Ascomycota",
              "Bacillariophyta","Basidiomycota","Bigyra","Bolidophyceae",
              "Chlorophyta","Choanoflagellata","Chordata","Chytridiomycota","Ciliophora",
              "Cnidaria","Cryptophyceae","Dictyochophyceae","Dinophyceae", 
              "Discoba","Echinodermata","Eustigmatophyceae","Haptophyta",
              "Metamonada","Microsporidia","Mollusca","Mucoromycota","Nematoda","Oomycota",
              "Pelagophyceae","Phaeophyceae","Pinguiophyceae","Placozoa",
              "Platyhelminthes","Porifera",
              "Raphidophyceae","Rhizaria","Rhodophyta","Rotifera",
              "Streptophyta","Tardigrada","Zoopagomycota"
)

euk_bray <- lapply(lineages, function(taxa) {
  dat <- extract_columns(euk, taxa)
  if (ncol(dat) == 0) return(NULL)
  vegdist(dat, method = "bray")
})

#partial_control_geo
r<-c()
p<-c()
for (i in 1:41){
  r[i]<-mantel.partial(votu_bray,euk_bray[[i]],geo,permutations = 999,method="pearson",na.rm=T)$statistic
  p[i]<-mantel.partial(votu_bray,euk_bray[[i]],geo,permutations = 999,method="pearson",na.rm=T)$signif
}

p<-as.data.frame(p)
r<-as.data.frame(r)
rownames(p)<-lineages
rownames(r)<-lineages

#alg_39_exist_species_remove_singletons
lineages <- c("Amoebozoa","Annelida","Apicomplexa", "Arthropoda","Ascomycota",
              "Bacillariophyta","Basidiomycota","Bigyra","Bolidophyceae",
              "Chlorophyta","Choanoflagellata","Chordata","Chytridiomycota","Ciliophora",
              "Cnidaria","Cryptophyceae","Dictyochophyceae","Dinophyceae", 
              "Discoba","Echinodermata","Eustigmatophyceae","Haptophyta",
              "Metamonada","Microsporidia","Mollusca","Mucoromycota","Nematoda","Oomycota",
              "Pelagophyceae","Phaeophyceae","Pinguiophyceae",
              "Platyhelminthes","Porifera",
              "Rhizaria","Rhodophyta","Rotifera",
              "Streptophyta","Tardigrada","Zoopagomycota"
)

euk_bray <- lapply(lineages, function(taxa) {
  dat <- extract_columns(euk, taxa)
  if (ncol(dat) == 0) return(NULL)
  vegdist(dat, method = "bray")
})

#partial_control_geo
r<-c()
p<-c()
for (i in 1:39){
  r[i]<-mantel.partial(votu_bray,euk_bray[[i]],geo,permutations = 999,method="spearman",na.rm=T)$statistic
  p[i]<-mantel.partial(votu_bray,euk_bray[[i]],geo,permutations = 999,method="spearman",na.rm=T)$signif
}

p<-as.data.frame(p)
r<-as.data.frame(r)
rownames(p)<-lineages
rownames(r)<-lineages

#spearman correlation
#correlation between imi_NCLDV-euk pairs (occurrence >= 10% samples)
combine<-read.delim('combine_Imitervirales.txt',row.names = 1)
combine<-as.data.frame(t(combine))
cor<-rcorr(as.matrix(combine),type='spearman')
r<-cor$r
p<-cor$P
r<-r[1:25,26:282]
p<-p[1:25,26:282]
p_adjust<- matrix(p.adjust(as.vector(p), method="bonferroni"),
                  nrow=nrow(p), ncol=ncol(p))
data1 = melt(r, value.name="r")
data2 = melt(p_adjust, value.name="p_adjust")
data = cbind(data1, p_adjust = data2$p_adjust)
filtered_rows <- data[abs(data$r)>0.6&data$p_adjust< 0.05, ]
clean_data<-na.omit(filtered_rows)

#correlation between Pim_NCLDV-euk pairs (occurrence >= 10% samples)
combine<-read.delim('combine_Pimascovirales.txt',row.names = 1)
combine<-as.data.frame(t(combine))
cor<-rcorr(as.matrix(combine),type='spearman')
r<-cor$r
p<-cor$P
r<-r[1:9,10:131]
p<-p[1:9,10:131]
p_adjust<- matrix(p.adjust(as.vector(p), method="bonferroni"),
                  nrow=nrow(p), ncol=ncol(p))
data1 = melt(r, value.name="r")
data2 = melt(p_adjust, value.name="p_adjust")
data = cbind(data1, p_adjust = data2$p_adjust)
filtered_rows <- data[abs(data$r)>0.6&data$p_adjust< 0.05, ]
clean_data<-na.omit(filtered_rows)

#correlation between Pan_NCLDV-euk pairs (occurrence >= 10% samples)
combine<-read.delim('combine_Pandoravirales.txt',row.names = 1)
combine<-as.data.frame(t(combine))
cor<-rcorr(as.matrix(combine),type='spearman')
r<-cor$r
p<-cor$P
r<-r[1:11,12:215]
p<-p[1:11,12:215]
p_adjust<- matrix(p.adjust(as.vector(p), method="bonferroni"),
                  nrow=nrow(p), ncol=ncol(p))
data1 = melt(r, value.name="r")
data2 = melt(p_adjust, value.name="p_adjust")
data = cbind(data1, p_adjust = data2$p_adjust)
filtered_rows <- data[abs(data$r)>0.6&data$p_adjust< 0.05, ]
clean_data<-na.omit(filtered_rows)

#plot
setwd('Fig. 5/Fig/')
#Fig.5a
#bar
data<-read.delim('Fig. 5a.txt',row.names = 1)
ggplot(data, aes(x=correlation_r,y=..density..))+  
  geom_histogram(binwidth = 0.02,alpha=0.3,colour="black",size=0.2,fill="gray")+
  geom_density(alpha=.45,fill="gray",color="gray",lwd=0.22)+
  theme_test()+
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(hjust =0.5,size=7,colour = 'black'),
        axis.text.y=element_text(size=7,colour = 'black'),
        axis.ticks = element_line(size = 0.25),
        axis.ticks.length = unit(0.05, "cm"),
        panel.border = element_rect(size=0.45))+
  labs(x="Mantel correlation R",y="Density")+ylim(0,6)+ 
  scale_x_continuous(breaks = seq(-1, 1, by = 0.2))+
  geom_vline(xintercept=c(0.53), linetype=2,color="#66B3FF",size=0.4)

#pie
sales <- c(41.9,58.1)
names<-c("Known","Unknown")
share<-sales/sum(sales)*100
data <- data.frame(
  sales,share,names)
ggplot()+
  geom_arc_bar(data=data,aes(x0 = 0, y0 = 0, r0 = 0, r = 1,amount=sales,explode=c(0,0),fill=names),
               stat="pie",alpha=0.6)+
  coord_fixed()+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("Known","Unknown"),values=c("#66B3FF","lightgray"))

#Fig.5b
data<-read.delim('Fig. 5b.txt',row.names = 1)
ggplot(data, mapping = aes(x=Euk,y=GV))+
  geom_point(fill="lightgray",size = 1.2, alpha = 0.85,shape=21,stroke=0.22)+
  labs(x="Abundance of eukaryotes",y="Abundance of GVs")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(color="#66B3FF",method = 'lm',se=T,size=0.4,fullrange=T,fill="lightgray",alpha=0.25)+
  theme_test()+ theme(legend.position ="none")+
  theme( axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         axis.text.x = element_text(hjust =0.5,size=7,colour = 'black'),
         axis.text.y=element_text(size=7,colour = 'black'),
         axis.ticks = element_line(size = 0.25),
         axis.ticks.length = unit(0.05, "cm"),
         panel.border = element_rect(size=0.45))+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),size=2.3,
           method="pearson")

#Fig.5c
node <- read.csv("node.csv")
edge <- read.csv("edge.csv")
network <- tbl_graph(nodes = node, edges = edge, directed = FALSE)
node_colors <- c(
  "Virus" = "#66B3FF",   
  "Algae" = "#4CAF50",   
  "Animal" = "#FF9224",  
  "Fungi" = "#d3a4ff",   
  "Protozoa" = "#FFD306")

edge_colors <- c(
  "Known" = "#84C1FF",   
  "Unknown" = "lightgray")

ggraph(network, layout = "circle") + 
  geom_edge_arc(aes(edge_width = number,edge_color=lineage), curvature = 0.1,
                alpha = 0.4) +
  scale_edge_width(range = c(0.7, 2.2)) +
  geom_node_point(aes(fill=type,shape=shape), size = 7) + 
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = node_colors) + 
  scale_edge_color_manual(values = edge_colors) + 
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void()