#R_code (v4.2.0) for Figure 4
#loading package
library(gplots)  
library(pheatmap) 
library(readxl) 
library(vegan)
library(ggpubr)
library(ggpmisc)
library(ggplot2)
library(SoDA)
library(splines)
library(ape)
library(ggridges)
library(betapart)
library(NST)
#Fig. 4a
#PCoA_national
setwd("E:/LDVs/Rdata/R_for_github/Fig. 4/national_PCoA/")
otu <- read.delim('NCLDV_OTU_national_bray.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu1 <- data.frame(t(otu))
otu2 <- log(otu1+1)
group <- read.delim('group_national.txt', sep = '\t', stringsAsFactors = FALSE)
dist <- vegdist(otu2, method = "bray")
pcoa_result <- pcoa(dist)
sample_site <- data.frame(
  PCoA1 = pcoa_result$vectors[, 1],
  PCoA2 = pcoa_result$vectors[, 2],
  names = rownames(pcoa_result$vectors))

sample_site <- cbind(sample_site, group[, -1])
variance_explained <- round(pcoa_result$values$Relative_eig[1:2] * 100, 2)
PERMANOVA <- adonis2(dist ~ Sample, data = group, permutations = 999)

otu <- read.csv('PCOA.csv')
ggplot(data = otu, mapping = aes(PCoA1, PCoA2)) + 
  geom_point(aes(fill = Site,shape=Tax), size = 2, alpha = 0.9) + 
  scale_shape_manual(breaks = c('Nucleocytoviricota','Uroviricota'),values=c(21,24))+
  scale_fill_manual(breaks=c('DD','DY','QD','LYG','YC','NB','WZ','XM','ST','ZH','BH','SY'),
                    values=c("#6A5ACD","#6495ED","#1E90FF","#7CCD7C","#B4EEB4","#EE9A00","#FFAF60",
                             "#DDA0DD", "#EAC100","#FFD700","#EE7942","#20B2AA"))+  
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
        axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+facet_wrap(~Tax)

#PCoA_temporal
setwd("temporal_PCoA/")
otu <- read.delim('Phage_OTU_time.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu1 <- data.frame(t(otu))
otu2 <- log(otu1+1)
group <- read.delim('group_time.txt', sep = '\t', stringsAsFactors = FALSE)
dist <- vegdist(otu2, method = "bray")
pcoa_result <- pcoa(dist)
sample_site <- data.frame(
  PCoA1 = pcoa_result$vectors[, 1],
  PCoA2 = pcoa_result$vectors[, 2],
  names = rownames(pcoa_result$vectors))

sample_site <- cbind(sample_site, group[, -1])
variance_explained <- round(pcoa_result$values$Relative_eig[1:2] * 100, 2)
PERMANOVA<-adonis2(dist ~ Sample, data = group, permutations = 999)

otu<-read.csv('PCOA.csv')
ggplot(data = otu, mapping = aes(PCoA1, PCoA2)) + 
  geom_point(aes(fill = Time,shape=Tax), size = 2, alpha = 0.9) + 
  scale_shape_manual(breaks = c('Nucleocytoviricota','Uroviricota'),values=c(21,24))+
  scale_fill_manual(breaks=c("Aug","Oct","Dec","Feb","April","June"),
                    values=c("#D4EDFF", "#8CC0FF", "#5AA9FF", "#2894FF", "#0068C6", "#003D79"))+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
        axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+facet_wrap(~Tax)+xlim(-0.5,0.5)+ylim(-0.4,0.4)

#PCOA_vertical
setwd("vertical_PCoA/")
otu <- read.delim('phage_OTU_depth.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu1 <- data.frame(t(otu))
otu2 <- log(otu1+1)
group <- read.delim('group_depth.txt', sep = '\t', stringsAsFactors = FALSE)
dist <- vegdist(otu2, method = "bray")
pcoa_result <- pcoa(dist)
sample_site <- data.frame(
  PCoA1 = pcoa_result$vectors[, 1],
  PCoA2 = pcoa_result$vectors[, 2],
  names = rownames(pcoa_result$vectors))

sample_site <- cbind(sample_site, group[, -1])
variance_explained <- round(pcoa_result$values$Relative_eig[1:2] * 100, 2)
PERMANOVA<-adonis2(dist ~ Sample, data = group, permutations = 999)

otu<-read.csv('PCOA.csv')
ggplot(data = otu, mapping = aes(PCoA1, PCoA2)) + 
  geom_point(aes(fill = Depth,shape=Tax), size = 2, alpha = 0.9) + 
  scale_shape_manual(breaks = c('Nucleocytoviricota','Uroviricota'),values=c(21,24))+
  scale_fill_manual(breaks=c("Twenty","Forty","Sixty","Eighty","Hundred"),
                    values=c("#E8F7F9","#C1EBF0","#92D3D6","#55B5B7","#3F8D93"))+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
        axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+facet_wrap(~Tax)+ylim(-0.25,0.3)

#Fig. 4b DDR
#distance_calculate_national
setwd("national_DDR/")
otu<-read.delim('NCLDV_OTU_national_bray.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
loca<-read_excel('DDR.xlsx')

distv<-vegdist(t(otu),method = 'bray')
distv_num<-1-as.numeric(distv)
distv_num<-log10(distv_num)
distv_num<-as.numeric(distv_num)
loca_geo<-geoXY(loca$Latitude,loca$Longitude)
rownames(loca_geo)<-rownames(loca)
dist_loca<-vegdist(loca_geo,method = 'euclidean')
dist_loca_log<-log10(dist_loca)
dist_loca_num<-as.numeric(dist_loca_log)
dist<-data.frame(distv_num,dist_loca_num)
#dist_df <- melt(as.matrix(distv))
#dist_df <- melt(as.matrix(dist_loca))

#Plot_DDR_national
data<-read.csv('bray_national.csv')
data$Tax=factor(data$Tax, levels=c("Nucleocytoviricota","Uroviricota"))
data<-data[order(data$Tax, decreasing = TRUE), ]
ggplot(data,aes(x=dist_loca_num,y=distv_num))+
  geom_point(aes(color=Tax),size = 1, alpha = 0.5,shape=16)+
  theme(axis.line = element_line(color="black"))+
  labs(x="",y="")+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#6ED0CF"))+ 
  scale_color_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                     values=c("#6BB5FF","#6ED0CF"))+
  geom_smooth(aes(color=Tax,fill=Tax),method = 'lm',se=T,size=0.8,fullrange=F,alpha=0.2)+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=15),
        legend.position = "none")+ylim(-2.5,0.3)
  stat_poly_eq(aes(color=Tax,label = paste(..eq.label..,..p.value.., sep = "~`,`~")),size=5,
               formula = y ~ x,parse = TRUE)

#DDR_local
setwd("local_DDR/")
#distance_calculate_local
otu<-read.delim('Phage_local_OTU.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
distv<-vegdist(t(otu),method = 'bray')
distv_num<-1-as.numeric(distv)
distv_num<-log10(distv_num)
distv_num<-as.numeric(distv_num)
loca_geo<-read.delim('DDR.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
dist_loca<-vegdist(loca_geo,method = 'euclidean')
dist_loca_log<-log10(dist_loca)
dist_loca_num<-as.numeric(dist_loca_log)
dist<-data.frame(distv_num,dist_loca_num)

#Plot_DDR_local
data<-read.csv('bray_local.csv')
data$Tax=factor(data$Tax, levels=c("Nucleocytoviricota","Uroviricota"))
data<-data[order(data$Tax, decreasing = TRUE), ]
ggplot(data,aes(x=dist_loca_num,y=distv_num))+
  geom_point(aes(color=Tax),size = 1, alpha = 0.5,shape=16)+
  theme(axis.line = element_line(color="black"))+
  labs(x="",y="")+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#6ED0CF"))+ 
  scale_color_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                     values=c("#6BB5FF","#6ED0CF"))+
  geom_smooth(aes(color=Tax,fill=Tax),method = 'lm',se=T,size=0.8,fullrange=F,alpha=0.2)+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=15),
        legend.position = "none")+
  ylim(-0.5,0.07)+
  stat_poly_eq(aes(color=Tax,label = paste(..eq.label..,..p.value.., sep = "~`,`~")),size=5,
               formula = y ~ x,parse = TRUE)

#Fig. 4c TDR
setwd("temporal_TDR/")
#distance_calculate_time
otu<-read.delim('NCLDV_OTU_time.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
distv<-vegdist(t(otu),method = 'bray')
distv_num<-as.numeric(1-distv)
Time<-read.delim('Time.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
dist_loca<-vegdist(Time,method = 'euclidean')
dist_loca_num<-as.numeric(dist_loca)
dist<-data.frame(distv_num,dist_loca_num)

data<-read.csv('bray_time.csv')
data$Tax=factor(data$Tax, levels=c("Uroviricota","Nucleocytoviricota"))
data<-data[order(data$Tax, decreasing = TRUE), ]
ggplot(data,aes(x=dist_loca_num,y=distv_num))+
  geom_point(aes(color=Tax),size = 1, alpha = 0.3,shape=16)+
  theme(axis.line = element_line(color="black"))+
  labs(x="",y="")+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#6ED0CF"))+ 
  scale_color_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                     values=c("#6BB5FF","#6ED0CF"))+
  geom_smooth(aes(color=Tax,fill=Tax),method = 'lm',se=T,size=0.8,fullrange=F,alpha=0.2)+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=15),
        legend.position = "none")+scale_x_continuous(limits = c(0, 10),breaks = seq(0, 10, 2))+
  scale_y_continuous(limits = c(-0.05, 1.1),breaks = seq(0, 1, 0.2))+
  stat_poly_eq(aes(color=Tax,label = paste(..eq.label..,..p.value..,sep = "~`,`~")),size=5,
               formula = y ~ x,parse = TRUE)

#vertical_DDR
setwd("vertical_DDR/")
#Fig. 4c DDR_depth
#distance_calculate_depth
otu<-read.delim('phage_OTU_depth.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
distv<-vegdist(t(otu),method = 'bray')
distv_num<-as.numeric(1-distv)
Depth<-read.delim('depth.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
dist_loca<-vegdist(Depth,method = 'euclidean')
dist_loca_num<-as.numeric(dist_loca)
dist<-data.frame(distv_num,dist_loca_num)
  
data<-read.csv('bray_depth.csv')
data$Tax=factor(data$Tax, levels=c("Nucleocytoviricota","Uroviricota"))
data<-data[order(data$Tax, decreasing = TRUE), ]
ggplot(data,aes(x=dist_loca_num,y=distv_num))+
    geom_point(aes(color=Tax),size = 1, alpha = 0.2,shape=16)+
    theme(axis.line = element_line(color="black"))+
    labs(x="",y="")+
    scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                      values=c("#6BB5FF","#6ED0CF"))+ 
    scale_color_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                       values=c("#6BB5FF","#6ED0CF"))+
    geom_smooth(aes(color=Tax,fill=Tax),method = 'lm',se=T,size=0.8,fullrange=F,alpha=0.2)+
    theme_test()+ theme(legend.position ="none")+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
          axis.text.y=element_text(size=12,colour = 'black'),
          panel.border = element_rect(size=0.8),
          axis.ticks = element_line(size = 0.4),
          axis.ticks.length = unit(0.08, "cm"),
          strip.text=element_text(size=15),
          legend.position = "none")+scale_x_continuous(limits = c(0, 0.8),breaks = seq(0, 0.8, 0.2))+
  scale_y_continuous(limits = c(0.15, 1.1),breaks = seq(0, 1, 0.2))+
    stat_poly_eq(aes(color=Tax,label = paste(..eq.label..,..p.value.., sep = "~`,`~")),size=5,
                 formula = y ~ x,parse = TRUE)
       
#Fig. 4d
#Pi value_time
setwd("Nucleotide_variations/")
data<-read_excel("Fig. 4d.xlsx")
data$Time <- factor(data$Time,levels=c("Aug","Oct","Dec","Feb","April","June"))
ggplot(data, aes(x = Pi, y = Time,fill=Tax)) +
  geom_density_ridges(alpha=0.8,scale=0.6,linewidth = 0.3)+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#6ED0CF"))+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=15),
        legend.position = "none")+xlim(-0.0004,0.0055)

data<-read_excel("Fig. 4d.xlsx")
data$Time <- factor(data$Time,levels=c("Aug","Oct","Dec","Feb","April","June"))
ggplot(data,aes(x=Time,y=Pi,fill=Tax,group=Tax))+
  stat_summary(alpha=0.8,fun = mean,geom="bar",color="black",width=0.7,size=0.3,position = position_dodge(0.7))+
  stat_summary(fun.data = mean_se,geom="errorbar",width=.2,size=0.3,position = position_dodge(0.7))+
  geom_jitter(aes(fill=Tax),shape=21,size=1,alpha=0.9,stroke = 0.2,position=position_jitterdodge(jitter.width = 0.25,dodge.width = 0.8))+
  labs(x="",y="")+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#6ED0CF"))+
  scale_color_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#6ED0CF"))+
  theme(axis.line = element_line(color="black"))+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_blank(),
        panel.border = element_rect(size=0.9),
        axis.ticks.y = element_blank(),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+coord_flip()+scale_y_continuous(limits=c(0,0.0012),breaks = c(0,0.0006,0.0012))

#Pi value_depth
data<-read_excel("Fig. 4d.xlsx")
data$Depth <- factor(data$Depth,levels=c("20","40","60","80","100"))
ggplot(data, aes(x = Pi, y = Depth,fill=Tax)) +
  geom_density_ridges(alpha=0.8,scale=0.5,linewidth = 0.3)+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#6ED0CF"))+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=15),
        legend.position = "none")

data<-read_excel("Fig. 4d.xlsx")
data$Depth <- factor(data$Depth,levels=c("20","40","60","80","100"))
ggplot(data,aes(x=Depth,y=Pi,fill=Tax,group=Tax))+
  stat_summary(alpha=0.8,fun = mean,geom="bar",color="black",width=0.6,size=0.3,position = position_dodge(0.6))+
  stat_summary(fun.data = mean_se,geom="errorbar",width=.2,size=0.3,position = position_dodge(0.7))+
  geom_jitter(aes(fill=Tax),shape=21,size=1,alpha=0.9,stroke = 0.2,position=position_jitterdodge(jitter.width = 0.25,dodge.width = 0.7))+
  labs(x="",y="")+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#6ED0CF"))+
  scale_color_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                     values=c("#6BB5FF","#6ED0CF"))+
  theme(axis.line = element_line(color="black"))+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_blank(),
        panel.border = element_rect(size=0.9),
        axis.ticks.y = element_blank(),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+coord_flip()+scale_y_continuous(limits=c(0,0.0014),breaks = c(0,0.0007,0.0014))
