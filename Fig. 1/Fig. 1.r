#R_code (v4.2.0) for Figure 1
#load all package
library(geojsonsf)
library(sf)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(ggspatial)
library(cowplot)
library(tidyverse)
library(scales)
library(ggforce)
library(ggsci)
library(dplyr)
library(readxl)
library(reshape2)
library(ggpubr)
library(ggridges)

#work dir
setwd("D:/Process_work/huge_virus/Rdata/R_for_github/Fig. 1/")

#Fig. 1a
#Qingdao map 
QD <- sf::read_sf("qingdao.json")
mydata <- read.delim('Fig. 1a_QD.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
ggplot()+
  geom_sf(data=QD,aes(fill=NULL),size=0.1,color="#363636")+
  labs(x=NULL,y=NULL)+
  geom_sf(data=QD,fill="#F0F0F0",size=0.1,color="#363636")+
  xlim(120.5,121)+ylim(36.1,36.6)+
  theme_bw()+
  theme(aspect.ratio = 1.3,panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black",linetype=5,size=0.4),
        legend.title = element_blank(),
        axis.ticks = element_line(size = 0.25),
        axis.ticks.length = unit(0, "cm"),
        axis.text.x = element_text(size=7,colour = 'black'),
        axis.text.y=element_text(size=7,colour = 'black'))+
  geom_point(data=mydata,aes(x=Latitude,y=Longitude),
             fill="white",alpha=0.8,shape=22,size=4)

#China map 
#read map data
API_pre = "http://xzqh.mca.gov.cn/data/"
China = st_read(dsn = paste0(API_pre, "quanguo.json"), stringsAsFactors=FALSE) 
st_crs(China) = 4326
China_line = st_read(dsn = paste0(API_pre, "quanguo_Line.geojson"), stringsAsFactors=FALSE) 
st_crs(China_line) = 4326
gjx <- China_line[China_line$QUHUADAIMA == "guojiexian",]
province <- read.delim("province.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#read map colour
colour <- read.csv("colour.csv")
#read minimap data
nine_lines = read_sf('nanhai.geojson') 
#define colour
colour$new_colour <- rep(0,nrow(colour))
colour$new_colour[which(colour$shengfen=="Jiangsu")] <- 1
colour$new_colour[which(colour$shengfen=="Guangdong")] <- 1
colour$new_colour[which(colour$shengfen=="Guangxi")] <- 1
colour$new_colour[which(colour$shengfen=="Liaoning")] <- 1
colour$new_colour[which(colour$shengfen=="Zhejiang")] <- 1
colour$new_colour[which(colour$shengfen=="Fujian")] <- 1
colour$new_colour[which(colour$shengfen=="Shandong")] <- 1
colour$new_colour[which(colour$shengfen=="Hainan")] <- 1
colour$QUHUADAIMA <- as.character(colour$QUHUADAIMA)
China <- dplyr::left_join(China,colour,by= "QUHUADAIMA")

data <- read.delim('Fig. 1a.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#gglot China map
ggplot()+
  geom_sf(data =China,aes(fill = factor(new_colour)),size=1)+
  geom_sf(data = gjx)+
  theme_test()+ 
  theme(aspect.ratio = 1,axis.text = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color="white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,color="white",linetype=1,size=0),
        legend.text = element_text(size=7),plot.margin=unit(c(0,0,0,0),"mm"),
        legend.title = element_text(size=7),)+
  scale_fill_manual(values = c("white","#E0E0E0"))+
  geom_point(data=data,aes(x=Longitude,y=Latitude,size=Number),
             fill="#6495ED",alpha=0.8,shape=21)+
  scale_size_continuous(limits = c(0, 15), breaks = seq(0, 15, 5))

#time-series map
data<- read.delim('Fig. 1a_time.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
data$Month <- factor(data$Month,levels=c("Aug","Oct","Dec","Feb","April","June"))
custom_colors <- colorRampPalette(c("#F5FAF4", "#C9E6C4", "#8CCD82", "#5AB56A"))(90)
ggplot(data) +
  geom_density_ridges(aes(x = Temp, y = Month,fill=Number),alpha=1,scale=0.85,linewidth = 0.4)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_gradientn(colors = custom_colors)+
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
        strip.text=element_text(size=15))

#Fig. 1c
#genome size
data<-read_excel("Fig. 1c.xlsx")
data$quality <- factor(data$quality,level=c('High','Medium','Low'))
ggplot(data, aes(x=quality, y=length)) +  
  stat_boxplot(geom = 'errorbar',width=0.2,linetype=1,size=0.3)+
  geom_boxplot(width=0.45,outlier.color="white",notch=F,size = 0.3,fill="lightgray")+
  geom_jitter(color="black",width=0.2,alpha=0.2,size=1.1)+
  theme(axis.text.x=element_text(angle=0,hjust=0.2, vjust=0.5)) +
  labs(x="",y="Genome size (kb)")+
  theme_classic()+
  theme(axis.title.x = element_text(size=16.5),axis.title.y = element_text(size=16.5),
        axis.text.x = element_text(hjust =0.5,size=14,colour = 'black'),
        axis.text.y=element_text(size=14,colour = 'black'),
        axis.line = element_line(colour = "black", size = 0.5))+
  ylim(100,650)+
  theme(legend.position = "none")

#hallmark
ggplot(data, aes(x=quality, y=hallmark)) +  
  stat_boxplot(geom = 'errorbar',width=0.2,linetype=1,size=0.3)+
  geom_boxplot(width=0.45,outlier.color="white",notch=F,size = 0.3,fill="lightgray")+
  geom_jitter(color="black",width=0.2,alpha=0.2,size=1.1)+
  theme(axis.text.x=element_text(angle=0,hjust=0.2, vjust=0.5)) +
  labs(x="",y="# of viral hallmarks")+
  theme_classic()+
  theme(axis.title.x = element_text(size=16.5),axis.title.y = element_text(size=16.5),
        axis.text.x = element_text(hjust =0.5,size=14,colour = 'black'),
        axis.text.y=element_text(size=14,colour = 'black'),
        axis.line = element_line(colour = "black", size = 0.5))+
  ylim(0,230)+
  theme(legend.position = "none")

#Contamination
ggplot(data, aes(x=quality, y=contamination)) +  
  stat_boxplot(geom = 'errorbar',width=0.2,linetype=1,size=0.3)+
  geom_boxplot(width=0.45,outlier.color="white",notch=F,size = 0.3,fill="lightgray")+
  geom_jitter(color="black",width=0.2,alpha=0.2,size=1.1)+
  theme(axis.text.x=element_text(angle=0,hjust=0.2, vjust=0.5)) +
  labs(x="",y="Host contamination (%)")+
  theme_classic()+
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text.x = element_text(hjust =0.5,size=13.5,colour = 'black'),
        axis.text.y=element_text(size=13.5,colour = 'black'),
        axis.line = element_line(colour = "black", size = 0.5))+
  theme(legend.position = "none")