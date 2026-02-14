#R_code (v4.2.0) for Figure 6
#loading package
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpmisc)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(SoDA)
library(linkET)
library(Hmisc)
library(minpack.lm)
library(stats4)
library(grid)
library(ape)
library(parallel)
library(bigmemory)
library(NST)
library(permute)
library(DirichletReg)
library(iCAMP)
library(reshape2)

#Fig. 6a environment factors
#NCLDV
setwd("Fig. 6/Mantel/NCLDV/")
env<-read.csv('Environment_ncldv.csv',row.names = 1) 
OTUs<-t(read.delim('NCLDV_OTU.txt',row.names = 1))
geo<-read.csv('distance_ncldv.csv',row.names = 1) 
samp.ck=match.name(rn.list=list(env,OTUs))
samp.ck=match.name(rn.list=list(geo,OTUs))
env<-apply(env[,1:11],2,as.numeric)

OTUs<-vegdist(OTUs,method = "bray")
geo<-vegdist(geo,method = "euclidean")

#mantel test for each env with OTUs
mantel<-list()
for (i in 1:11){
  mantel[[i]]<-vegdist(env[,i],method='euclidean',upper = FALSE,na.rm = T)
}
names(mantel)<-colnames(env)[1:11]

#mantel test for env and vOTU
r<-c()
p<-c()
for (i in 1:11){
  r[i]<-mantel.partial(OTUs,mantel[[i]],geo,permutations = 999,method="pearson",na.rm=T)$statistic
  p[i]<-mantel.partial(OTUs,mantel[[i]],geo,permutations = 999,method="pearson",na.rm=T)$signif
}
p<-as.data.frame(p)
r<-as.data.frame(r)
rownames(p)<-colnames(env)
rownames(r)<-colnames(env)

#Phage
setwd("Fig. 6/Mantel/Phage/")
env<-read.csv('Environment_phage.csv',row.names = 1) 
OTUs<-t(read.delim('Phage_OTU.txt',row.names = 1))
geo<-read.csv('distance_phage.csv',row.names = 1) 
samp.ck=match.name(rn.list=list(env,OTUs))
samp.ck=match.name(rn.list=list(geo,OTUs))

env<-apply(env[,1:11],2,as.numeric)
OTUs<-vegdist(OTUs,method = "bray")
geo<-vegdist(geo,method = "euclidean")

#mantel test for each env with OTUs
mantel<-list()
for (i in 1:11){
  mantel[[i]]<-vegdist(env[,i],method='euclidean',upper = FALSE,na.rm = T)
}
names(mantel)<-colnames(env)[1:11]

#mantel test for env and vOTU
r<-c()
p<-c()
for (i in 1:11){
  r[i]<-mantel.partial(OTUs,mantel[[i]],geo,permutations = 999,method="pearson",na.rm=T)$statistic
  p[i]<-mantel.partial(OTUs,mantel[[i]],geo,permutations = 999,method="pearson",na.rm=T)$signif
}
p<-as.data.frame(p)
r<-as.data.frame(r)
rownames(p)<-colnames(env)
rownames(r)<-colnames(env)

#plot_mantel
setwd("Fig. 6/Mantel/")
mantel<-read.delim('mantel.txt')
env<-read.csv('Environment.csv',row.names = 1) 
env<-apply(env[,1:11],2,as.numeric)
qcorrplot(rcorr(env), type = "upper", diag = F, grid_size = 0.4,grid_col = "lightgray") +
  geom_square(linetype=0) +
  geom_couple(data = mantel,aes(colour = p_value, size = r_value), curvature = 0.1) +
  set_corrplot_style(colours = c("#6ED0CF", "white", "#6BB5FF")) +
  scale_size_manual(breaks = c("<0.1","0.1-0.2",">0.2"),
                    values = c(0.3, 0.6, 0.9))+
  scale_colour_manual(breaks = c(">0.05","0.001-0.05","<0.001"),
                      values=c("lightgray","#9ACD32","#87CEEB")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 2), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

#Fig. 6b,c
#NST NCLDV
setwd("Fig. 6/NST/NCLDV/")
com.file="NCLDV_OTU_all.txt" 
group.file="group_all.txt"
comm=t(read.delim(com.file, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE))
group=read.delim(group.file, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
samp.ck=NST::match.name(rn.list=list(comm=comm,group=group))
comm=samp.ck$comm
comm=comm[,colSums(comm)>0,drop=FALSE]
group=samp.ck$group

tnst=tNST(comm=comm, group=group,
          dist.method='bray', abundance.weighted=TRUE, rand=5000,
          output.rand=TRUE, nworker=28, LB=FALSE, null.model="PF",dirichlet=F,
          between.group=F, SES=T, RC=T)

#NST phage
setwd("Fig. 6/NST/phage/")
com.file="Phage_OTU_all.txt" 
group.file="group_all.txt"
comm=t(read.delim(com.file, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE))
group=read.delim(group.file, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
samp.ck=NST::match.name(rn.list=list(comm=comm,group=group))
comm=samp.ck$comm
comm=comm[,colSums(comm)>0,drop=FALSE]
group=samp.ck$group

tnst=tNST(comm=comm, group=group,
          dist.method='bray', abundance.weighted=TRUE, rand=5000,
          output.rand=TRUE, nworker=28, LB=FALSE, null.model="PF",dirichlet=F,
          between.group=F, SES=T, RC=T)

#plot Fig. 6b,c
setwd("Fig. 6/NST/")
data<-read_excel("NST.xlsx") 
data$class=factor(data$class, levels=c("Giant virus","Large phage"))
ggplot(data,aes(x=class,y=pairwise_NST))+
  stat_summary(aes(fill=class), alpha=0.8,fun = mean,geom="bar",color="NA",width=0.45)+
  stat_summary(fun.data = mean_se,geom="errorbar",width=.08,size=0.2)+
  labs(x="",y="Normalized stochastic ratio (%)")+
  scale_fill_manual(breaks=c("Giant virus","Large phage"),
                    values=c("#6BB5FF","#6ED0CF"))+ 
  theme_test()+theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13),
        axis.text.x = element_text(hjust = 0.5,size=12,angle = 0,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        legend.position = "none")+ 
  coord_cartesian(ylim = c(0,100))

data<-read_excel("Fig. 6/NST.xlsx") 
data$class=factor(data$class, levels=c("Giant virus","Large phage"))
ggplot(data,aes(x=class,y=Ratio,fill=level))+
  geom_bar(stat="identity", color="black",width=0.45,size=0.2,position="stack", aes(fill=level),alpha=0.8)+
  theme_bw()+labs(x="",y="Raup-Crick proportion (%)")+
  scale_fill_manual(breaks=c("more","mid","less"),values = c("#B3E2CD","#A8CEFF","#DDBBEA"))+
  theme_test()+theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13),
        axis.text.x = element_text(hjust = 0.5,size=12,angle = 0,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        legend.position = "none")+ 
  coord_cartesian(ylim = c(0,100))

#Fig. 6d iCAMP
#Phage
setwd("Fig. 6/iCAMP/Phage/")
otu <- read.delim('Phage_OTU.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
comm <- t(otu)
tree <- read.tree("Terl210.iqtree.contree")
tax <- read.delim('taxanomy.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
save.wd=("Fig. 6/iCAMP/Phage/")

icamp.out=icamp.big(comm=comm,tree=tree,pd.wd=save.wd,rand = 1000,ds = 0.2,
                    bin.size.limit=5,phylo.rand.scale ="across.all",
                    taxa.rand.scale = "across.all",
                    sig.index = "SES.RC",phylo.metric ="bMPD",
                    nworker = 28,ses.cut = 1.96, rc.cut = 0.95,
                    taxo.metric="bray")

#summarize each bin
icampbin=icamp.bins(icamp.detail = icamp.out, treat = NULL,
                    clas = tax, boot=F,
                    rand.time = 1000, between.group = F)

#summarize each taxon
cate=icamp.cate(icamp.bins.result = icampbin,comm = comm,cate =tax,
                  silent = FALSE,between.group = F)

#NCLDV
setwd("Fig. 6/iCAMP/NCLDV/")
otu <- read.delim('NCLDV_OTU.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
comm <- t(otu)
tree <- read.tree("GV64.iqtree.contree")
tax <- read.delim('taxanomy.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
save.wd=("Fig. 6/iCAMP/NCLDV/")

icamp.out=icamp.big(comm=comm,tree=tree,pd.wd=save.wd,rand = 1000,ds = 0.2,
                    bin.size.limit=5,phylo.rand.scale ="across.all",
                    taxa.rand.scale = "across.all",
                    sig.index = "SES.RC",phylo.metric ="bMPD",
                    nworker = 28,ses.cut = 1.96, rc.cut = 0.95,
                    taxo.metric="bray")

#summarize each bin
icampbin=icamp.bins(icamp.detail = icamp.out, treat = NULL,
                    clas = tax, boot=F,
                    rand.time = 1000, between.group = F)

#summarize each taxon
cate=icamp.cate(icamp.bins.result = icampbin,comm = comm,cate =tax,
                silent = FALSE,between.group = F)

#plot Fig. 6d summary iCAMP circle
#NCLDV
#overall
setwd("Fig. 6/iCAMP/")
data<-data.frame(variable=c(13.2,5.4,19.3,3.1,59), group = paste0("a", 1:5))
ggplot(data, aes(x = 3, y = variable, fill = group))+ geom_col() +
  coord_polar(theta = "y") +xlim(c(1.5, 4.5))+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("a1","a2","a3","a4","a5"),
                    values=c("#FFD9A6","#FEC3B9","#B3E2CD","#DDBBEA","#A8CEFF"))

#group
data<-read_excel("Fig. 6d.xlsx")
data$Process <- factor(data$Process,level=c("Hes","Hos","DL","HD","DR"))
data$Tax <- factor(data$Tax,level=c("Asfuvirales","Imitervirales","Algavirales","Pandoravirales","Pimascovirales"))
ggplot(data,aes(x=Tax,y=Proportion,fill=Process))+
  geom_bar(stat="identity", color="NA",width=0.8,size=0.2,position="stack", aes(fill=Process),alpha=0.8)+
  theme_bw()+labs(x="",y="Relative importance (%)")+
  scale_fill_manual(breaks=c("Hes","Hos","DL","HD","DR"),values = c("#FFD9A6","#FEC3B9","#B3E2CD","#DDBBEA","#A8CEFF"))+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(hjust = 1,size=15,angle = 20,colour = 'black'),
        axis.text.y=element_text(size=15,colour = 'black'),
        panel.border = element_rect(size=1),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.15, "cm"),
        legend.position = "none")

#Phage
data<-data.frame(variable=c(7.8,2.4,10.4,4.4,75), group = paste0("a", 1:5))
ggplot(data, aes(x = 3, y = variable, fill = group))+ geom_col() +
  coord_polar(theta = "y") +xlim(c(1.5, 4.5))+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("a1","a2","a3","a4","a5"),
                    values=c("#FFD9A6","#FEC3B9","#B3E2CD","#DDBBEA","#A8CEFF"))

#group
data<-read_excel("Fig. 6d.xlsx")
data$Process <- factor(data$Process,level=c("Hes","Hos","DL","HD","DR"))
data$Tax <- factor(data$Tax,level=c("Biggiephage","Jabbarphage","Judaphage","Whopperphage","Mahaphage"))
ggplot(data,aes(x=Tax,y=Proportion,fill=Process))+
  geom_bar(stat="identity", color="NA",width=0.8,size=0.2,position="stack", aes(fill=Process),alpha=0.8)+
  theme_bw()+labs(x="",y="Relative importance (%)")+
  scale_fill_manual(breaks=c("Hes","Hos","DL","HD","DR"),values = c("#FFD9A6","#FEC3B9","#B3E2CD","#DDBBEA","#A8CEFF"))+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(hjust = 1,size=14,angle = 20,colour = 'black'),
        axis.text.y=element_text(size=14,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        legend.position = "none")