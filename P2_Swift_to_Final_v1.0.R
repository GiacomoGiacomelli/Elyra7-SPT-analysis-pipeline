###PART 2######

#TARDIS
##Input
###"output_TARDIS.csv"
##BasicTab
###MaxJD->2µm
###MaxDt->3
###FrameTime->24ms
###LongestTrackLength->100 frames
###LocalizationUncert->30nm
###BGsubUnc->Slow&exact
###Auto-choose pops **** or 2 populations
###Fit w/bleach
##StorageTab
###Store SWIFT parameters

#SWIFT
##Input
###"output_TARDIS.csv"
##Tracking --> Parameters
###Load from TARDIS output folder
##PerformanceParameters
###Max Displacement->3000nm
###Max Displacement(postprocessing)->3000nm
###Maximum Memory ->2000.00 MB
##Analysis
###Calculate JD
###Calculate MJD
###Calculate MSD
###Calculate D by MSD
###Calculate Start time and duration
###Identify Motion Type
###Identify Noise

##Output
###Save Swift file
###Export all info about localizations, segments, tracks (in this order, as a single file)
#########


###PART 3######
library(spatstat)
library(ggplot2)
library(stringr)
library(fields)  
library(tidyverse)
library(conflicted)
library(patchwork)

folder<-"F:/Giacomo/Dips_paper/Analysis/"
setwd(folder)  #Choose the folder containing all strains folders
MAIN<-getwd()
SUBS<-list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

for (z in SUBS) {
  fileTname <- Sys.glob(paste(MAIN,"/",z,"/output_TARDIS.export.csv",sep=""))
  fileTname1 <- Sys.glob(paste(MAIN,"/",z,"/output.csv",sep=""))
  file_tracks<-read.csv(fileTname, header=T)
  file_together<-read.csv(fileTname1, header=T)
  
  file_tracks$id<-NULL
  file_tracks$chi2<-NULL
  file_tracks$sigma..nm. <- NULL
  file_tracks$uncertainty_xy..nm.<- NULL
  file_tracks$intensity..photon.<-NULL
  file_tracks$offset..photon.<-NULL
  file_tracks$bkgstd..photon.<-NULL
  nc1<-ncol(file_together)
  nc2<-ncol(file_tracks)
  
  file_tracks[(nc2+1):(nc1+nc2+1)]<-NA
  ncol(file_tracks)
  colnames(file_tracks)[(nc2+1):(nc1+nc2)]<-names(file_together)
  colnames(file_tracks)[(nc1+nc2+1)]<-"Repeated"
  file_tracks[(nc1+nc2+1)]<-"no"
  file_tracks$CellName<-"noCell"
  file_tracks[c(1:(length(file_together$id))),c((nc2+1):(nc1+nc2))]<-file_together  

if (length(file_together$id)!=length(file_tracks$track.id)){
  file_tracks[c((length(file_together$id)+1):length(file_tracks$track.id)),]$Repeated<-"yes"
}
  
for (k in (length(file_together$id)+1):length(file_tracks$track.id)){
  if (is.na(file_together[round(file_together$x..nm.)==round(file_tracks[k,]$x) & round(file_together$y..nm.)==round(file_tracks[k,]$y) & file_together$frame==file_tracks[k,]$t,][,15][1])){
    print(k)
  }
  else {
    file_tracks[k,c((nc2+1):(nc1+nc2))]<-file_together[round(file_together$x..nm.)==round(file_tracks[k,]$x) & round(file_together$y..nm.)==round(file_tracks[k,]$y) & file_together$frame==file_tracks[k,]$t,c(1:(nc1-1))]
  }
}
file_tracks[c((length(file_together$id)+1):length(file_tracks$track.id)),]$Unique<-"no"

####all this should be removed ahead of time
file_tracks$x..nm. <- NULL
file_tracks$y..nm. <- NULL

####This is the file that contains all the info you need
write.csv(file_tracks[file_tracks$CellName!="noCell",], file=paste(MAIN,"/",z,"/FinalOutput.csv",sep=""), row.names = FALSE)

####From here on we will get plots and calculate bubbles and diffusion coefficients.

c<-ggplot(data=file_tracks[file_tracks$Protein_Position=="Betwixt" & file_tracks$Repeated=="no",], aes(x=Position.X.betwixt, y=Position.Y.N))+
  geom_bin2d(binwidth=25)+
  scale_fill_continuous(type = "viridis")+  #betwixt get double signal compared to the rest, so it needs to be visually normalized
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank())+
  xlab("Betwixt (nm)")

a<-ggplot(data=file_tracks[file_tracks$Protein_Position=="Major_pole"& file_tracks$Repeated=="no",], aes(x=Position.X.polar, y=Position.Y.N))+
  geom_bin2d(binwidth=25)+
  scale_fill_continuous(type = "viridis")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Pole_M (nm)")+
  ylab("Width (nm)")

d<-ggplot(data=file_tracks[file_tracks$Protein_Position=="Midcell" & file_tracks$Repeated=="no",], aes(x=Position.X.N, y=Position.Y.N))+
  geom_bin2d(binwidth=25)+
  scale_fill_continuous(type = "viridis")+
  coord_fixed()+
  theme_bw()+
  theme(legend.text = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())+
  xlab("Midcell (nm)")

b<-ggplot(data=file_tracks[file_tracks$Protein_Position=="Minor_pole" & file_tracks$Repeated=="no",], aes(x=Position.X.polar, y=Position.Y.N))+
  geom_bin2d(binwidth=25)+
  scale_fill_continuous(type = "viridis")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank())+
  xlab("Pole_m (nm)")

e<-ggplot(data=file_tracks, aes(x=Position.X.N, y=Position.Y.N))+
  geom_bin2d(binwidth=25)+
  scale_fill_continuous(type = "viridis")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank())+
  xlab("Pole_m (nm)")

f<-ggplot(data=file_tracks[file_tracks$seg.motion=="diffusion" & file_tracks$seg.loc_count>4,], aes(x=Position.X.N, y=Position.Y.N))+
  geom_bin2d(binwidth=25)+
  scale_fill_continuous(type = "viridis")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank())+
  xlab("Pole_m (nm)")

g<-ggplot(data=file_tracks[file_tracks$seg.motion=="immobile" & file_tracks$seg.loc_count>4,], aes(x=Position.X.N, y=Position.Y.N))+
  geom_bin2d(binwidth=25)+
  scale_fill_continuous(type = "viridis")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank())+
  xlab("Pole_m (nm)")

h<-ggplot(data=file_tracks[file_tracks$CellLength_SMLM<2000 & file_tracks$seg_pos==5 & file_tracks$seg.mjd>0 & file_tracks$seg.motion =="diffusion",], aes(x=seg.mjd/1000, weight=seg.mjd_n, y=(length(file_tracks[file_tracks$seg_pos==4 & file_tracks$seg.mjd>0 & file_tracks$seg.motion =="diffusion",]$track.id)/length(file_tracks[file_tracks$seg_pos==5 & file_tracks$seg.mjd>0,]$track.id))*(after_stat(density))))+ 
  geom_histogram(binwidth=0.01, fill="#03045e", alpha=0.5)+
  theme_bw()+
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.05))+
  ylab("Density")+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none")

i<-ggplot(data=file_tracks[file_tracks$seg_pos==5 & file_tracks$seg.mjd>0 & file_tracks$seg.motion =="immobile",], aes(x=seg.mjd/1000, weight=seg.mjd_n, y=(length(file_tracks[file_tracks$seg_pos==4 & file_tracks$seg.mjd>0 & file_tracks$seg.motion =="immobile",]$track.id)/length(file_tracks[file_tracks$seg_pos==5 & file_tracks$seg.mjd>0,]$track.id))*(after_stat(density))))+ 
  geom_histogram(binwidth=0.01, fill="red", alpha=0.5)+
  theme_bw()+
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.05))+
  ylab("Density")+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none")

f1<-ggplot(data=file_tracks[file_tracks$seg.motion=="diffusion" & file_tracks$seg.mjd<=125 & file_tracks$seg.loc_count>4,], aes(x=Position.X.N, y=Position.Y.N))+
  geom_bin2d(binwidth=25)+
  scale_fill_continuous(type = "viridis")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank())+
  xlab("Pole_m (nm)")

f3<-ggplot(data=file_tracks[file_tracks$seg.motion=="diffusion" & file_tracks$seg.mjd>250 & file_tracks$seg.loc_count>4,], aes(x=Position.X.N, y=Position.Y.N))+
  geom_bin2d(binwidth=25)+
  scale_fill_continuous(type = "viridis")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank())+
  xlab("Pole_m (nm)")

f2<-ggplot(data=file_tracks[file_tracks$seg.motion=="diffusion" & file_tracks$seg.mjd>125 & file_tracks$seg.mjd<=250 & file_tracks$seg.loc_count>4,], aes(x=Position.X.N, y=Position.Y.N))+
  geom_bin2d(binwidth=25)+
  scale_fill_continuous(type = "viridis")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank())+
  xlab("Pole_m (nm)")

pdf(paste(MAIN,"/",z,"/HeatMap_subcellular_",z,".pdf",sep=""), height = 5, width = 5)
print(a|b|c|d)
dev.off()

pdf(paste(MAIN,"/",z,"/HeatMap_",z,".pdf",sep=""), height = 5, width = 5)
print(e)
dev.off()

pdf(paste(MAIN,"/",z,"/HeatMap_diff_imm",z,".pdf",sep=""), height = 5, width = 5)
print((f/g)|(h/i))
dev.off()

pdf(paste(MAIN,"/",z,"/HeatMap_diff_split",z,".pdf",sep=""), height = 5, width = 5)
print(f1/f2/f3)
dev.off()
}


##PART 6
#After you have all conditions analysed up to PART5 ###
library("ggplot2")
library("plotrix")
library("fields")                                           
library("spatstat")
library("pgirmess")
library("RColorBrewer")
library("patchwork")
library("ggbeeswarm")
library("remotes")
#remotes::install_github("finnlindgren/StatCompLab") #run only the first time
library("StatCompLab")
######FITTING
library("minpack.lm")
library("Rmixmod")
# Generate gamma rvs
library("broom")

folder<-"F:/Giacomo/Dips_paper/Analysis/"
setwd(folder)  #Choose the folder containing all strains folders
MAIN<-getwd()
SUBS<-list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
dt<-0.024  ###integration + transfer (total frame time)


for (z in 1:length(SUBS)) {
  if (z==1){
    jdT<-read.table(paste(MAIN,"/",SUBS[z],"/FinalOutput.csv",sep=""), header=T, sep=",")
    cat(c(paste("strain:",z,sep=""),capture.output(SUBS[z])), file= "Strains_order.txt")
  } else {
    cat(c(paste("\n","strain:",z,sep=""),capture.output(SUBS[z])), file = "Strains_order.txt", append = TRUE)
  addition<-read.table(paste(MAIN,"/",SUBS[z],"/FinalOutput.csv",sep=""), header=T, sep=",")
  jdT<-rbind(jdT,addition)
  }
}

jdT<-jdT[!is.na(jdT$Repeated),]
j<-jdT[jdT$seg.loc_count>=4,]
j<-as.data.frame(j[j$j>0 & !is.na(j$j),]$j)
colnames(j)<-"Jump"
j$Jump<-j$Jump/1000

P = ecdf(j$Jump)
z = seq(0, 2.5, by=0.001)
p = P(z)
cj<-as.data.frame(z)
cj[2]<-p
colnames(cj)[1]<-"r"
colnames(cj)[2]<-"CD"

#plot jump distance
hist(j$Jump, xlim=c(0,2.5), breaks = 250)

#plot cumulative distribution
plot(cj)


###################################################################################

histdataT<-hist(j$Jump, xlim=c(0,2.5), breaks = 250)
forfitT<-as.data.frame(histdataT$mids)
colnames(forfitT)[1]<-"x"
forfitT[2]<-histdataT$density
colnames(forfitT)[2]<-"y"
plot(forfitT)

hist(j$Jump, xlim=c(0,2.5), breaks = 250)

fun1_T = function(x,sigma){(x/sigma^2)*exp((-x^2)/(2*(sigma^2)))}
fun2_T = function(x,sigma1,sigma2,f1){f1*(x/sigma1^2)*exp((-x^2)/(2*(sigma1^2)))+(1-f1)*(x/sigma2^2)*exp((-x^2)/(2*(sigma2^2)))}
fun3_T = function(x,sigma1,sigma2,sigma3,f1,f2,f3){
  f3<-1-f1-f2
  if (f3<0) return(rep(NA,length(x)))
  (1-f2-f3)*(x/sigma1^2)*exp((-x^2)/(2*(sigma1^2)))+(1-f1-f3)*(x/sigma2^2)*exp((-x^2)/(2*(sigma2^2)))+(1-f1-f2)*(x/sigma3^2)*exp((-x^2)/(2*(sigma3^2)))}

fun4_T = function(x,sigma1,sigma2,sigma3,sigma4,f1,f2,f3,f4){
  f4<-1-f1-f2-f3
  if (f4<0) return(rep(NA,length(x)))
  (1-f2-f3-f4)*(x/sigma1^2)*exp((-x^2)/(2*(sigma1^2)))+(1-f1-f3-f4)*(x/sigma2^2)*exp((-x^2)/(2*(sigma2^2)))+(1-f1-f2-f4)*(x/sigma3^2)*exp((-x^2)/(2*(sigma3^2)))+(1-f1-f2-f3)*(x/sigma4^2)*exp((-x^2)/(2*(sigma4^2)))}



wT<-mixmodCluster(j$Jump,1)
fit1s1_T<-as.numeric(wT["bestResult"][6][2])
fit1_T <- nlsLM(y~fun1_T(x, sigma),forfitT,
                start=list(sigma=fit1s1_T))
hist(wT, xlim=c(0,2.5), breaks = 250)

wT<-mixmodCluster(j$Jump,2)
fit2p1_T<-as.numeric(wT["bestResult"][6][1][1])
fit2p2_T<-as.numeric(wT["bestResult"][6][1][2])
fit2s1_T<-as.numeric(wT["bestResult"][6][2][1])
fit2s2_T<-as.numeric(wT["bestResult"][6][2][2])
fit2_T <- nlsLM(y~fun2_T(x, sigma1,sigma2,f1),forfitT,
                start=list(sigma1=fit2s1_T, sigma2=fit2s2_T, f1=fit2p1_T), lower=c(sigma1=0, sigma2=0, f1=0), upper=c(sigma1=10, sigma2=10, f1=1))
hist(wT, xlim=c(0,2.5), breaks = 250)

wT<-mixmodCluster(j$Jump,3)
fit3p1_T<-as.numeric(wT["bestResult"][6][1][1])
fit3p2_T<-as.numeric(wT["bestResult"][6][1][2])
fit3p3_T<-as.numeric(wT["bestResult"][6][1][3])
fit3s1_T<-as.numeric(wT["bestResult"][6][2][1])
fit3s2_T<-as.numeric(wT["bestResult"][6][2][2])
fit3s3_T<-as.numeric(wT["bestResult"][6][2][3])
fit3_T <- nlsLM(y~fun3_T(x, sigma1,sigma2,sigma3,f1,f2),forfitT,
                start=list(sigma1=fit3s1_T, sigma2=fit3s2_T, sigma3=fit3s3_T, f1=fit3p1_T, f2=fit3p2_T), lower=c(sigma1=0, sigma2=0, sigma3=0, f1=0, f2=0), upper=c(sigma1=10, sigma2=10, sigma3=10, f1=1, f2=1))
hist(wT, xlim=c(0,2.5), breaks = 250)


wT<-mixmodCluster(j$Jump,4)
fit4p1_T<-as.numeric(wT["bestResult"][6][1][1])
fit4p2_T<-as.numeric(wT["bestResult"][6][1][2])
fit4p3_T<-as.numeric(wT["bestResult"][6][1][3])
fit4p4_T<-as.numeric(wT["bestResult"][6][1][4])
fit4s1_T<-as.numeric(wT["bestResult"][6][2][1])
fit4s2_T<-as.numeric(wT["bestResult"][6][2][2])
fit4s3_T<-as.numeric(wT["bestResult"][6][2][3])
fit4s4_T<-as.numeric(wT["bestResult"][6][2][4])
fit4_T <- nlsLM(y~fun4_T(x, sigma1,sigma2,sigma3,sigma4,f1,f2,f3),forfitT,
                start=list(sigma1=fit4s1_T, sigma2=fit4s2_T, sigma3=fit4s3_T, sigma4=fit4s4_T, f1=fit4p1_T, f2=fit4p2_T, f3=fit4p3_T), lower=c(sigma1=0, sigma2=0, sigma3=0, sigma4=0, f1=0, f2=0, f3=0), upper=c(sigma1=10, sigma2=10, sigma3=10, sigma4=10, f1=1, f2=1, f3=1))
hist(wT, xlim=c(0,2.5), breaks = 250)

###absolute values for sigmas
sigma1T_1<-abs(summary(fit1_T)$parameters[1])
sigma2T_1<-abs(summary(fit2_T)$parameters[1])
sigma2T_2<-abs(summary(fit2_T)$parameters[2])
sigma3T_1<-abs(summary(fit3_T)$parameters[1])
sigma3T_2<-abs(summary(fit3_T)$parameters[2])
sigma3T_3<-abs(summary(fit3_T)$parameters[3])
sigma4T_1<-abs(summary(fit4_T)$parameters[1])
sigma4T_2<-abs(summary(fit4_T)$parameters[2])
sigma4T_3<-abs(summary(fit4_T)$parameters[3])
sigma4T_4<-abs(summary(fit4_T)$parameters[4])

eq1_T = function(x){(x/sigma1T_1^2)*exp((-x^2)/(2*(sigma1T_1^2)))}

eq2_T = function(x){summary(fit2_T)$parameters[3]*(x/sigma2T_1^2)*exp((-x^2)/(2*(sigma2T_1^2)))+(1-summary(fit2_T)$parameters[3])*(x/sigma2T_2^2)*exp((-x^2)/(2*(sigma2T_2^2)))}
eq2_1_T = function(x){summary(fit2_T)$parameters[3]*(x/sigma2T_1^2)*exp((-x^2)/(2*(sigma2T_1^2)))}
eq2_2_T = function(x){(1-summary(fit2_T)$parameters[3])*(x/sigma2T_2^2)*exp((-x^2)/(2*(sigma2T_2^2)))}

eq3_T = function(x){(1-summary(fit3_T)$parameters[5]-(1-summary(fit3_T)$parameters[4]-summary(fit3_T)$parameters[5]))*(x/sigma3T_1^2)*exp((-x^2)/(2*(sigma3T_1^2)))+(1-summary(fit3_T)$parameters[4]-(1-summary(fit3_T)$parameters[4]-summary(fit3_T)$parameters[5]))*(x/sigma3T_2^2)*exp((-x^2)/(2*(sigma3T_2^2)))+(1-summary(fit3_T)$parameters[5]-summary(fit3_T)$parameters[4])*(x/sigma3T_3^2)*exp((-x^2)/(2*(sigma3T_3^2)))}
eq3_1_T = function(x){(1-summary(fit3_T)$parameters[5]-(1-summary(fit3_T)$parameters[4]-summary(fit3_T)$parameters[5]))*(x/sigma3T_1^2)*exp((-x^2)/(2*(sigma3T_1^2)))}
eq3_2_T = function(x){(1-summary(fit3_T)$parameters[4]-(1-summary(fit3_T)$parameters[4]-summary(fit3_T)$parameters[5]))*(x/sigma3T_2^2)*exp((-x^2)/(2*(sigma3T_2^2)))}
eq3_3_T = function(x){(1-summary(fit3_T)$parameters[5]-summary(fit3_T)$parameters[4])*(x/sigma3T_3^2)*exp((-x^2)/(2*(sigma3T_3^2)))}

eq4_T = function(x){(1-summary(fit4_T)$parameters[6]-summary(fit4_T)$parameters[7]-(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[6]-summary(fit4_T)$parameters[7]))*(x/sigma4T_1^2)*exp((-x^2)/(2*(sigma4T_1^2)))+(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[7]-(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[6]-summary(fit4_T)$parameters[7]))*(x/sigma4T_2^2)*exp((-x^2)/(2*(sigma4T_2^2)))+(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[6]-(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[6]-summary(fit4_T)$parameters[7]))*(x/sigma4T_3^2)*exp((-x^2)/(2*(sigma4T_3^2)))+(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[6]-summary(fit4_T)$parameters[7])*(x/sigma4T_4^2)*exp((-x^2)/(2*(sigma4T_4^2)))}
eq4_1_T = function(x){(1-summary(fit4_T)$parameters[6]-summary(fit4_T)$parameters[7]-(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[6]-summary(fit4_T)$parameters[7]))*(x/sigma4T_1^2)*exp((-x^2)/(2*(sigma4T_1^2)))}
eq4_2_T = function(x){(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[7]-(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[6]-summary(fit4_T)$parameters[7]))*(x/sigma4T_2^2)*exp((-x^2)/(2*(sigma4T_2^2)))}
eq4_3_T = function(x){(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[6]-(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[6]-summary(fit4_T)$parameters[7]))*(x/sigma4T_3^2)*exp((-x^2)/(2*(sigma4T_3^2)))}
eq4_4_T = function(x){(1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[6]-summary(fit4_T)$parameters[7])*(x/sigma4T_4^2)*exp((-x^2)/(2*(sigma4T_4^2)))}



#c("#03045e", "#353f84", "#6477a9", "#96b2cf", "#caf0f8")

pop1_Tot<-ggplot(data=forfitT, aes(x=x, y=y))+ 
  geom_line(col="red", linewidth=1)+
  stat_function(fun=eq1_T, col="#03045e", linewidth =1)+
  theme_bw()+
  xlim(0,1.5)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  ylab("Density")

pop2_Tot<-ggplot(data=forfitT, aes(x=x, y=y))+ 
  geom_line(col="red", linewidth=1)+
  stat_function(fun=eq2_T, col="#03045e", linewidth =1)+
  stat_function(fun=eq2_1_T, col="#353f84", linewidth =1, linetype = "dashed")+
  stat_function(fun=eq2_2_T, col="#6477a9", linewidth =1, linetype = "dashed")+
  theme_bw()+
  xlim(0,1.5)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  ylab("Density")

pop3_Tot<-ggplot(data=forfitT, aes(x=x, y=y))+ 
  geom_line(col="red", linewidth=1)+
  stat_function(fun=eq3_T, col="#03045e", linewidth =1)+
  stat_function(fun=eq3_1_T, col="#353f84", linewidth =1, linetype = "dashed")+
  stat_function(fun=eq3_2_T, col="#6477a9", linewidth =1, linetype = "dashed")+
  stat_function(fun=eq3_3_T, col="#96b2cf", linewidth =1, linetype = "dashed")+
  theme_bw()+
  xlim(0,1.5)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  ylab("Density")

pop4_Tot<-ggplot(data=forfitT, aes(x=x, y=y))+ 
  geom_line(col="red", linewidth=1)+
  stat_function(fun=eq4_T, col="#03045e", linewidth =1)+
  stat_function(fun=eq4_1_T, col="#353f84", linewidth =1, linetype = "dashed")+
  stat_function(fun=eq4_2_T, col="#6477a9", linewidth =1, linetype = "dashed")+
  stat_function(fun=eq4_3_T, col="#96b2cf", linewidth =1, linetype = "dashed")+
  stat_function(fun=eq4_4_T, col="#caf0f8", linewidth =1, linetype = "dashed")+
  theme_bw()+
  xlim(0,1.5)+
  ylab("Density")+
  xlab("Jump distance [µm]")

pdf("Total_JD_fitting_B4.pdf", height = 5, width = 5)
pop1_Tot/pop2_Tot/pop3_Tot/pop4_Tot
dev.off()


A_1t4<-capture.output(anova(fit1_T,fit2_T,fit3_T,fit4_T))
A_1t2<-capture.output(anova(fit1_T,fit2_T))
A_2t3<-capture.output(anova(fit2_T,fit3_T))
A_3t4<-capture.output(anova(fit3_T,fit4_T))
A_bic<-capture.output(BIC(fit1_T,fit2_T,fit3_T,fit4_T))
A_aic<-capture.output(AIC(fit1_T,fit2_T,fit3_T,fit4_T))
A_g1<-capture.output(glance(fit1_T))
A_g2<-capture.output(glance(fit2_T))
A_g3<-capture.output(glance(fit3_T))
A_g4<-capture.output(glance(fit4_T))

writeLines(c("Anova 1 vs 4 pops.",A_1t4,"","Anova 1 vs 2 pops.",A_1t2,"","Anova 2 vs 3 pops.",A_2t3,"","Anova 3 vs 4 pops.",A_3t4,"","BIC analysis",A_bic,"","AIC analysis",A_aic,"","Summary 1 pop",A_g1,"","Summary 2 pop",A_g2,"","Summary 3 pop",A_g3,"","Summary 4 pop",A_g4), "anova_output_ALL.txt")

Fit1<-capture.output(fit1_T)
Fit2<-capture.output(fit2_T)
Fit3<-capture.output(fit3_T)
fit1<-capture.output(fit4_T)
writeLines(c(Fit1, Fit2, Fit3, fit1), "fits.txt")

####Bubbles general
loc_error<-mean(jdT$uncertainty_xy)
sigma1<-abs(summary(fit1_T)$parameters[1])
D1<-(9/(13*dt))*(sigma1^2-(2*((loc_error/1000)^2)))
MJ1<-sigma1*sqrt(pi/2)
proportion1<-1

sigma2<-abs(summary(fit2_T)$parameters[1:2])
D2<-(9/(13*(dt)))*(sigma2^2-(2*(((loc_error)/1000)^2)))
MJ2<-sigma2*sqrt(pi/2)
proportion2<-c(summary(fit2_T)$parameters[3], 1-summary(fit2_T)$parameters[3])

sigma3<-abs(summary(fit3_T)$parameters[1:3])
D3<-(9/(13*(dt)))*(sigma3^2-(2*(((loc_error)/1000)^2)))
MJ3<-sigma3*sqrt(pi/2)
proportion3<-c(summary(fit3_T)$parameters[4:5],1-summary(fit3_T)$parameters[4]-summary(fit3_T)$parameters[5])

sigma4<-abs(summary(fit4_T)$parameters[1:4])
D4<-(9/(13*(dt)))*(sigma4^2-(2*(((loc_error)/1000)^2)))
MJ4<-sigma4*sqrt(pi/2)
proportion4<-c(summary(fit4_T)$parameters[5:7],1-summary(fit4_T)$parameters[5]-summary(fit4_T)$parameters[6]-summary(fit4_T)$parameters[7])

for (i in 1:length(D1)){
  if (D1[i]<0)
  {
    D1[i]<-0.0001
  }  
}

for (i in 1:length(D2)){
  if (D2[i]<0)
  {
    D2[i]<-0.0001
  }  
}

for (i in 1:length(D3)){
  if (D3[i]<0)
  {
    D3[i]<-0.0001
  }  
}

for (i in 1:length(D4)){
  if (D4[i]<0)
  {
    D4[i]<-0.0001
  }  
}

bubbles1<-as.data.frame(cbind(D1,sigma1,proportion1,MJ1))
bubbles1[5]<-1
colnames(bubbles1)<-c("Diffusion_Adjusted","Sigma","proportion","meanMJD","Components")
bubbles2<-as.data.frame(cbind(D2,sigma2,proportion2,MJ2))
bubbles2[5]<-2
colnames(bubbles2)<-c("Diffusion_Adjusted","Sigma","proportion","meanMJD","Components")
bubbles3<-as.data.frame(cbind(D3,sigma3,proportion3,MJ3))
bubbles3[5]<-3
colnames(bubbles3)<-c("Diffusion_Adjusted","Sigma","proportion","meanMJD","Components")
bubbles4<-as.data.frame(cbind(D4,sigma4,proportion4,MJ4))
bubbles4[5]<-4
colnames(bubbles4)<-c("Diffusion_Adjusted","Sigma","proportion","meanMJD","Components")
jump_bubbles<-rbind(bubbles1,bubbles2, bubbles3, bubbles4)
jump_bubbles[6]<-"MERGED"
colnames(jump_bubbles)[6]<-"Strain"
###bubble general end

for (z in 1:(length(SUBS))) {
  j<-jdT[jdT$Strain==SUBS[z] & jdT$seg.loc_count>=4,]
  j<-as.data.frame(j[j$j>0,]$j)
  colnames(j)<-"Jump"
  j$Jump<-j$Jump/1000 # only if jump is expressed in nanometers
  
  P1 = ecdf(j$Jump)
  z1 = seq(0, 2.5, by=0.001)
  p1 = P1(z1)
  cj1<-as.data.frame(z1)
  cj1[2]<-p1
  colnames(cj1)[1]<-"r"
  colnames(cj1)[2]<-"CD"
  
  histdata1<-hist(j$Jump, xlim=c(0,2.5), breaks = 250)
  forfit1<-as.data.frame(histdata1$mids)
  colnames(forfit1)[1]<-"x"
  forfit1[2]<-histdata1$density
  colnames(forfit1)[2]<-"y"
  plot(forfit1)
  
  fun1A = function(x,f1){f1*(x/(sigma1T_1)^2)*exp((-x^2)/(2*(sigma1T_1^2)))}
  fun2A = function(x,f1){f1*(x/sigma2T_1^2)*exp((-x^2)/(2*(sigma2T_1^2)))+(1-f1)*(x/sigma2T_2^2)*exp((-x^2)/(2*(sigma2T_2^2)))}
  fun3A = function(x,f1,f2,f3){
    f3<-1-f1-f2
    if (f3<0) return(rep(NA,length(x)))
    (1-f2-f3)*(x/sigma3T_1^2)*exp((-x^2)/(2*(sigma3T_1^2)))+(1-f1-f3)*(x/sigma3T_2^2)*exp((-x^2)/(2*(sigma3T_2^2)))+(1-f1-f2)*(x/sigma3T_3^2)*exp((-x^2)/(2*(sigma3T_3^2)))}
  
  fun4A = function(x,f1,f2,f3,f4){
    f4<-1-f1-f2-f3
    if (f4<0) return(rep(NA,length(x)))
    (1-f2-f3-f4)*(x/sigma4T_1^2)*exp((-x^2)/(2*(sigma4T_1^2)))+(1-f1-f3-f4)*(x/sigma4T_2^2)*exp((-x^2)/(2*(sigma4T_2^2)))+(1-f1-f2-f4)*(x/sigma4T_3^2)*exp((-x^2)/(2*(sigma4T_3^2)))+(1-f1-f2-f3)*(x/sigma4T_4^2)*exp((-x^2)/(2*(sigma4T_4^2)))}
  
  
  fit2A <- nlsLM(y~fun2A(x, f1),forfit1,
                 start=list(f1=summary(fit2_T)$parameters[3]), lower=c(f1=0), upper=c(f1=1))
  
  fit3A <- nlsLM(y~fun3A(x,f1,f2),forfit1,
                 start=list(f1=summary(fit3_T)$parameters[4], f2=summary(fit3_T)$parameters[5]), lower=c( f1=0, f2=0), upper=c(f1=1, f2=1))
  
  fit1A <- nlsLM(y~fun4A(x,f1,f2,f3),forfit1,
                 start=list(f1=summary(fit4_T)$parameters[5], f2=summary(fit4_T)$parameters[6], f3=summary(fit4_T)$parameters[7]), lower=c( f1=0, f2=0, f3=0), upper=c(f1=1, f2=1, f3=1))
  
  eq2A = function(x){summary(fit2A)$parameters[1]*(x/sigma2T_1^2)*exp((-x^2)/(2*(sigma2T_1^2)))+(1-summary(fit2A)$parameters[1])*(x/sigma2T_2^2)*exp((-x^2)/(2*(sigma2T_2^2)))}
  eq2_1A = function(x){summary(fit2A)$parameters[1]*(x/sigma2T_1^2)*exp((-x^2)/(2*(sigma2T_1^2)))}
  eq2_2A = function(x){(1-summary(fit2A)$parameters[1])*(x/sigma2T_2^2)*exp((-x^2)/(2*(sigma2T_2^2)))}
  
  eq3A   = function(x){(1-summary(fit3A)$parameters[2]-(1-summary(fit3A)$parameters[1]-summary(fit3A)$parameters[2]))*(x/sigma3T_1^2)*exp((-x^2)/(2*(sigma3T_1^2)))+(1-summary(fit3A)$parameters[1]-(1-summary(fit3A)$parameters[1]-summary(fit3A)$parameters[2]))*(x/sigma3T_2^2)*exp((-x^2)/(2*(sigma3T_2^2)))+(1-summary(fit3A)$parameters[2]-summary(fit3A)$parameters[1])*(x/sigma3T_3^2)*exp((-x^2)/(2*(sigma3T_3^2)))}
  eq3_1A = function(x){(1-summary(fit3A)$parameters[2]-(1-summary(fit3A)$parameters[1]-summary(fit3A)$parameters[2]))*(x/sigma3T_1^2)*exp((-x^2)/(2*(sigma3T_1^2)))}
  eq3_2A = function(x){(1-summary(fit3A)$parameters[1]-(1-summary(fit3A)$parameters[1]-summary(fit3A)$parameters[2]))*(x/sigma3T_2^2)*exp((-x^2)/(2*(sigma3T_2^2)))}
  eq3_3A = function(x){(1-summary(fit3A)$parameters[2]-summary(fit3A)$parameters[1])*(x/sigma3T_3^2)*exp((-x^2)/(2*(sigma3T_3^2)))}
  
  eq4A = function(x){(1-summary(fit1A)$parameters[2]-summary(fit1A)$parameters[3]-(1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[2]-summary(fit1A)$parameters[3]))*(x/sigma4T_1^2)*exp((-x^2)/(2*(sigma4T_1^2)))+
      (1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[3]-(1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[2]-summary(fit1A)$parameters[3]))*(x/sigma4T_2^2)*exp((-x^2)/(2*(sigma4T_2^2)))+
      (1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[2]-(1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[2]-summary(fit1A)$parameters[3]))*(x/sigma4T_3^2)*exp((-x^2)/(2*(sigma4T_3^2)))+
      (1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[2]-summary(fit1A)$parameters[3])*(x/sigma4T_4^2)*exp((-x^2)/(2*(sigma4T_4^2)))}
  eq4_1A = function(x){(1-summary(fit1A)$parameters[2]-summary(fit1A)$parameters[3]-(1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[2]-summary(fit1A)$parameters[3]))*(x/sigma4T_1^2)*exp((-x^2)/(2*(sigma4T_1^2)))}
  eq4_2A = function(x){(1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[3]-(1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[2]-summary(fit1A)$parameters[3]))*(x/sigma4T_2^2)*exp((-x^2)/(2*(sigma4T_2^2)))}
  eq4_3A = function(x){(1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[2]-(1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[2]-summary(fit1A)$parameters[3]))*(x/sigma4T_3^2)*exp((-x^2)/(2*(sigma4T_3^2)))}
  eq4_4A = function(x){(1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[2]-summary(fit1A)$parameters[3])*(x/sigma4T_4^2)*exp((-x^2)/(2*(sigma4T_4^2)))}
  
  pop1_A<-ggplot(data=forfit1, aes(x=x, y=y))+ 
    geom_line(col="red", linewidth=1)+
    stat_function(fun=eq1_T, col="#03045e", linewidth =1)+
    theme_bw()+
    xlim(0,1.5)+
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
    ylab("Density")
  
  pop2_A<-ggplot(data=forfit1, aes(x=x, y=y))+ 
    geom_line(col="red", linewidth=1)+
    stat_function(fun=eq2A, col="#03045e", linewidth=1)+
    stat_function(fun=eq2_1A, col="#353f84", linewidth=1, linetype = "dashed")+
    stat_function(fun=eq2_2A, col="#6477a9", linewidth=1, linetype = "dashed")+
    theme_bw()+
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
    xlim(0,1.5)+
    ylab("Density")
  
  pop3_A<-ggplot(data=forfit1, aes(x=x, y=y))+ 
    geom_line(col="red", linewidth=1)+
    stat_function(fun=eq3A, col="#03045e", linewidth=1)+
    stat_function(fun=eq3_1A, col="#353f84", linewidth=1, linetype = "dashed")+
    stat_function(fun=eq3_2A, col="#6477a9", linewidth=1, linetype = "dashed")+
    stat_function(fun=eq3_3A, col="#96b2cf", linewidth=1, linetype = "dashed")+
    theme_bw()+
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
    xlim(0,1.5)+
    ylab("Density")
  
  pop4_A<-ggplot(data=forfit1, aes(x=x, y=y))+ 
    geom_line(col="red", linewidth=1)+
    stat_function(fun=eq4A, col="#03045e", linewidth=1)+
    stat_function(fun=eq4_1A, col="#353f84", linewidth=1, linetype = "dashed")+
    stat_function(fun=eq4_2A, col="#6477a9", linewidth=1, linetype = "dashed")+
    stat_function(fun=eq4_3A, col="#96b2cf", linewidth=1, linetype = "dashed")+
    stat_function(fun=eq4_4A, col="#caf0f8", linewidth=1, linetype = "dashed")+
    theme_bw()+
    xlim(0,1.5)+
    ylab("Density")+
    xlab("Jump distance [µm]")
  
  pdf(paste("s",unique(SUBS)[z],".pdf",sep=""), height = 5, width = 5)
  print(pop1_A/pop2_A/pop3_A/pop4_A)
  dev.off()
  
  proportion2A<-c(summary(fit2A)$parameters[1], 1-summary(fit2A)$parameters[1])
  proportion3A<-c(summary(fit3A)$parameters[1:2],1-summary(fit3A)$parameters[1]-summary(fit3A)$parameters[2])
  proportion4A<-c(summary(fit1A)$parameters[1:3],1-summary(fit1A)$parameters[1]-summary(fit1A)$parameters[2]-summary(fit1A)$parameters[3])

  bubbles2A<-as.data.frame(cbind(D2,sigma2,proportion2A,MJ2))
  bubbles2A[5]<-2
  colnames(bubbles2A)<-c("Diffusion_Adjusted","Sigma","proportion","meanMJD","Components")
  bubbles3A<-as.data.frame(cbind(D3,sigma3,proportion3A,MJ3))
  bubbles3A[5]<-3
  colnames(bubbles3A)<-c("Diffusion_Adjusted","Sigma","proportion","meanMJD","Components")
  bubbles4A<-as.data.frame(cbind(D4,sigma4,proportion4A,MJ4))
  bubbles4A[5]<-4
  colnames(bubbles4A)<-c("Diffusion_Adjusted","Sigma","proportion","meanMJD","Components")
  jump_bubblesA<-rbind(bubbles2A, bubbles3A, bubbles4A)
  jump_bubblesA[6]<-unique(SUBS)[z]
  colnames(jump_bubblesA)[6]<-"Strain"
  jump_bubbles<-rbind(jump_bubbles,jump_bubblesA)
  
}
write.csv(jump_bubbles, file=paste(MAIN,"/bubbles.csv", sep=""), row.names = FALSE, quote=FALSE)

#jump_bubbles<-read.csv(file="F:/Giacomo/Dips_paper/Analysis/bubbles.csv", header=TRUE)

pdf("Sigma_bubbles.pdf", height = 5, width = 5)
ggplot(data=jump_bubbles[jump_bubbles$Strain=="MERGED",], aes(x=as.character(Components), y=Sigma, size=100*proportion))+
  geom_point(alpha=0.5, shape=21, color="black", fill="#6477a9")+
  scale_size_area(max_size = 20, name="%")+
  theme_bw()+
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10), limits=c(0.0001,10))+
  xlab("# Subpopulations")+
  ylab(expression(paste("Sigma",sep="")))
dev.off()

pdf("Diffusion_bubbles.pdf", height = 5, width = 5)
ggplot(data=jump_bubbles[jump_bubbles$Strain=="MERGED",], aes(x=as.character(Components), y=Diffusion_Adjusted, size=100*proportion))+
  geom_point(alpha=0.5, shape=21, color="black", fill="#6477a9")+
  scale_size_area(max_size = 20, name="%")+
  theme_bw()+
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10), limits=c(0.0001,10))+
  xlab("# Subpopulations")+
  ylab(expression(paste("Diffusion coefficient D [",µm^2,"/s]",sep="")))
dev.off()

#jump_bubbles$Strain<-factor(jump_bubbles$Strain, levels = c("MERGED","DdipD_pSG_HaloDipD", "DdipD_pSG_HaloDipD_MMC", "DdipABCD_pSG_HaloDipD", "DdipABCD_pSG_HaloDipD_MMC"))

pdf("Comparison_bubbles.pdf", height = 5, width = 7)
ggplot(data=jump_bubbles[jump_bubbles$Components==4,], aes(x=Strain, y=Diffusion_Adjusted, size=100*proportion, fill=Strain)) +
  theme(axis.text.x = element_text(angle=45)) +
  geom_point(alpha=0.5, shape=21, color="black")+
  scale_size_area(max_size = 20, name="%")+
  scale_fill_manual(values=c("darkgrey","#03045e", "#6477a9", "red", "orange"))+
  theme_bw()+
  scale_y_log10(breaks = c(0.0001,0.001, 0.01, 0.1, 1,10), limits=c(0.0001,10))
dev.off()







####from here on is manual analysis and data mining
###
###
##
##
##

DD<-read.csv("F:/Giacomo/Dips_paper/Analysis/20260114_HaloDipD_MMC_IBIDI/FinalOutput.csv", header=TRUE)
AD<-read.csv("F:/Giacomo/Dips_paper/Analysis/20260116_HaloDipD_deltadipA_IBIDI/FinalOutput.csv", header=TRUE)
WT<-read.csv("F:/Giacomo/Dips_paper/Analysis/20260120_WT_MMC_TMR_IBIDI/FinalOutput.csv", header=TRUE)
DD<-DD[!is.na(DD$Repeated),]
AD<-AD[!is.na(AD$Repeated),]
WT<-WT[!is.na(WT$Repeated),]
St<-rbind(DD,AD,WT)

s1<-DD
s2<-AD
s3<-WT
#s4<-fileWI
#s5<-fileWN
#s6<-file1001

length(s3[s3$Unique=="yes",1])

y1<-ggplot(data=s1[s1$seg_pos==5 & s1$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..count../329))+
  geom_histogram(binwidth =10,fill="royalblue4")+
  xlim(0,800)+
  ylim(0,4)+
  theme_bw()

y2<-ggplot(data=s2[s2$seg_pos==5 & s2$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..count../424))+
  geom_histogram(binwidth =10,fill="red")+
  xlim(0,800)+
  ylim(0,4)+
  theme_bw()

y3<-ggplot(data=s3[s3$seg_pos==5 & s3$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..count../383))+
  geom_histogram(binwidth =10,fill="royalblue")+
  xlim(0,800)+
  ylim(0,4)+
  theme_bw()
# 
# y4<-ggplot(data=s4[s4$seg_pos==5 & s4$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..density..))+
#   geom_histogram(binwidth =10,fill="orange")+
#   xlim(0,1500)+
#   ylim(0,0.025)+
#   theme_bw()
# 
# y6<-ggplot(data=s5[s5$seg_pos==5 & s5$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..density..))+
#   geom_histogram(binwidth =10,fill="yellow")+
#   xlim(0,1500)+
#   ylim(0,0.025)+
#   theme_bw()
# 
# y5<-ggplot(data=s6[s6$seg_pos==5 & s6$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..density..))+
#   geom_histogram(binwidth =10,fill="darkgreen")+
#   xlim(0,1500)+
#   ylim(0,0.025)+
#   theme_bw()
(y1/y2/y3)
pdf("Comparison_mjd_weight.pdf", height = 5, width = 8)
(y1/y2/y3)
dev.off()

y1<-ggplot(data=s1[s1$seg_pos==5 & s1$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..count..))+
  geom_histogram(binwidth =10,fill="royalblue4")+
  xlim(0,1500)+
  ylim(0,4000)+
  theme_bw()

y2<-ggplot(data=s2[s2$seg_pos==5 & s2$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..count..))+
  geom_histogram(binwidth =10,fill="orange")+
  xlim(0,1500)+
  ylim(0,4000)+
  theme_bw()

y3<-ggplot(data=s3[s3$seg_pos==5 & s3$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..count..))+
  geom_histogram(binwidth =10,fill="royalblue")+
  xlim(0,1500)+
  ylim(0,4000)+
  theme_bw()

y4<-ggplot(data=s4[s4$seg_pos==5 & s4$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..count..))+
  geom_histogram(binwidth =10,fill="red")+
  xlim(0,1500)+
  ylim(0,4000)+
  theme_bw()

y6<-ggplot(data=s5[s5$seg_pos==5 & s5$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..count..))+
  geom_histogram(binwidth =10,fill="yellow")+
  xlim(0,1500)+
  ylim(0,4000)+
  theme_bw()

y5<-ggplot(data=s6[s6$seg_pos==5 & s6$seg.motion=="diffusion",], aes(x=seg.mjd, weight=seg.mjd_n, y=..count..))+
  geom_histogram(binwidth =10,fill="darkgreen")+
  xlim(0,1500)+
  ylim(0,4000)+
  theme_bw()
(y1/y2)|(y3/y4)|(y5/y6)
pdf("Comparison_mjd_weight_count.pdf", height = 5, width = 8)
(y1/y2)|(y3/y4)|(y5/y6)
dev.off()