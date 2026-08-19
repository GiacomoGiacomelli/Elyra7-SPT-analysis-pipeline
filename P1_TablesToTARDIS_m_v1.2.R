library("spatstat")
library(ggplot2)
library(stringr)
library(fields)  
library(tidyverse)
library(conflicted)
library(patchwork)
library(RANN)
library(ggbeeswarm)

library("dbscan")
library("mixtools")
library("MASS")
library("agricolae")

####clockwise vs counterclockwise

#' @title Check whether points for an owin are clockwise
#' @param x a dataframe with x coordinates in the first column and y coordinates in the second. 
#' @details Similarly to owin, the polygon should not be closed
#' @return A logical telling whether the polygon is arranged clockwise.
#' @author The idea has been scavenged from https://stackoverflow.com/a/1165943/1082004

clockwise <- function(x) {
  
  x.coords <- c(x[[1]], x[[1]][1])
  y.coords <- c(x[[2]], x[[2]][1])
  
  double.area <- sum(sapply(2:length(x.coords), function(i) {
    (x.coords[i] - x.coords[i-1])*(y.coords[i] + y.coords[i-1])
  }))
  
  double.area > 0
}

polyarea <- function(x, y) {
  0.5 * abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y))
}

###PART 1######
#INPUTS
NPmax<-60000  #Maximum number of photons
NPmin<-400     #Minimum number of photons
PSFmax<-30000 ##Maximum PSF half width (nm)
PSFmin<-100   ##Minimum PSF half width (nm)
FF<-0      ##Frames to be removed due to bleaching phase
folder<-"Z:/Paul Molis/Experimente/P1_Elongasome_Divisom_Taggen/E7_SPT_Elongsom_tracking_Parameters_testing"  ##Main folder, it contains all the separate strains, in their respective folders, that need to be analyzed
Av<-"outline" ##use "outline" if you want to use the cell outline for the rotation of the cells, use "points" if you want to use the points obtained during SMLM        
grouped<-"no" ##will change the y flipping behaviour if using the Elyra grouping function
elyradrift<-"yes" ##will use the "drift.txt" file to drift the data
autodrift<-"nyet"  ##will activate automated lateral shifting if active
FShift<-5000 ##Number of frames to be checked to apply shift (the more the better, the more the slower)
Tab_pos<-"/tables"
##Membrane<-"yes"  ###if yes, it will attempt to have areas with only membrane signal

setwd(folder)  #Choose the folder containing all strains folders
MAIN<-getwd()
SUBS<-list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
#z<-SUBS[1]
for (z in SUBS) {
  tablenames <- Sys.glob(paste(MAIN,"/",z,Tab_pos,"/series*.txt",sep=""))
  foldernames <- Sys.glob(paste(MAIN,"/",z,"/rois/cells*",sep=""))
  foldernames1 <- Sys.glob(paste(MAIN,"/",z,"/rois/smallcells*",sep=""))
  dir.create(paste(MAIN,"/",z,"/areas/",sep=""))
  F_table<-data.frame()
  F_TARDIS<-data.frame()
  
  if (elyradrift=="yes") {
    a<-read.table(paste(z,"/drift.txt",sep=""), sep="\t",header=TRUE)
    files <- list.files(
      path = paste(MAIN,"/",z,Tab_pos, sep=""),
      pattern = "^series[0-9]{3}\\.txt$"
    )
  }
  #n<-1
  for (n in 1:length(tablenames)) {
    Table<-tablenames[n]
    print(Table)
    CellsOutlineFolder<-paste(foldernames[n],"/",sep="")
    SmallCellsOutlineFolder<-paste(foldernames1[n],"/",sep="")#####
    print(CellsOutlineFolder)
    
    #####a<-read.table(paste(z,"/drift.txt",sep=""), sep="\t",header=TRUE)
    
    TTemp<-read.table(Table, header = TRUE,sep="\t")
    T <- na.omit(TTemp)
    if (length(TTemp$Index)-length(T$Index)>0) {
      ccc<-length(TTemp$Index)-length(T$Index)
      print(paste(ccc," rows contained NA and needed to be removed", sep = ""))
    }
    
    if (grouped=="yes"){
    T$Position.Y..nm.<-61940.6-T[,6]
    }
    
    #test parameters for filtering
    # ggplot(data=T[T$Number.Photons>=0 & T$PSF.half.width..nm.>=70 & T$PSF.half.width..nm.<=1000 & T$Number.Photons<=2000,], aes(x=PSF.half.width..nm.))+
    #   geom_histogram()
    # 
    # ggplot(data=T[T$Number.Photons>=0 & T$PSF.half.width..nm.>=70 & T$PSF.half.width..nm.<=1000 & T$Number.Photons<=2000,], aes(x=Number.Photons))+
    #   geom_histogram()
    #}
    
    ##determine shift based on filtered data removing data, this is also where you decide how to filter your data
    T<-T[T$Number.Photons<NPmax & T$PSF.half.width..nm.>PSFmin & T$Number.Photons>NPmin & T$PSF.half.width..nm.<PSFmax & T$First.Frame>FF,]
    #plot(T[T$First.Frame<2000,]$Position.Y..nm.~T[T$First.Frame<2000,]$Position.X..nm.)
    if (elyradrift=="yes") {
      T$Position.X..nm.<-T$Position.X..nm.+a[a$file==files[n],]$x
      T$Position.Y..nm.<-61940.6-(T$Position.Y..nm.+a[a$file==files[n],]$y)
    } else {
      T$Position.Y..nm.<-61940.6-(T$Position.Y..nm.+a[a$file==files[n],]$y)
    }
    #plot(T[T$First.Frame<2000,]$Position.Y..nm.~T[T$First.Frame<2000,]$Position.X..nm.)
    
    outlines <- Sys.glob(paste(CellsOutlineFolder,"Cell*",sep=""))
    smalloutlines <- Sys.glob(paste(SmallCellsOutlineFolder,"Cell*",sep="")) ####
    
    listC<-c()
    listS<-c()
    
    for (i in 1:length(smalloutlines)){
      cell<-read.table(outlines[i], header=TRUE) 
      smallcell<-read.table(smalloutlines[i], header=TRUE) 
      
      if (clockwise(data.frame(x=cell$X, y=cell$Y))) {
        cellV<-list(list(x=rev(cell$X*1000), y=rev(cell$Y*1000)))
        cellVS<-list(list(x=smallcell$X*1000, y=smallcell$Y*1000))
        listC<-c(listC,cellV)
        listS<-c(listS,cellVS)
      } else {
        cellV<-list(list(x=cell$X*1000, y=cell$Y*1000))
        cellVS<-list(list(x=rev(smallcell$X*1000), y=rev(smallcell$Y*1000)))
        listC<-c(listC,cellV)
        listS<-c(listS,cellVS)
      }

      W1 <- owin(poly = c(cellV, cellVS))
      #plot(W1)
      
    }
      # Fixed window
    W <- owin(poly=c(listC, listS))
    plot(W)
    
    #ggplot(data=T[T$Number.Photons>=0 & T$PSF.half.width..nm.>=70 & T$PSF.half.width..nm.<=1000 & T$Number.Photons<=2000,], aes(x=Position.X..nm., y=Position.Y..nm.))+
    #  geom_point()+
    #  coord_fixed()
    
    x<-T[c(1:FShift),5]
    y<-T[c(1:FShift),6]

    # Points outside window
    X <- ppp(x,y,owin(xrange=c(0,61940.6), yrange=c(0,61940.6)))
    best_dx <- 0
    best_dy <- 0
    
    if (autodrift=="yes") {
    # Plot original points
    plot(W, main = "Original Points and Window")
    plot(X, add = TRUE, col = "red")

    # Define a grid of shifts
    dx_vals <- seq(-500, 500, by = 25)
    dy_vals <- seq(-500, 500, by = 25)

    # Store best result
    max_inside0 <- sum(inside.owin(X$x, X$y, W))
    max_inside <- sum(inside.owin(X$x, X$y, W))

    # Try all shift combinations
    for (dx in dx_vals) {
      for (dy in dy_vals) {
        X_shifted <- shift(X, vec = c(dx, dy))
        count <- sum(inside.owin(X_shifted$x, X_shifted$y, W))
    
        if (count > max_inside) {
        max_inside <- count
        best_dx <- dx
        best_dy <- dy
        }
      }
    }

    # Apply best shift
    X_best <- shift(X, vec = c(best_dx, best_dy))
    # Results
    cat("Best shift:\n")
    cat("dx =", best_dx, ", dy =", best_dy, "\n")
    cat("Points inside window before:", max_inside0, "\n")
    cat("Points inside window after:", max_inside, "\n")
    
    # Plot result
    plot(W, main = "Shifted Points in Window")
    plot(X_best, add = TRUE, col = "blue")
    
    T$Position.X..nm.<-T$Position.X..nm.+best_dx
    T$Position.Y..nm.<-T$Position.Y..nm.+best_dy
    
    } else {
    
    plot(W, main = "Original Points and Window")
    plot(X, add = TRUE, col = "red")
    
    }

    LISTofAREAS<-c()
    LISTofAREAS0<-c()
    LISTofAREASF<-c()
    
    LISTofAREAS_m<-c()
    LISTofAREAS0_m<-c()
    LISTofAREASF_m<-c()
    
    T <- T %>%
      mutate(FOV = n,
             CellName = NA,
             CellLength_Outline = NA, 
             CellWidth_Outline = NA, 
             CellArea = NA,
             Position.X.R = NA,
             Position.Y.R = NA,
             Position.X.0 = NA,
             Position.X.polar = NA,
             Position.X.betwixt = NA,
             Position.X.midcell = NA,
             Position.Y.0 = NA,
             Unique = "no",
             CellLength_SMLM = NA,
             CellWidth_SMLM = NA,
             Position.X.N = NA,
             Position.Y.N = NA,
             Midcell_portion = 0,
             Pole_M_portion = 0,
             Pole_m_portion = 0,
             Polar_portion = 0,
             Poles_Mm_ratio = NA,
             Midcell_portion_m = 0,
             Pole_M_portion_m = 0,
             Pole_m_portion_m = 0,
             Polar_portion_m = 0,
             Poles_Mm_ratio_m = NA,
             Strain = z,   #####basename(getwd())    check if that is correct
             MaxProteinPosition = NA,
             MaxProteinPosition_N = NA,
             MaxPP_DistFromMid = NA,
             CellArea_SMLM = NA,
             DistFromMid = NA,
             Length_Width_Ratio_Outline = NA,
             Length_Width_Ratio_SMLM = NA,
             Circularity = NA,
             Perimeter = NA,
             Protein_Position = NA,
             DX = best_dx,
             DY = best_dy,
             ProtCount = NA,
             NearestNeighbour = NA,
             NearestIndex = NA,
             MembraneAssociated = "no",
             ProtCount_m = 0
             )
    
    
    for (o in 1:length(outlines)){
      cell<-read.table(outlines[o], header=TRUE) 
      smallcell<-read.table(smalloutlines[o], header=TRUE) 
      
      if (clockwise(data.frame(x=cell$X, y=cell$Y))) {
        cellV<-list(list(x=rev(cell$X*1000), y=rev(cell$Y*1000)))
        cellVS<-list(list(x=smallcell$X*1000, y=smallcell$Y*1000))
        listC<-c(listC,cellV)
        listS<-c(listS,cellVS)
      } else {
        cellV<-list(list(x=cell$X*1000, y=cell$Y*1000))
        cellVS<-list(list(x=rev(smallcell$X*1000), y=rev(smallcell$Y*1000)))
        listC<-c(listC,cellV)
        listS<-c(listS,cellVS)
      }
      
      cm1 <- owin(poly = c(cellV))
      pp3<-as.ppp(T[c(5:6)],cm1)
      cm1_m <- owin(poly = c(cellV, cellVS))
      pp3_m<-as.ppp(T[c(5:6)],cm1_m)
      realsA<-subset(pp3,cm1)
      plot(realsA)
      realsA_m<-subset(pp3_m,cm1_m)
      plot(realsA_m)
    
      if (length(realsA$x)>1){
        realsA_DF<-as.data.frame(realsA)
        T[which(paste(T$Position.X..nm.,T$Position.Y..nm.) %in% paste(realsA_DF$x,realsA_DF$y)),]$CellName<-paste(n,"_",str_sub(outlines[o], -11,-5),sep="")
        T[T$Position.X..nm.==realsA[1]$x & T$Position.Y..nm.==realsA[1]$y,]$CellLength_Outline<-diameter(cm1)
        T[T$Position.X..nm.==realsA[1]$x & T$Position.Y..nm.==realsA[1]$y,]$CellArea<-area.owin(cm1)
        LISTofAREAS<-c(LISTofAREAS,cellV)
        
        if (length(realsA_m$x)>1){
          realsA_mDF<-as.data.frame(realsA_m)
          T[which(paste(T$Position.X..nm.,T$Position.Y..nm.) %in% paste(realsA_mDF$x,realsA_mDF$y)),]$MembraneAssociated<-"yes"
          LISTofAREAS_m<-c(LISTofAREAS_m,c(cellV,cellVS))
        }
        
        CC<-T[T$CellName==paste(n,"_",str_sub(outlines[o], -11,-5),sep="") & !is.na(T$CellName),]
        CC$CellLength_Outline<-CC[1,]$CellLength_Outline
        CC$CellArea<-CC[1,]$CellArea
        CC[1,]$Unique<-"yes"
        
        if (Av=="outline") {
        c1<-cell  #by using "cell", we are deciding to align based on the outline, by using "CC[,c(5,6)]" we are using the points collected during SMLM 
        } else {
        c1<-CC[,c(5,6)]  
        }
        
        plot(c1)
        plot(CC$Position.Y..nm.~CC$Position.X..nm.)
        
        D <- rdist(c1)
        max_idx <- which(D == max(D), arr.ind = TRUE)[1, ]
        M <- rbind(c1[max_idx[1], ], c1[max_idx[2], ])
        #plot data
        plot(M[,1],M[,2])
        #calculate rotation angle
        alpha <- -atan((M[1,2]-tail(M,1)[,2])/(M[1,1]-tail(M,1)[,1]))
        #rotation matrix
        rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
        #shift, rotate, shift back
        M2 <- t(rotm %*% (t(c1)-c(M[1,1],M[1,2]))+c(M[1,1],M[1,2]))
        #plot(M2)
      
        M3 <- t(rotm %*% (t(CC[c(5,6)])-c(M[1,1],M[1,2])*1000)+c(M[1,1],M[1,2])*1000)
        #plot
        #plot(M3)
        
        cell0<-list(list(x=rev(M2[,1]-min(M2[,1])), y=rev(M2[,2]-min(M2[,2]))))
        LISTofAREAS0<-c(LISTofAREAS0,cell0)
        
        cellF<-list(list(x=rev(M2[,1]-min(M2[,1])+ (as.numeric(str_sub(outlines[o], -7,-5))-1)*10), y=rev(M2[,2]-min(M2[,2]) + (n-1)*10)))
        LISTofAREASF<-c(LISTofAREASF,cellF)
        
        CC$Position.X.R <- M3[,1]
        CC$Position.Y.R <- M3[,2]
        
        mmmX <- min(M2[,1])*1000
        mmmY <- min(M2[,2])*1000
        #MMMY <- max(M2[,2])*1000-min(M2[,2])*1000
        
        CC <- CC %>%
          mutate(Position.X.0 = Position.X.R - mmmX,
                 Position.Y.0 = Position.Y.R - mmmY,
                 )
        
        CC$Perimeter<-perimeter(cm1)
        
        # w1<--1
        # w2<--1
        # w1<-(-CC$CellLength_Outline[1]-sqrt((CC$CellLength_Outline[1]^2)+4*CC$CellArea[1]*((pi/4)-1)))/(2*((pi/4)-1))
        # w2<-(-CC$CellLength_Outline[1]+sqrt((CC$CellLength_Outline[1]^2)+4*CC$CellArea[1]*((pi/4)-1)))/(2*((pi/4)-1))
        
        nn <- RANN::nn2(CC[,c(5,6)], k = 2)
        CC$NearestNeighbour <- nn$nn.dists[, 2]
        CC$NearestIndex <- nn$nn.idx[, 2]
        
        #CC$CellWidth_Outline<-ifelse(w1 > 0 & w1 < CC$CellLength_Outline[1], w1, w2)
        CC$CellWidth_Outline<-max(M2[,2])*1000-min(M2[,2])*1000
        
        CC$CellLength_SMLM<-max(CC$Position.X.0)
        CC$CellWidth_SMLM<-max(CC$Position.Y.0)
        CC$Position.X.N<-(CC$Position.X.0 / CC$CellLength_Outline[1])*2000
        CC$Position.Y.N<-(CC$Position.Y.0 / CC$CellWidth_Outline[1])*800
        
        CC <- CC %>%
          mutate(Position.X..nm. = Position.X..nm. - min(Position.X..nm.) + (as.numeric(str_sub(outlines[o], -7,-5))-1)*10000.0,
                 Position.Y..nm. = Position.Y..nm. - min(Position.Y..nm.) + (n-1)*10000.0,
                 DistFromMid = (Position.X.0 / CellLength_Outline)-0.5)

        CC$ProtCount<-length(CC$Position.X.0)
        CC$ProtCount_m<-length(CC[CC$MembraneAssociated=="yes",]$Position.X.0)
        
        if(length(CC[CC$Position.X.0>=0.4*max(CC$CellLength_Outline) & CC$Position.X.0<=0.6*max(CC$CellLength_Outline),]$Protein_Position)>0){
          CC[CC$Position.X.0>=0.4*max(CC$CellLength_Outline) & CC$Position.X.0<=0.6*max(CC$CellLength_Outline),]$Protein_Position<-"Midcell"
        }
        
        if(length(CC[CC$Position.X.0>0.15*max(CC$CellLength_Outline) & CC$Position.X.0<0.4*max(CC$CellLength_Outline),]$Protein_Position)>0){
          CC[CC$Position.X.0>0.15*max(CC$CellLength_Outline) & CC$Position.X.0<0.4*max(CC$CellLength_Outline),]$Protein_Position<-"Betwixt"
        }
        
        if(length(CC[CC$Position.X.0>0.6*max(CC$CellLength_Outline) & CC$Position.X.0<0.85*max(CC$CellLength_Outline),]$Protein_Position)>0){
          CC[CC$Position.X.0>0.6*max(CC$CellLength_Outline) & CC$Position.X.0<0.85*max(CC$CellLength_Outline),]$Protein_Position<-"Betwixt"
        }
        
        CC$Midcell_portion<-length(CC[CC$Protein_Position=="Midcell",]$Index)/CC$ProtCount[1]
        CC$Midcell_portion_m<-length(CC[CC$Protein_Position=="Midcell" & CC$MembraneAssociated=="yes",]$Index)/CC$ProtCount_m[1]  
        
        CC[CC$Position.X.0>=0.4*max(CC$CellLength_Outline) & CC$Position.X.0<=0.5*max(CC$CellLength_Outline),]$Position.X.midcell<-CC[CC$Position.X.0>=0.4*max(CC$CellLength_Outline) & CC$Position.X.0<=0.5*max(CC$CellLength_Outline),]$Position.X.N
        CC[CC$Position.X.0>0.5*max(CC$CellLength_Outline) & CC$Position.X.0<=0.6*max(CC$CellLength_Outline),]$Position.X.midcell<-2000-CC[CC$Position.X.0>0.5*max(CC$CellLength_Outline) & CC$Position.X.0<=0.6*max(CC$CellLength_Outline),]$Position.X.N
        
        CC[CC$Position.X.0>0.15*max(CC$CellLength_Outline) & CC$Position.X.0<0.4*max(CC$CellLength_Outline),]$Position.X.betwixt<-CC[CC$Position.X.0>0.15*max(CC$CellLength_Outline) & CC$Position.X.0<0.4*max(CC$CellLength_Outline),]$Position.X.N
        CC[CC$Position.X.0>0.6*max(CC$CellLength_Outline) & CC$Position.X.0<0.85*max(CC$CellLength_Outline),]$Position.X.betwixt<-2000-CC[CC$Position.X.0>0.6*max(CC$CellLength_Outline) & CC$Position.X.0<0.85*max(CC$CellLength_Outline),]$Position.X.N
        
        CC[CC$Position.X.0<=0.15*max(CC$CellLength_Outline),]$Position.X.polar<-CC[CC$Position.X.0<=0.15*max(CC$CellLength_Outline),]$Position.X.N
        CC[CC$Position.X.0>=0.85*max(CC$CellLength_Outline),]$Position.X.polar<-2000-CC[CC$Position.X.0>=0.85*max(CC$CellLength_Outline),]$Position.X.N
        
          
        if (length(CC[CC$Position.X.0>=0.85*max(CC$CellLength_Outline),]$Index) >= length(CC[CC$Position.X.0<=0.15*max(CC$CellLength_Outline),]$Index) & length(CC[CC$Position.X.0>=0.85*max(CC$CellLength_Outline),]$Index)>0){
          CC[CC$Position.X.0>=0.85*max(CC$CellLength_Outline),]$Protein_Position<-"Major_pole"
          if (length(CC[CC$Position.X.0<=0.15*max(CC$CellLength_Outline),]$Protein_Position)>0){
          CC[CC$Position.X.0<=0.15*max(CC$CellLength_Outline),]$Protein_Position<-"Minor_pole"
          }
        }
        if (length(CC[CC$Position.X.0>=0.85*max(CC$CellLength_Outline),]$Index) < length(CC[CC$Position.X.0<=0.15*max(CC$CellLength_Outline),]$Index) & length(CC[CC$Position.X.0<=0.15*max(CC$CellLength_Outline),]$Index)>0){
          if (length(CC[CC$Position.X.0>=0.85*max(CC$CellLength_Outline),]$Protein_Position)>0){
          CC[CC$Position.X.0>=0.85*max(CC$CellLength_Outline),]$Protein_Position<-"Minor_pole"
          }
          CC[CC$Position.X.0<=0.15*max(CC$CellLength_Outline),]$Protein_Position<-"Major_pole"
        }
        
        CC$Polar_portion<-length(CC[CC$Protein_Position=="Major_pole" | CC$Protein_Position=="Minor_pole",]$Index)/CC$ProtCount[1]
        CC$Polar_portion_m<-length(CC[(CC$Protein_Position=="Major_pole" | CC$Protein_Position=="Minor_pole") & CC$MembraneAssociated=="yes",]$Index)/CC$ProtCount_m[1]
        
        CC$Pole_M_portion<-length(CC[CC$Protein_Position=="Major_pole",]$Index)/CC$ProtCount[1]
        CC$Pole_M_portion_m<-length(CC[CC$Protein_Position=="Major_pole" & CC$MembraneAssociated=="yes",]$Index)/CC$ProtCount_m[1]
        
        CC$Pole_m_portion<-length(CC[CC$Protein_Position=="Minor_pole",]$Index)/CC$ProtCount[1]
        CC$Pole_m_portion_m<-length(CC[CC$Protein_Position=="Minor_pole" & CC$MembraneAssociated=="yes",]$Index)/CC$ProtCount_m[1]
        
        CC$Poles_Mm_ratio<- CC$Pole_m_portion[1]/CC$Pole_M_portion[1]
        CC$Poles_Mm_ratio_m<- CC$Pole_m_portion_m[1]/CC$Pole_M_portion_m[1]
        
 
        H <- hist(CC$Position.X.0, plot = FALSE, breaks=seq(0,CC$CellLength_Outline[1],l=round((CC$CellLength_Outline[1]/25)+1, digits=0)))
        CC$MaxProteinPosition<-H$breaks[which.max(H$counts)]   #not sure if useful
        CC$MaxPP_DistFromMid<-abs((H$breaks[which.max(H$counts)]/max(CC$CellLength_Outline))-0.5)
 
        HN <- hist(CC$Position.X.N, plot = TRUE, breaks=seq(0,2000,l=81))
        CC$MaxProteinPosition_N<-HN$breaks[which.max(HN$counts)]  #not sure if useful
      
        points <- cbind(CC$Position.X..nm. , CC$Position.Y..nm.)
        hull_indices <- chull(points)
        hull_points <- points[c(hull_indices, hull_indices[1]), ]
        
        CC$CellArea_SMLM<-abs(polyarea(hull_points[,1], hull_points[,2]))
        CC$Length_Width_Ratio_SMLM<-CC$CellLength_SMLM[1]/CC$CellWidth_SMLM[1]   ####need to think if there is a better alternative
        CC$Length_Width_Ratio_Outline<-CC$CellLength_Outline[1]/CC$CellWidth_Outline[1]   ####need to think if there is a better alternative
        CC$Circularity<-(4* pi *CC$CellArea[1]) / (CC$Perimeter[1]^2)
        
        F_table<-rbind(F_table,CC)
      }
    }
  saveRDS(LISTofAREAS, file=paste(MAIN,"/",z,"/areas/output_area_",n,".RData", sep=""))
  saveRDS(LISTofAREAS0, file=paste(MAIN,"/",z,"/areas/output_area_0_",n,".RData", sep=""))
  saveRDS(LISTofAREASF, file=paste(MAIN,"/",z,"/areas/output_area_F_",n,".RData", sep=""))
  }
colnames(F_table)[c(1,2,5,6,11,8,12,9,10,7)]<-c("id","frame","x [nm]","y [nm]","sigma [nm]","intensity [photon]","offset [photon]","bkgstd [photon]","chi2","uncertainty_xy [nm]")
write.csv(F_table, file=paste(MAIN,"/",z,"/output.csv", sep=""), row.names = FALSE, quote=FALSE)
F_TARDIS<-F_table[,c(1,2,5,6,11,8,12,9,10,7)]
write.csv(F_TARDIS, file=paste(MAIN,"/",z,"/output_TARDIS.csv", sep=""), row.names = FALSE, quote=FALSE)
F_TARDIS_test<-F_table[c(1:100000),c(1,2,5,6,11,8,12,9,10,7)]
write.csv(F_TARDIS, file=paste(MAIN,"/",z,"/output_TARDIS.csv", sep=""), row.names = FALSE, quote=FALSE)
write.csv(F_TARDIS_test, file=paste(MAIN,"/",z,"/output_TARDIS_test.csv", sep=""), row.names = FALSE, quote=FALSE)

}
mean(F_TARDIS$`uncertainty_xy [nm]`)

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


warnings()
