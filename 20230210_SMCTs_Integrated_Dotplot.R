
library(Seurat)
library(ggplot2)
library(dplyr)

PTS1_ctrl <- read.table("/path/to/SpeedTest_ctrl_PTS1.txt")
PTS1_4h <- read.table("/path/to/SpeedTest_4hours_PTS1.txt")
PTS1_12h <- read.table("/path/to/SpeedTest_12hours_PTS1.txt")

I1_ctrl <- CreateSeuratObject(PTS1_ctrl)
I1_4h <- CreateSeuratObject(PTS1_4h)
I1_12h <- CreateSeuratObject(PTS1_12h)

IRI <- merge(I1_ctrl, c(I1_4h,I1_12h), add.cell.ids = c("ctrl","4h","12h"))
IRI <- NormalizeData(IRI)
IRI <- ScaleData(IRI,features = c("Slc1a5","Slc5a8","Slc5a12"))
IRI <- RenameIdents(object = IRI,
                     'IRIsham1b1' = 'ctrl', 
                     'IRIsham1b2' = 'ctrl', 
                     'IRIsham2' = 'ctrl',
                     "IRIsham3" = 'ctrl',
                     "IRI12h1b1" = "12hours", 
                     "IRI12h1b2" = "12hours",
                     "IRI12h2" = "12hours",
                     "IRI12h3" = "12hours",
                     "IRI4h1" = "4hours",
                     "IRI4h2" = "4hours",
                     "IRI4h3" = "4hours")
levels(x = IRI) <- c("ctrl","4hours","12hours")

PTS2_ctrl <- read.table("/path/to/SpeedTest_ctrl_PTS2.txt")
PTS2_4h <- read.table("/path/to/SpeedTest_4hours_PTS2.txt")
PTS2_12h <- read.table("/path/to/SpeedTest_12hours_PTS2.txt")

I2_ctrl <- CreateSeuratObject(PTS2_ctrl)
I2_4h <- CreateSeuratObject(PTS2_4h)
I2_12h <- CreateSeuratObject(PTS2_12h)

IRI2 <- merge(I2_ctrl, c(I2_4h,I2_12h), add.cell.ids = c("ctrl","4h","12h"))
IRI2<- NormalizeData(IRI2)
IRI2<- ScaleData(IRI2,features = c("Slc1a5","Slc5a8","Slc5a12"))
IRI2 <- RenameIdents(object = IRI2,
                     'IRIsham1b1' = 'ctrl', 
                     'IRIsham1b2' = 'ctrl', 
                     'IRIsham2' = 'ctrl',
                     "IRIsham3" = 'ctrl',
                     "IRI12h1b1" = "12hours", 
                     "IRI12h1b2" = "12hours",
                     "IRI12h2" = "12hours",
                     "IRI12h3" = "12hours",
                     "IRI4h1" = "4hours",
                     "IRI4h2" = "4hours",
                     "IRI4h3" = "4hours")
levels(x = IRI2) <- c("ctrl","4hours","12hours")


PTS3_ctrl <- read.table("/path/to/SpeedTest_ctrl_PTS3.txt")
PTS3_4h <- read.table("/path/to/SpeedTest_4hours_PTS3.txt")
PTS3_12h <- read.table("/path/to/SpeedTest_12hours_PTS3.txt")

I3_ctrl <- CreateSeuratObject(PTS3_ctrl)
I3_4h <- CreateSeuratObject(PTS3_4h)
I3_12h <- CreateSeuratObject(PTS3_12h)

IRI3 <- merge(I3_ctrl, c(I3_4h,I3_12h), add.cell.ids = c("ctrl","4h","12h"))
IRI3 <- NormalizeData(IRI3)
IRI3 <- ScaleData(IRI3, features = c("Slc1a5","Slc5a8","Slc5a12"))
IRI3 <- RenameIdents(object = IRI3,
                     'IRIsham1b1' = 'ctrl', 
                     'IRIsham1b2' = 'ctrl', 
                     'IRIsham2' = 'ctrl',
                     "IRIsham3" = 'ctrl',
                     "IRI12h1b1" = "12hours", 
                     "IRI12h1b2" = "12hours",
                     "IRI12h2" = "12hours",
                     "IRI12h3" = "12hours",
                     "IRI4h1" = "4hours",
                     "IRI4h2" = "4hours",
                     "IRI4h3" = "4hours")

levels(x = IRI3) <- c("ctrl","4hours","12hours")

marker = c("Slc5a8","Slc5a12","Slc1a5")

p <- DotPlot(IRI, features = marker,cols = c("blue", "red")) + scale_radius(range=c(1,15),breaks=c(25,50,75)) +scale_color_gradient(low="blue", high="red",limits = c(-1.25, 1.25)) + ggtitle("PTS1")
q <- DotPlot(IRI2, features = marker, cols = c("blue", "red"))+ scale_radius(range=c(1,15),breaks=c(5,15,25)) +scale_color_gradient(low="blue", high="red",limits = c(-1.25, 1.25)) + ggtitle("PTS2")
r <- DotPlot(IRI3, features = marker, cols = c("blue", "red"))+ scale_radius(range=c(1,15),breaks=c(10,30,50)) +scale_color_gradient(low="blue", high="red",limits = c(-1.25, 1.25)) + ggtitle("PTS3")

marker = c("Slc1a5")

s <- DotPlot(IRI, features = rev(marker),cols = c("blue", "red")) + scale_radius(limits = c(0, 0.9),range=c(1,15),breaks=c(0.2,0.4,0.6,0.8)) +scale_color_gradient(low="blue", high="red",limits = c(-1.25, 1.25)) +ggtitle("PTS1")
t <- DotPlot(IRI2, features = rev(marker), cols = c("blue", "red"))+ scale_radius(limits = c(0, 0.9),range=c(1,15),breaks=c(0.2,0.4,0.6,0.8)) +scale_color_gradient(low="blue", high="red",limits = c(-1.25, 1.25)) + ggtitle("PTS2")
u <- DotPlot(IRI3, features = rev(marker), cols = c("blue", "red"))+ scale_radius(limits = c(0, 0.9),range=c(1,15),breaks=c(0.2,0.4,0.6,0.8)) +scale_color_gradient(low="blue", high="red",limits = c(-1.25, 1.25)) + ggtitle("PTS3")
