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
IRI <- RenameIdents(object = IRI, 'X4hours' = "4hours",  'X12hours' = "12hours")
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
IRI2 <- RenameIdents(object = IRI2, 'X4hours' = "4hours",  'X12hours' = "12hours")

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
IRI3 <- RenameIdents(object = IRI3, 'X4hours' = "4hours",'X12hours' = "12hours")

levels(x = IRI3) <- c("ctrl","4hours","12hours")

marker = c("Slc1a5","Slc5a8","Slc5a12")

DotPlot(IRI, features = rev(marker),cols = c("blue", "red"), dot.scale = 15) +ggtitle("PTS1")
DotPlot(IRI2, features = rev(marker), cols = c("blue", "red"), dot.scale = 15) + ggtitle("PTS2")
DotPlot(IRI3, features = rev(marker), cols = c("blue", "red"), dot.scale = 15) + ggtitle("PTS3")
