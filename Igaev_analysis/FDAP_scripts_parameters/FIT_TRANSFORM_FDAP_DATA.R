#!/usr/bin/env Rscript
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
library(ggplot2)
require(minpack.lm)

filenames <- Sys.glob(file.path("*.txt"))
filenames
n <- length(filenames)
timelist = list()
datalist = list()
dir.create("jpeg_fits")
dir.create("dat")
for (i in filenames){  
  i
  x <- read.table(i, sep = "\t", header = TRUE, comment.char="")
  name <- gsub(".txt","",i)
  x$Name <- name
  x <- cbind(x[,c("Name","Time..s.","ND.T")], x[grep("X" , names(x))])
  x$ND.T <- x$ND.T - 1
  #x <- x[c("Name","Time..s.","ND.T","X")]
  #cols.num <- c("Time..s.","ND.T","X")
  #x[cols.num] <- sapply(x[cols.num],as.numeric)
  
  tim <- data.frame(matrix(ncol=0,nrow=113))
  tim$Time <- x$Time..s.
  names(tim)[names(tim) == 'Time'] <- name
  timelist[[name]] <- tim
  
  Int_columns <- colnames(x[grep("X", names(x))])
  Int_columns
  dat <- data.frame(matrix(ncol=0,nrow=113))
  for (j in Int_columns) {
    Int0_bg <- x[j][,1][1]
    x[j] <- x[j] - Int0_bg
    time <- x$ND.T
    y <- x[j]
    colnames(y) <- "y"
    str <- cbind(time,y)[-c(1),]
    data_to_fit <- structure(list(x = str$time, y=str$y), class = "data.frame", row.names = c(NA, -112L), .Names = c("x", "y"))
    
    fit <- nlsLM(y ~ y0 + a1*exp(-x/b1) + a2*exp(-x/b2), data=data_to_fit, start = list(y0=10, a1=1000, b1=1, a2=10, b2=5), control = list(maxiter = 100))
    print(coef(fit))
    Int0_extrap <- coef(fit)[1] + coef(fit)[2] + coef(fit)[4]
    x[j][,1][1] <- Int0_extrap
    
    jpeg(paste("jpeg_fits/",paste(paste(name,j,sep="_"),"jpg",sep="."),sep=""))
    plot(x[j][,1]~Time..s., data = x, main =paste(name,j,sep="_"), xlab = "Time (s)", ylab = "Intensity",ylim=c(0,max(x[,4])))
    curve(predict(fit, newdata = data.frame(x)), col = "red", add = TRUE)
    dev.off()
    
    x[j][,1] <- x[j][,1]/Int0_extrap
    
    dat[name] <- x[j]
    print(colnames(dat))
    write.table(dat, paste("dat/",paste(paste(name,j,sep="_"),"dat",sep="."),sep=""), row.names = FALSE, col.names = FALSE, fileEncoding = "ASCII")
    datalist[paste(name,j,sep="_")] <- dat
    }
  }

big_data = do.call(cbind, datalist)
big_time = do.call(cbind, timelist)
big_time <- big_time[-c(1),]
write.table(big_data, "test_final.csv", sep = ";", row.names = FALSE) #_time
write.table(big_time, "test_time_final.csv", sep = ";", row.names = FALSE)

