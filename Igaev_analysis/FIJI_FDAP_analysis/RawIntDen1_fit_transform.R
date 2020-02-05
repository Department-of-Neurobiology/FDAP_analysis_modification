setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
library(ggplot2)
require(minpack.lm)
filenames <- Sys.glob(file.path("*RawIntDen1.csv"))  #parameters
filenames
datalist = list()
dir.create("jpeg_fits")
dir.create("dat")
for (i in filenames){  
  x <- read.table(i, sep = ",", header = TRUE, fill = TRUE)
  x
  name <- gsub("_RawIntDen1.csv","",i)
  Int0_bg <- x$RawIntDen1[1]
  x$RawIntDen1 <- x$RawIntDen1 - Int0_bg
  x$X <- x$X - 1
  time <- x$X
  y <- x$RawIntDen1
  #colnames(y) <- "y"
  str <- cbind(time,y)[-c(1),]
  str <- as.data.frame(str)
  data_to_fit <- structure(list(x = str$time, y=str$y), class = "data.frame", row.names = c(NA, -112L), .Names = c("x", "y"))
  
  fit <- nlsLM(y ~ y0 + a1*exp(-x/b1) + a2*exp(-x/b2), data=data_to_fit, start = list(y0=10, a1=1000, b1=1, a2=10, b2=5), control = list(maxiter = 100))
  print(coef(fit))
  Int0_extrap <- coef(fit)[1] + coef(fit)[2] + coef(fit)[4]
  
  x$RawIntDen1[1] <- Int0_extrap
  jpeg(paste("jpeg_fits/",paste(paste(name,sep="_"),"jpg",sep="."),sep=""))
  plot(RawIntDen1~X, data = x, main =paste(name,sep="_"), xlab = "Intensity", ylab = "Time (s)")
  curve(predict(fit, newdata = data.frame(x)), col = "red", add = TRUE)
  dev.off()
  x$RawIntDen1 <- x$RawIntDen1/Int0_extrap
  x$RawIntDen1
  dat <- data.frame(matrix(ncol=0,nrow=113))
  dat[name] <- x$RawIntDen1
  print(colnames(dat))
  write.table(dat, paste("dat/",paste(paste(name,sep="_"),"dat",sep="."),sep=""), row.names = FALSE, col.names = FALSE, fileEncoding = "ASCII")
  datalist[name] <- dat
}

big_data = do.call(cbind, datalist)
write.table(big_data, "all_RawIntDen1_fit_transfromed.csv", sep = ";", row.names = FALSE)
