setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
filenames <- Sys.glob(file.path("*RawIntDen1.csv"))  #parameters
filenames
datalist = list()
for (i in filenames){  
  x <- read.table(i, sep = ",", header = TRUE, fill = TRUE, row.names = "X")
  name <- gsub(".csv","",i)
  dat <- data.frame(matrix(ncol=0,nrow=113))
  dat[name] <- x$RawIntDen1
  print(colnames(dat))
  datalist[name] <- dat
}
big_data = do.call(cbind, datalist)
write.table(big_data, "all_RawIntDen1_assembled.csv", sep = ";", row.names = FALSE)

