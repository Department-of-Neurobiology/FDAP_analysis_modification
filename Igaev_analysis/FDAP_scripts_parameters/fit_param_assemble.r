library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

old_filenames <- Sys.glob(file.path("*fit_params.dat"))  #parameters
old_filenames

for (i in old_filenames){  
  better_name <- gsub(str_extract(i,"\\d+"),sprintf("%02d",as.numeric(str_extract(i,"\\d+"))),i)
  better_name
  file.rename(i, better_name)
}

filenames <- Sys.glob(file.path("*fit_params.dat"))  #parameters
filenames
chisqlist = list()
konlist = list()
kofflist = list()
boundlist = list()

for (i in filenames){  
  x <- read.table(i, sep = " ", col.names=paste("column", 0:3, sep="_"),fill = TRUE, row.names = "column_0")
  name <- gsub(".dat","",i)
  chisqlist[name] <- x["chisq/dof","column_1"]
  konlist[name] <- x["kon_fit","column_1"]
  kofflist[name] <- x["koff_fit","column_1"]
  boundlist[name] <- x["bound","column_1"]
}
chisq = t(do.call(cbind, chisqlist))
colnames(chisq) <- "chisq"
kon = t(do.call(cbind, konlist))
colnames(kon) <- "kon"
koff = t(do.call(cbind, kofflist))
colnames(koff) <- "koff"
bound = t(do.call(cbind, boundlist))
colnames(bound) <- "bound"

final_table <- cbind(kon,koff,bound,chisq)
summary(final_table)
write.table(final_table, "fit_params_table.csv",dec = ".", sep = ";", col.names = NA) 



 