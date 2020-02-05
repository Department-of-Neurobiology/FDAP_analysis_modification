library(ggplot2)
library(Rmisc)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
old_filenames <- Sys.glob(file.path("*best_fit.dat"))  #parameters
old_filenames

for (i in old_filenames){  
  better_name <- gsub(str_extract(i,"\\d+"),sprintf("%02d",as.numeric(str_extract(i,"\\d+"))),i)
  better_name
  file.rename(i, better_name)
}

filenames <- Sys.glob(file.path("*best_fit.dat"))
filenames
datalist = list()
for (i in filenames){  
  x <- read.table(i, sep = " ", col.names=paste("fit"),fill = TRUE)
  name <- gsub(".dat","",i)
  dat <- data.frame(matrix(ncol=0,nrow=113))
  dat[name] <- x$fit
  print(colnames(dat))
  datalist[name] <- dat
}
big_data = do.call(cbind, datalist)
write.table(big_data, "best_fits_final.csv", sep = ",", row.names = FALSE)

library(tidyr)
plot_fits <- read.table("best_fits_final.csv", sep = ",", header = TRUE)
plot_fits <- cbind(plot_fits, "type"="fits","time"=1:nrow(plot_fits))
long_plot_fits <- gather(plot_fits, cell, best_fit, -c(time:type))
long_plot_fits$cell <- as.factor(long_plot_fits$cell)
long_plot_fits$cell_names <- long_plot_fits$cell
levels(long_plot_fits$cell_names)
levels(long_plot_fits$cell) <- as.factor(c(1:nlevels(long_plot_fits$cell_names)))


### check file name and separators
plot_Int <- read.table("../test_final.csv", sep = ",", header = TRUE)
###
plot_Int <- cbind(plot_Int, "type"="Int","time"=1:nrow(plot_fits))
long_plot_Int <- gather(plot_Int, cell, measured_intensity, -c(time:type))
long_plot_Int$cell <- as.factor(long_plot_Int$cell)
long_plot_Int$cell_names <- long_plot_Int$cell
levels(long_plot_Int$cell) <- as.factor(c(1:nlevels(long_plot_Int$cell_names)))

plot_final <- merge(long_plot_fits, long_plot_Int, by=c("cell","time"))

### check file name and separators
plot_final$construct <- "hTau352"
###
plot_final_summary <- summarySE(plot_final, measurevar="measured_intensity", groupvars=c("construct","time"))
### change se to sd for further plotting
plot_final_summary$CI_lower <- plot_final_summary$measured_intensity + plot_final_summary$se
plot_final_summary$CI_upper <- plot_final_summary$measured_intensity - plot_final_summary$se
###

#svg("ala.svg", width=10, height=6) 
ggplot(plot_final, aes(x = time, y = measured_intensity)) + 
  geom_point(aes(color = cell)) +
  geom_line(aes(y=best_fit, color=cell)) +
  geom_ribbon(data = plot_final_summary, aes(ymin=CI_lower,ymax=CI_upper,fill=construct),alpha=0.4) +#color="grey70"
  theme_classic() +
  ggtitle("Plot of mean measured intensity and best fits by time") +
  xlab("Time (s)") + 
  ylab("Normalized intensity") +
  scale_x_continuous(breaks=seq(0, 110, 10)) +
  scale_y_continuous(limits = c(0,1.1),breaks=seq(0, 1.1, 0.1))
#dev.off() 

### comment after finishing to twick plot
write.table(plot_final_summary, "for_plotting_mean_and_se.csv", sep = ";", row.names = FALSE)
###

plot_together <- read.table("for_plotting_mean_and_se.csv", sep = ",", header = TRUE)
#svg("wt_ala_mean_se.svg", width=10, height=6) 
ggplot(plot_together, aes(x = time, y = measured_intensity)) + 
  geom_line(aes(y=measured_intensity, color=construct)) +
  geom_ribbon(data = plot_together, aes(ymin=CI_lower,ymax=CI_upper,fill=construct),alpha=0.4) +
  theme_classic() +
  ggtitle("Plot of mean measured intensity +-se by time") +
  xlab("Time (s)") + 
  ylab("Normalized intensity") +
  scale_x_continuous(breaks=seq(0, 110, 10)) +
  scale_y_continuous(limits = c(0,1.1),breaks=seq(0, 1.1, 0.1))
#dev.off() 