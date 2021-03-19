#######
#Libraries
#######
#install.packages("jpeg")
library(jpeg)

#######
#Set working directory to the file location in RStudio
#######
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
output_dir <- paste("output_Reyher_fits", sep="")
dir.create(output_dir)

#######
#Functions
#######
plot_jpeg = function(path, add=FALSE)
{
  require('jpeg')
  jpg = readJPEG(path, native=T) # read the file
  res = dim(jpg)[2:1] # get the resolution, [x, y]
  if (!add) # initialize an empty plot area if add==FALSE
    plot(1,1,xlim=c(1,res[1]),ylim=c(1,res[2]),asp=1,type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  rasterImage(jpg,1,1,res[1],res[2])
}

#######
#Showing figures from subfolders
#######
#check on one
#plot_jpeg('20200819_tub_orange_d1_cell014/ROI.jpg')

#read subfolders
subfolders_with_home <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
subfolders <- subfolders_with_home[-c(1)]

#loop through all
for (subfolder in subfolders) {
  path_to_file <- paste(subfolder, '/NewFitlambdafix.jpg', sep="")
  if (file.exists(path_to_file)) {
    plot_jpeg(path_to_file)
    file.copy(path_to_file, output_dir)
    path_to_output_file <- paste(output_dir, '/NewFitlambdafix.jpg', sep="")
    subfolder_name <- gsub("\\.", "",subfolder)
    file.rename(from = path_to_output_file, to = paste(output_dir,subfolder_name,'_ROI.jpg',sep=""))
  }
  else {
    print(paste("The image file does not exist in", subfolder))
  }
}  

