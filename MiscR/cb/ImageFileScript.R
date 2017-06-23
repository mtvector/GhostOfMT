img_path <- "~/Desktop/ss HEP2 jan 27-B08_2016012700304_x000000y000000-04x-FL"
dest_path <- "~/code/data/"
oneUp = img_path
bottom = F
queue <- list(img_path)
dir.create(dest_path)
while(length(queue)>0){
  setwd(img_path)
  dfnames <- dir(img_path)[ !grepl(".csv", dir())]
  for(i in 1:length(dfnames)){
    queue[[length(queue)+1]] <- file.path(img_path, dfnames[i])
  }
  img <- dfnames[grepl(".png", dir())]
  if(length(img)>0){
    dir.create(paste0(dest_path,oneUp))
    file.copy(file.path(img_path,img),to = paste0(dest_path,oneUp))
  }
  if(length(queue)>0){
    img_path <- queue[[1]]
  }
  oneUp = strsplit(img_path,"/")[[1]][length(strsplit(img_path,"/")[[1]])]
  while(!dir.exists(queue[[1]])& length(queue)>1){
    queue <- queue[-1]
  }
  queue <- queue[-1]
}
