get_covariates <- function(mo=1,lat=20){
  require(maptools);require(dplyr);require(sp);require(raster);require(reshape2)
  data("wrld_simpl")
  mymonths=c('January','February','March','April','May','June',
             'July','August','September','October','November','December')
  flist <- list.files('ERA5/',pattern='.RDS')
  if('sst.RDS' %in% flist){
    index <- which(flist=='sst.RDS')
    flist <- flist[-index]
  }
  
  var_list <- list()
  for(j in 1:length(flist)){
    vars <- readRDS(paste0('ERA5/',
                           flist[j]))
    
    mydates <- as.Date(attributes(vars)$dimnames[[3]])
    mylats <- as.numeric(attributes(vars)$dimnames[[2]])
    vars <- vars[,which(mylats>=lat),which(months(mydates)==mymonths[mo] & format(mydates,'%Y')<2019)]
    
    rm(temp2)
    for(i in 1:dim(vars)[3]){
      temp1 <- melt(vars[,,i],varname=c('lon','lat'))
      if(exists('temp2')){
        temp2 <- left_join(temp2,temp1,by=c('lon','lat'))
      } else {
        temp2 <- temp1
      }
    }
    
    colnames(temp2)[3:42] <- seq(1979,2018)
    temp2$lon[which(temp2$lon>180)] <- (360-temp2$lon[which(temp2$lon>180)])*(-1)
    coordinates(temp2)=~lon+lat
    gridded(temp2)=TRUE
    var_stack <- stack(temp2)
    crs(var_stack) <- crs(wrld_simpl)
    var_stack <- projectRaster(var_stack,crs=CRS(ease2))
    var_list[[j]] <- projectRaster(var_stack,crs=CRS(ease2))
    
  }
  names(var_list) <- flist
  return(var_list)
  
}








