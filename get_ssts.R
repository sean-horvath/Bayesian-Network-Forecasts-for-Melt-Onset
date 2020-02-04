get_ssts <- function(mo=1){
  require(maptools);require(dplyr);require(sp);require(raster);require(reshape2)
  data("wrld_simpl")
  mymonths=c('January','February','March','April','May','June',
             'July','August','September','October','November','December')
  vars <- readRDS('ERA5/sst.RDS')
  
  mydates <- as.Date(attributes(vars)$dimnames[[3]])
  mylats <- as.numeric(attributes(vars)$dimnames[[2]])
  vars <- vars[,which(mylats>=(-20)),which(months(mydates)==mymonths[mo] & format(mydates,'%Y')<2019)]
  
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
  

  sst_list <- list()
  # Pacific, north of -20
  box1 <- which(temp2[,1]<=-70 & temp2[,1]>=-80 & temp2[,2]>=(-20) & temp2[,2]<=9)
  box2 <- which(temp2[,1]<=-80 & temp2[,1]>=-84 & temp2[,2]>=(-20) & temp2[,2]<=9)
  box3 <- which(temp2[,1]<=-84 & temp2[,1]>=-90 & temp2[,2]>=(-20) & temp2[,2]<=14)
  box4 <- which(temp2[,1]<=-90 & temp2[,1]>=-100 & temp2[,2]>=(-20) & temp2[,2]<=18)
  box5 <- which(temp2[,1]<=-100 & temp2[,1]>=-145 & temp2[,2]>=(-20) & temp2[,2]<=66)
  box6 <- which(temp2[,1]<=145 & temp2[,1]>=100 & temp2[,2]>=(-20) & temp2[,2]<=66)
  box7 <- which(temp2[,1]<=-145 & temp2[,2]>=(-20) & temp2[,2]<=66)
  box8 <- which(temp2[,1]>=145 & temp2[,2]>=(-20) & temp2[,2]<=66)
  
  pacific_index <- c(box1,box2,box3,box4,box5,box6,box7,box8)
  temp <- temp2[pacific_index,]
  
  coordinates(temp)=~lon+lat
  gridded(temp)=TRUE
  sst_pac_stack <- stack(temp)
  crs(sst_pac_stack) <- crs(wrld_simpl)
  sst_list[[1]] <- sst_pac_stack
  
  # Atlantic, north of -20
  box1 <- which(temp2[,1]<=45 & temp2[,1]>=30 & temp2[,2]>=31 & temp2[,2]<=60)
  box2 <- which(temp2[,1]<=30 & temp2[,1]>=20 & temp2[,2]>=31 & temp2[,2]<=65)
  box3 <- which(temp2[,1]<=20 & temp2[,1]>=-65 & temp2[,2]>=(-20) & temp2[,2]<=65)
  box4 <- which(temp2[,1]<=-65 & temp2[,1]>=-70 & temp2[,2]>=(-20) & temp2[,2]<=55)
  box5 <- which(temp2[,1]<=-70 & temp2[,1]>=-75 & temp2[,2]>=(-20) & temp2[,2]<=55)
  box6 <- which(temp2[,1]<=-75 & temp2[,1]>=-84 & temp2[,2]>=9 & temp2[,2]<=40)
  box7 <- which(temp2[,1]<=-84 & temp2[,1]>=-90 & temp2[,2]>=14 & temp2[,2]<=40)
  box8 <- which(temp2[,1]<=-90 & temp2[,1]>=-100 & temp2[,2]>=18 & temp2[,2]<=40)
  
  atlantic_index <- c(box1,box2,box3,box4,box5,box6,box7,box8)
  temp <- temp2[atlantic_index,]
  
  coordinates(temp)=~lon+lat
  gridded(temp)=TRUE
  sst_atl_stack <- stack(temp)
  crs(sst_atl_stack) <- crs(wrld_simpl)
  sst_list[[2]] <- sst_atl_stack

  names(sst_list) <- c('Pacific','Atlantic')
  return(sst_list)
  
}








