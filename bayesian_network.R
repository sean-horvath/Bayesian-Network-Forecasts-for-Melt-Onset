################################################################################
# This script uses a Bayesian Network to forecast the mean melt onset date in  #
# a selection of seas in the Arctic using a variety of atmospheric and oceanic #
# principle components as predictors.  Separate models are constructed for     #
# lead times ranging from 1- to 5-months.  Covariates are selected based on    #
# statistically significant correlations with melt onset.                      #
################################################################################


library(rasterVis)
library(reshape2)
library(dplyr)
library(lineup)
library(corrplot)
library(tidyverse)


# Melt Onset Data ---------------------------------------------------------


# Read in the melt onset data and a mask raster to group melt onset by sea
melt <- readRDS('Melt.RDS')
masks <- raster('Masks/SSMI_Regions2.tif')

# Get values of both in a data.frame
melt_vals <- as.data.frame(getValues(melt))
seas <- getValues(masks)

# Relabel the mask factors
seas <- factor(seas,
               levels=c(seq(0,16),20,21),
               labels=c('Lakes','Open Ocean','Okhotsk',
                            'Bering','Hudson','North American',
                            'Labrador','Greenland','Barents',
                            'Kara','Laptev','East Siberian',
                            'Chukchi','Beaufort','CAA',
                            'Central','Baffin','Ignore','Ignore2'))

# Then we'll merge the two data.frames, put in long format, and remove seas
# we're not interested in
melt_vals <- cbind(melt_vals,seas)
colnames(melt_vals) <- c(seq(1979,2018),'Seas')

melt_vals <- melt(melt_vals,
                  id.vars='Seas',
                  variable.name='Year',
                  value.name='Melt')
melt_vals <- melt_vals[!is.na(melt_vals$Melt),]

melt_summary <- melt_vals %>%
  group_by(Year,Seas) %>%
  summarise_all(mean,na.rm=T) %>%
  filter(!Seas%in%c('Lakes','Open Ocean','Ignore','Ignore2'))


# Snow Melt ---------------------------------------------------------------


# Here we're going to test if the timing of snow melt in the northern 
# hemisphere is correlated with melt onset in any of the seas.  For each sea
# we correlate every location of snow melt, then produce maps to visualize.
snowmelt <- stack('snowretreat.tif')
names(snowmelt) <- paste0('YR',seq(1979,2018))

snow_vals <- t(getValues(snowmelt))
unique_seas <- unique(melt_summary$Seas)
rlist <- list()
for(i in 1:length(unique_seas)){
  x <- melt_summary$Melt[which(melt_summary$Seas==unique_seas[i])]
  y <- cor(x,snow_vals,use='pairwise.complete.obs')
  rlist[[i]] <- setValues(raster(snowmelt,layer=1),y)
}

# Create a raster stack with correlation maps, remove statistically
# insignificant values, and plot
snow_cor <- stack(rlist)
snow_cor[snow_cor<0.312 & snow_cor>(-0.312)] <- NA
names(snow_cor) <- unique_seas

mapTheme <- rasterTheme(region=rev(brewer.pal(11,'RdBu')),
                        axis.line=list(col='transparent'))
levelplot(snow_cor,margin=F,
          par.settings=mapTheme,
          main='Correlations with Melt Onset',
          at=seq(-1,1,length.out=21),
          scales=list(draw=FALSE),
          between=list(x=0,y=0),
          colorkey=list(labels=list(cex=1.5)))


# Principle Components ----------------------------------------------------


# Source two functions that will load in all of our potential predictor
# variables
source('get_covariates.R')
source('get_ssts.R')

# Load in a NetCDF file to get the CRS we're interested in
library(ncdf4)

latlon <- nc_open('EASE2_N25km.geolocation.v0.9.nc')
ease2 <- ncatt_get(latlon,'crs','proj4text')$value
nc_close(latlon)

# Now select which month we want to use as a predictor and load in the
# covraiates.  This can take some time (~15 minutes)
mo=10
var_list <- get_covariates(mo=mo,lat=20)
sst_list <- get_ssts(mo=mo)

# Renaming list items
all_list <- var_list
all_list[[length(all_list)+1]] <- sst_list[[1]]
all_list[[length(all_list)+1]] <- sst_list[[2]]
names(all_list)[16:17] <- names(sst_list)

# If we're using a predictor month that is from the previous year, we need to 
# adjusts the datasets so years match up
if(mo %in% c(9,10,11,12)){
  all_list <- lapply(all_list,function(x){
    temp <- dropLayer(x,dim(x)[3])
  })
}

# Now go through each list item and extract the first 3 principle components
all_pcs <- lapply(all_list,function(x){
  y <- getValues(x)
  index <- complete.cases(y)
  y1 <- y[index,]
  pcs <- prcomp(t(y1))
  df <- data.frame(pc1=pcs$x[,1]/pcs$sdev[1],
                   pc2=pcs$x[,2]/pcs$sdev[2],
                   pc3=pcs$x[,3]/pcs$sdev[3])
  return(df)
})

# And combine them all into one matrix
mypcs <- do.call(cbind,all_pcs)

# Here we can take a look at any covariability between PCs
M <- cor(mypcs)
corrplot(M, p.mat = res1$p, sig.level = 0.05)

# Now we'll look at correlations between each PC and melt onset in each sea
# First convert melt data to wide format
melt_wide <- spread(melt_summary,Seas,Melt)
if(mo %in% c(9,10,11,12)){
  melt_wide <- melt_wide[-1,]
  melt_summary <- filter(melt_summary,Year!=1979)
}

# Then we can correlate each column of PC matrix with each column of melt matrix
cor_mat <- corbetw2mat(melt_wide[,-1],mypcs,what='all')
cor_mat[cor_mat<0.312 & cor_mat>(-0.312)] <- 0
corrplot(cor_mat,method='square')

# Put correlations into long format for plotting
df <- melt(cor_mat,varnames=c('Sea','Var'),
           value.name='Corr') %>%
  mutate_if(is.factor,as.character) %>%
  mutate(Var=str_remove(Var,'.RDS'))

# Save our correlation data.frame
df2 <- df %>%
  filter(Corr!=0) %>%
  arrange(Sea,Corr)
write.csv(df2,file='C:/Users/seanm/Documents/Research/BayesianNetworks/pc_correlations_oct.csv')


# Field Means -------------------------------------------------------------

# We can also test correlations with average values of atmospheric variables
# over the same area as each sea melt onsest
# Correlations with field values:
all_field_corr <- lapply(all_list,function(x){
  y <- getValues(x)
  sampras <- raster(x,layer=1)
  unique_seas <- unique(melt_summary$Seas)
  rlist <- list()
  for(i in 1:length(unique_seas)){
    x1 <- melt_summary$Melt[which(melt_summary$Seas==unique_seas[i])]
    y1 <- cor(x1,t(y),use='pairwise.complete.obs')
    rlist[[i]] <- setValues(sampras,y1[1,])
  }
  all_corr <- stack(rlist)
  return(all_corr)
})

# Use "wrld_simpl" data set from maptools package to plot coastlines
data("wrld_simpl")
out <- wrld_simpl
out <- as(wrld_simpl, "SpatialLines")
out <- crop(out, extent(-180, 180, 20, 83.57027))
wrld_ice <- spTransform(out,crs(all_field_corr[[1]]))
detach(package:tidyverse,unload=T)
detach(package:ggplot2,unload=T)

# Loop through each variable and plot correlation maps, facetting by sea
for(i in 1:(length(all_field_corr)-2)){
  x <- all_field_corr[[i]]
  x_name <- names(all_field_corr)[i]
  x_name <- str_sub(x_name,1,(nchar(x_name)-4))
  x[x<0.312 & x>(-0.312)] <- NA
  names(x) <- unique_seas
  mapTheme <- rasterTheme(region=rev(brewer.pal(11,'RdBu')),
                          axis.line=list(col='transparent'))
  
  png(filename=paste0('C:/Users/seanm/Documents/Research/BayesianNetworks/Oct/',
                      x_name,'.png'),
      width=1113,
      height=774)
  p <- levelplot(x,margin=F,
            par.settings=mapTheme,
            main=paste0(x_name,' Correlations with Melt Onset'),
            at=seq(-1,1,length.out=21),
            scales=list(draw=FALSE),
            between=list(x=0,y=0),
            colorkey=list(labels=list(cex=1.5))) +
    layer(sp.polygons(wrld_ice, col="grey40",lwd=0.2),under=F)
  print(p)
  dev.off()
}


# Network Fitting ---------------------------------------------------------


# Create a network with PCs
library(bnlearn)

# Make names a bit neater
df2$Var <- str_replace(df2$Var,coll('.'),'_')
df3 <- df2 %>%
  filter(Sea%in%c('Laptev','Kara','Barents',
                  'Beaufort','Bering','CAA',
                  'Central'))

# Here we're going to create our network structure.  This is done based on 
# significant correlations between PCs and melt onset for each sea.
unique_seas <- unique(df3$Sea)
rlist <- list()
for(i in 1:length(unique_seas)){
  x <- filter(df3,Sea==unique_seas[i])
  temp <- paste(x$Var,sep=':',collapse=':')
  rlist[[i]] <- paste('[',x$Sea[1],'|',temp,']',sep='')
}
arc_string <- paste(unlist(rlist),collapse='')
arc_string <- paste(c(paste(paste0('[',unique(df3$Var),']'),collapse=''),
                arc_string),collapse='')
structure <- empty.graph(c(unique(df3$Sea),unique(df3$Var)))
modelstring(structure) <- arc_string

# Let's take a look at our network
# This plotting function was taken from an R-Bloggers post written by Daniel 
# Oehm (www.r-bloggers.com/bayesian-network-example-with-the-bnlearn-package/)
library(visNetwork)
plot.network <- function(structure, ht = "700px"){
  nodes.uniq <- unique(c(structure$arcs[,1], structure$arcs[,2]))
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "darkturquoise",
                      shadow = TRUE)
  
  edges <- data.frame(from = structure$arcs[,1],
                      to = structure$arcs[,2],
                      arrows = "to",
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "black")
  
  return(visNetwork(nodes, edges, height = ht, width = "100%"))
}

# observe structure
plot.network(structure)

# Clean up names
colnames(mypcs) <- str_replace(colnames(mypcs),
                               coll('.RDS.'),
                               '_')
colnames(mypcs) <- str_replace(colnames(mypcs),
                               coll('.'),
                               '_')

# Create our fitted data matrix
df_fit <- cbind(mypcs[,colnames(mypcs)%in%unique(df3$Var)],
                melt_wide[,colnames(melt_wide)%in%unique(df3$Sea)])

# Learn model parameters
bn.mod <- bn.fit(structure, data = df_fit)

# View our model
bn.mod
bn.fit.qqplot(bn.mod)
bn.fit.histogram(bn.mod)
bn.fit.xyplot(bn.mod)

# Or we can allow bnlearn to learn the model structure.  We'll neet to create
# a "blacck list" to prevent seas from being parent nodes
unique_seas <- unique(df3$Sea)
unique_var <- unique(df3$Var)
black_df <- data.frame(from=rep(unique_seas,each=length(unique_var)),
                       to=rep(unique_var,length(unique_seas)))
melt.hc <- hc(df_fit,blacklist=black_df)

plot.network(melt.hc)

# bn_cross_val <- bn.cv(df_fit,melt.hc)

# We'll test the predictive power by doing a rolling-validation.  We'll drop
# the last 10 years, fit the data, then predict the following year.  We'll then
# add that year to the fit, and predict the next year, continuing throught the
# final year
df_tot <- df_fit
df_fit <- df_tot[1:30,]
df_test <- df_tot[31,]
melt.hc <- hc(df_fit,blacklist=black_df)
fitted <- bn.fit(melt.hc,df_fit)
plot.network(melt.hc)

mynodes <- c('Central','CAA','Laptev','Beaufort',
             'Kara','Barents','Bering')
yrs <- seq(2009,2017)
df_pred <- data.frame(Sea=rep(mynodes,length(yrs)),
                      Year=rep(yrs,each=length(mynodes)),
                      Pred=NA,Obs=NA)

for(i in 1:length(mynodes)){
  for(j in 1:length(yrs)){
    df_fit <- df_tot[1:(29+j),]
    df_test <- df_tot[(30+j),]
    
    melt.hc <- hc(df_fit,blacklist=black_df)
    fitted <- bn.fit(melt.hc,df_fit)

    df_pred$Pred[which(df_pred$Sea==mynodes[i] & df_pred$Year==yrs[j])] <- predict(fitted,mynodes[i],df_test)
    df_pred$Obs[which(df_pred$Sea==mynodes[i] & df_pred$Year==yrs[j])] <- as.numeric(df_test[mynodes[i]])
  }
}

# View predicted and observed correlations
df_pred %>%
  group_by(Sea) %>%
  mutate(cor=cor(Pred,Obs)) %>%
  filter(Year==2009)

# Put back in long format for plotting
df_pred <- df_pred %>%
  gather('Type','Value',c(-Sea,-Year))

# Plot predicted and observed values
library(ggplot2)
ggplot(df_pred,aes(x=Year,y=Value,color=Sea)) +
  geom_line(aes(linetype=Type))

# And save our full dataset data.frame
saveRDS(df_tot,'C:/Users/seanm/Documents/Research/BayesianNetworks/df_tot_oct.RDS')



# Try a network with all months as predictors -----------------------------

# Read in data from all predictor months
df_files <- c('df_tot_oct.RDS',
              'df_tot_nov.RDS',
              'df_tot_dec.RDS',
              'df_tot_jan.RDS',
              'df_tot_feb.RDS',
              'df_tot_mar.RDS')
df_list <- list()
for(i in 1:length(df_files)){
  df_list[[i]] <- readRDS(df_files[i])
  names(df_list)[i] <- str_sub(df_files[i],8,10)
}

# create individual networks first:
df_fit <- df_list$oct
unique_var <- colnames(df_fit)[grep('_',colnames(df_fit))]
unique_seas <- colnames(df_fit)[!colnames(df_fit) %in% unique_var]

black_df <- data.frame(from=rep(unique_seas,each=length(unique_var)),
                       to=rep(unique_var,length(unique_seas)))

melt.hc <- hc(df_fit,blacklist=black_df)
plot.network(melt.hc)

# bn_cross_val <- bn.cv(df_fit,melt.hc)

# try predictions:
df_tot <- df_fit
mynodes <- c('Central','CAA','Laptev','Beaufort',
             'Kara','Barents','Bering')
yrs <- seq(2009,2018)
df_pred <- data.frame(Sea=rep(mynodes,length(yrs)),
                      Year=rep(yrs,each=length(mynodes)),
                      Pred=NA,Obs=NA)

for(i in 1:length(mynodes)){
  for(j in 1:length(yrs)){
    df_fit <- df_tot[1:(28+j),]
    df_test <- df_tot[(29+j),]
    
    melt.hc <- hc(df_fit,blacklist=black_df)
    fitted <- bn.fit(melt.hc,df_fit)
    
    df_pred$Pred[which(df_pred$Sea==mynodes[i] & df_pred$Year==yrs[j])] <- predict(fitted,mynodes[i],df_test)
    df_pred$Obs[which(df_pred$Sea==mynodes[i] & df_pred$Year==yrs[j])] <- as.numeric(df_test[mynodes[i]])
  }
}

df_pred %>%
  group_by(Sea) %>%
  mutate(cor=cor(Pred,Obs)) %>%
  filter(Year==2009)

df_pred <- df_pred %>%
  gather('Type','Value',c(-Sea,-Year))

library(ggplot2)
ggplot(df_pred,aes(x=Year,y=Value,color=Sea)) +
  geom_line(aes(linetype=Type)) +
  ggtitle('October')

# combine and do one big network:

df_list <- list()
for(i in 1:length(df_files)){
  temp <- readRDS(paste0('BayesianNetworks/',df_files[i]))
  if(dim(temp)[1]==40){temp <- temp[-1,]}
  if(i==length(df_files)){
    temp2 <- temp[,mynodes]
  }
  temp <- temp[,!colnames(temp) %in% mynodes]
  colnames(temp) <- paste0(colnames(temp),'_',str_sub(df_files[i],8,10))
  df_list[[i]] <- temp
  if(i==length(df_files)){
    df_list[[i+1]] <- temp2
  }
}

df_tot <- do.call(cbind,df_list)
col_names <- colnames(df_tot)[!colnames(df_tot) %in% mynodes]

black_df <- data.frame(from=rep(mynodes,each=length(col_names)),
                       to=rep(col_names,length(mynodes)))
mos <- str_sub(df_files,8,10)
black_df_list <- list()
for(i in 2:length(mos)){
  x <- col_names[grep(mos[i],col_names)]
  y <- col_names[grep(paste(mos[1:(i-1)],collapse='|'),col_names)]
  black_df_list[[(i-1)]] <- data.frame(from=rep(x,each=length(y)),
                                       to=rep(y,length(x)))
}
black_df <- rbind(black_df,
                  do.call(rbind,black_df_list))
melt.hc <- hc(df_tot,blacklist=black_df)
plot.network(melt.hc)

mynodes <- c('Central','CAA','Laptev','Beaufort',
             'Kara','Barents','Bering')
yrs <- seq(2009,2018)
df_pred <- data.frame(Sea=rep(mynodes,length(yrs)),
                      Year=rep(yrs,each=length(mynodes)),
                      Pred=NA,Obs=NA)

for(i in 1:length(mynodes)){
  for(j in 1:length(yrs)){
    df_fit <- df_tot[1:(28+j),]
    df_test <- df_tot[(29+j),]
    
    melt.hc <- hc(df_fit,blacklist=black_df)
    fitted <- bn.fit(melt.hc,df_fit)
    
    df_pred$Pred[which(df_pred$Sea==mynodes[i] & df_pred$Year==yrs[j])] <- predict(fitted,mynodes[i],df_test)
    df_pred$Obs[which(df_pred$Sea==mynodes[i] & df_pred$Year==yrs[j])] <- as.numeric(df_test[mynodes[i]])
  }
}

df_pred %>%
  group_by(Sea) %>%
  mutate(cor=cor(Pred,Obs)) %>%
  filter(Year==2009)

df_pred <- df_pred %>%
  gather('Type','Value',c(-Sea,-Year))

library(ggplot2)
ggplot(df_pred,aes(x=Year,y=Value,color=Sea)) +
  geom_line(aes(linetype=Type)) +
  ggtitle('Full Network')






