########################################
### Tree growth response to competition and health condition in a 
### subtropical savanna
###
### Neighborhood Calculations
###
### D. Alex Bowers
### Daniel J. Johnson
########################################

# required packages
library(data.table)
library(dplyr)

# needed data
trees <- read.csv("data/trees_neighborhood.csv")



# function for calculating neighborhood
neighborhood <- function(dt,r,f){ ## r = radius, dt neighborhood trees, f = focal stems
  dt <- data.table(dt)
  # matrices are a lot faster than data.frames so we create one to hold our results from loop
  temp <- NULL
  no.col <- 7   #  is the number of index values
  temp.mat <- matrix(nrow=nrow(f) , ncol=no.col)
  
  for(i in 1:nrow(f)){
    
    data.row <- f[i,]
    spp <- data.row$sp
    focal.tag <- data.row$StemTag
    x.row <- data.row$gx
    y.row <- data.row$gy
    status <- data.row$status_2024
    ba <- data.row$ba
    
    square <- dt[gx <= x.row+r & gx >= x.row-r & gy <= y.row+r & gy >= y.row-r,]
    
    square$dist <- sqrt((( square$gx - x.row )^2 + ( square$gy - y.row )^2))   # calculate distance between neighbors and focal tree  
    
    # this part uses data.table syntax notice the generic and familial references
    neighb <- square[dist<= r & square$StemTag != focal.tag,]   # remove focal tree and only keep individuals that are with "radius" meters of the focal tree 
    
    # remove neighbors
    neighb <- subset(neighb, neighb$dist != 0)
    
    neighb[,ch:= ifelse(sp==spp, 1,0)]  # this uses data.table notation
    
    neighb[,status:= ifelse(status=="A", 1,0)] ## stat
    
    neighb[,ch:= ifelse(sp == spp, "c","h")]  # this uses data.table notation (c = conspecificm h = heterospecific )
    BAcon <- neighb[ch == "c", sum(ba, na.rm = T)]    #sum basal area of conspecific neighbors
    BAhet <- neighb[ch == "h", sum(ba, na.rm=T)]    #sum BA of heterospec conmyco neighbors
    tBA <- sum(neighb$ba, na.rm=T)                                 # total basal area
    nT <- nrow(neighb)
    
    # distance weighting
    idw_con <- sum(ifelse(neighb$ch=="c", neighb$ba/neighb$dist, 0), na.rm = T) # distance weighted
    idw_het <- sum(ifelse(neighb$ch=="h", neighb$ba/neighb$dist, 0), na.rm = T) # distance weighted
    idw_tot <- sum(neighb$ba/neighb$dist, na.rm = T)                             # distance weighted
    
    temp.mat[i, ] <- c(BAcon, BAhet, tBA, nT, 
                       idw_con, idw_het, idw_tot) # put them in the matrix
    
    
    
  }
  temp.df <- data.frame(temp.mat) 
  names(temp.df) <- c("BAcon","BAhet","tBA","nT",
                      "idw_con","idw_het","idw_tot")
  
  return(data.frame(f, temp.df)) # return the matrix
}



# Determine what distance radius is the best model
## 5 m
trees_5 <- neighborhood(trees, r = 5, trees)
trees_5_mod <- lm(bai ~ idw_tot, data = trees_5)
## 10 m 
trees_10 <- neighborhood(trees, r = 10, trees)
trees_10_mod <- lm(bai ~ idw_tot, data = trees_10)
## 15 m
trees_15 <- neighborhood(trees, r = 15, trees)
trees_15_mod <- lm(bai ~ idw_tot, data = trees_15)
## 20 m
trees_20 <- neighborhood(trees, r = 20, trees)
trees_20_mod <- lm(bai ~ idw_tot, data = trees_20)
## 25 m
trees_25 <- neighborhood(trees, r = 25, trees)
trees_25_mod <- lm(bai ~ idw_tot, data = trees_25)

# find best AIC
AIC(trees_5_mod, trees_10_mod, trees_15_mod, trees_20_mod, trees_25_mod)

# Use best AIC and select trees that have neighbors inside the plot
trees_PCoA <- trees_15 %>% filter(gx >= 15 & gx <= 465, gy >= 15 & gy <= 465)

# data is now ready for the PCoA
trees_PCoA