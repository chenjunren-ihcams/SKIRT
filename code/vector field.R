library(readxl)
library(reshape2)
library(sqldf)


## Data preparation 

setwd("C:/Users/hp/Desktop/?????ؽ?-��??ͶӰ/upload files/")

immune.profiles.tsne <- read_excel("./data/immune.profiles.tsne.xlsx", guess_max = 11150)
patients <- unique(immune.profiles.tsne$ID)

immune.profiles.tsne$tSNE.1 <- -immune.profiles.tsne$tSNE.1
immune.profiles.tsne$tSNE.2 <- -immune.profiles.tsne$tSNE.2

x <- seq(-60,60,2)
y <- seq(-60,60,2)

## define function

point.in.pol <- function(pol.x, pol.y, point.x, point.y){
  #?ο?https://blog.csdn.net/lynon/article/details/82015834
  if (point.x < min(pol.x) | point.x > max(pol.x) | point.y < min(pol.y) | point.y > max(pol.y)) {
    return(FALSE)
  }
  c <- F
  for (i in 2:length(pol.x)) {
    if (
      ((point.y > pol.y[i - 1]) != (point.y > pol.y[i])) &
      point.x < ((pol.x[i] - pol.x[i - 1]) / (pol.y[i] - pol.y[i - 1]) * (point.y - pol.y[i - 1]) + pol.x[i - 1])
    ) {
      c <- !c
    }
  }
  return(c)
}

outline.func <- function(x, y, prob = 0.55, char = T){
  # input??coordinate of original points
  # output??coordinate of profile line
  if (!requireNamespace("MASS", quietly = T)) {
    install.packages("MASS")
  }
  kde2d.data <- MASS::kde2d(x, y, lims = c(floor(range(x)[1]) - ceiling(diff(range(x) / 10)), 
                                           ceiling(range(x)[2]) + ceiling(diff(range(x) / 10)), 
                                           floor(range(y)[1]) - ceiling(diff(range(y) / 10)),
                                           ceiling(range(y)[2]) + ceiling(diff(range(y) / 10))))
  lines.list <- contourLines(kde2d.data, nlevels = 1, levels = quantile(kde2d.data$z, probs = prob),)
  
  x <- round(lines.list[[1]]$x, 3)
  y <- round(lines.list[[1]]$y, 3)
  out <- cbind(x, y)
  out
}
profile <- outline.func(immune.profiles.tsne$tSNE.1,
                        immune.profiles.tsne$tSNE.2,
                        prob = 0.55,
                        char = T)




## vector field: from 0-1 to 1-3 (months)

from0to1_to_1to3_n <- array(0, dim=c(length(x), length(y)))
from0to1_to_1to3_x <- array(0, dim=c(length(x), length(y)))
from0to1_to_1to3_y <- array(0, dim=c(length(x), length(y)))

for (i in 1:length(patients)) {
  patient <- patients[i]
  I <- immune.profiles.tsne$ID==patient &
    immune.profiles.tsne$`days post-transplant` > 0*30 &
    immune.profiles.tsne$`days post-transplant` <= 1*30
  loc.0to1.x <- median(immune.profiles.tsne$tSNE.1[I])
  loc.0to1.y <- median(immune.profiles.tsne$tSNE.2[I])
  II <- immune.profiles.tsne$ID==patient &
    immune.profiles.tsne$`days post-transplant` > 1*30 &
    immune.profiles.tsne$`days post-transplant` <= 3*30
  loc.1to3.x <- median(immune.profiles.tsne$tSNE.1[II])
  loc.1to3.y <- median(immune.profiles.tsne$tSNE.2[II])
  
  ii <- floor((loc.0to1.x+60)/2)+1
  jj <- floor((loc.0to1.y+60)/2)+1
  
  if(sum(I)*sum(II) > 0){
    from0to1_to_1to3_n[ii,jj] <- from0to1_to_1to3_n[ii,jj] + 1
    from0to1_to_1to3_x[ii,jj] <- from0to1_to_1to3_n[ii,jj] + loc.1to3.x - loc.0to1.x
    from0to1_to_1to3_y[ii,jj] <- from0to1_to_1to3_n[ii,jj] + loc.1to3.y - loc.0to1.y
  }
  
}

from0to1_to_1to3_x <-  from0to1_to_1to3_x / from0to1_to_1to3_n
from0to1_to_1to3_y <-  from0to1_to_1to3_y / from0to1_to_1to3_n  
from0to1_to_1to3_x[is.na(from0to1_to_1to3_x) |from0to1_to_1to3_n < 4] <- NA
from0to1_to_1to3_y[is.na(from0to1_to_1to3_y) |from0to1_to_1to3_n < 4] <- NA  

## vector field: from 1-3 to 3-6  (months)

from1to3_to_3to6_n <- array(0, dim=c(length(x), length(y)))
from1to3_to_3to6_x <- array(0, dim=c(length(x), length(y)))
from1to3_to_3to6_y <- array(0, dim=c(length(x), length(y)))

for (i in 1:length(patients)) {
  patient <- patients[i]
  I <- immune.profiles.tsne$ID==patient &
    immune.profiles.tsne$`days post-transplant` > 0*30 &
    immune.profiles.tsne$`days post-transplant` <= 3*30
  loc.1to3.x <- median(immune.profiles.tsne$tSNE.1[I])
  loc.1to3.y <- median(immune.profiles.tsne$tSNE.2[I])
  II <- immune.profiles.tsne$ID==patient &
    immune.profiles.tsne$`days post-transplant` > 3*30 &
    immune.profiles.tsne$`days post-transplant` <= 6*30
  loc.3to6.x <- median(immune.profiles.tsne$tSNE.1[II])
  loc.3to6.y <- median(immune.profiles.tsne$tSNE.2[II])
  
  ii <- floor((loc.1to3.x+60)/2)+1
  jj <- floor((loc.1to3.y+60)/2)+1
  if(sum(I)*sum(II) > 0){
    from1to3_to_3to6_n[ii,jj] <- from1to3_to_3to6_n[ii,jj] + 1
    from1to3_to_3to6_x[ii,jj] <- from1to3_to_3to6_n[ii,jj] + loc.3to6.x - loc.1to3.x
    from1to3_to_3to6_y[ii,jj] <- from1to3_to_3to6_n[ii,jj] + loc.3to6.y - loc.1to3.y
  }
}

from1to3_to_3to6_x <-  from1to3_to_3to6_x / from1to3_to_3to6_n
from1to3_to_3to6_y <-  from1to3_to_3to6_y / from1to3_to_3to6_n  
from1to3_to_3to6_x[is.na(from1to3_to_3to6_x) |from1to3_to_3to6_n < 4] <- NA
from1to3_to_3to6_y[is.na(from1to3_to_3to6_y) |from1to3_to_3to6_n < 4] <- NA  


## vector field: from 3-6 to 6-12 (months)

from3to6_to_6to12_n <- array(0, dim=c(length(x), length(y)))
from3to6_to_6to12_x <- array(0, dim=c(length(x), length(y)))
from3to6_to_6to12_y <- array(0, dim=c(length(x), length(y)))

for (i in 1:length(patients)) {
  patient <- patients[i]
  I <- immune.profiles.tsne$ID==patient &
    immune.profiles.tsne$`days post-transplant` > 3*30 &
    immune.profiles.tsne$`days post-transplant` <= 6*30
  loc.3to6.x <- median(immune.profiles.tsne$tSNE.1[I])
  loc.3to6.y <- median(immune.profiles.tsne$tSNE.2[I])
  II <- immune.profiles.tsne$ID==patient &
    immune.profiles.tsne$`days post-transplant` > 6*30 &
    immune.profiles.tsne$`days post-transplant` <= 12*30
  loc.6to12.x <- median(immune.profiles.tsne$tSNE.1[II])
  loc.6to12.y <- median(immune.profiles.tsne$tSNE.2[II])
  
  ii <- floor((loc.3to6.x+60)/2)+1
  jj <- floor((loc.3to6.y+60)/2)+1
  
  if(sum(I)*sum(II) > 0){
    from3to6_to_6to12_n[ii,jj] <- from3to6_to_6to12_n[ii,jj] + 1
    from3to6_to_6to12_x[ii,jj] <- from3to6_to_6to12_n[ii,jj] + loc.6to12.x - loc.3to6.x
    from3to6_to_6to12_y[ii,jj] <- from3to6_to_6to12_n[ii,jj] + loc.6to12.y - loc.3to6.y
  }
}

from3to6_to_6to12_x <-  from3to6_to_6to12_x / from3to6_to_6to12_n
from3to6_to_6to12_y <-  from3to6_to_6to12_y / from3to6_to_6to12_n  
from3to6_to_6to12_x[is.na(from3to6_to_6to12_x) | from3to6_to_6to12_n < 3] <- NA
from3to6_to_6to12_y[is.na(from3to6_to_6to12_y) | from3to6_to_6to12_n < 3] <- NA  



## vector field: from 6-12 to 12-18 (months)

from6to12_to_12to18_n <- array(0, dim=c(length(x), length(y)))
from6to12_to_12to18_x <- array(0, dim=c(length(x), length(y)))
from6to12_to_12to18_y <- array(0, dim=c(length(x), length(y)))

for (i in 1:length(patients)) {
  patient <- patients[i]
  I <- immune.profiles.tsne$ID==patient &
    immune.profiles.tsne$`days post-transplant` > 6*30 &
    immune.profiles.tsne$`days post-transplant` <= 12*30
  loc.6to12.x <- median(immune.profiles.tsne$tSNE.1[I])
  loc.6to12.y <- median(immune.profiles.tsne$tSNE.2[I])
  II <- immune.profiles.tsne$ID==patient &
    immune.profiles.tsne$`days post-transplant` > 12*30 &
    immune.profiles.tsne$`days post-transplant` <= 18*30
  loc.12to18.x <- median(immune.profiles.tsne$tSNE.1[II])
  loc.12to18.y <- median(immune.profiles.tsne$tSNE.2[II])
  
  ii <- floor((loc.6to12.x+60)/2)+1
  jj <- floor((loc.6to12.y+60)/2)+1
  
  if(sum(I)*sum(II) > 0){
    from6to12_to_12to18_n[ii,jj] <- from6to12_to_12to18_n[ii,jj] + 1
    from6to12_to_12to18_x[ii,jj] <- from6to12_to_12to18_n[ii,jj] + loc.12to18.x - loc.6to12.x
    from6to12_to_12to18_y[ii,jj] <- from6to12_to_12to18_n[ii,jj] + loc.12to18.y - loc.6to12.y
  }
}

from6to12_to_12to18_x <-  from6to12_to_12to18_x / from6to12_to_12to18_n
from6to12_to_12to18_y <-  from6to12_to_12to18_y / from6to12_to_12to18_n
from6to12_to_12to18_x[is.na(from6to12_to_12to18_x) | from6to12_to_12to18_n < 3] <- NA
from6to12_to_12to18_y[is.na(from6to12_to_12to18_y) | from6to12_to_12to18_n < 3] <- NA



## vector field: from 12-18 to 18-24 (months)

from12to18_to_18to24_n <- array(0, dim=c(length(x), length(y)))
from12to18_to_18to24_x <- array(0, dim=c(length(x), length(y)))
from12to18_to_18to24_y <- array(0, dim=c(length(x), length(y)))

for (i in 1:length(patients)) {
  patient <- patients[i]
  I <- immune.profiles.tsne$ID==patient &
    immune.profiles.tsne$`days post-transplant` > 12*30 &
    immune.profiles.tsne$`days post-transplant` <= 18*30
  loc.12to18.x <- median(immune.profiles.tsne$tSNE.1[I])
  loc.12to18.y <- median(immune.profiles.tsne$tSNE.2[I])
  II <- immune.profiles.tsne$ID==patient &
    immune.profiles.tsne$`days post-transplant` > 18*30 &
    immune.profiles.tsne$`days post-transplant` <= 24*30
  loc.18to24.x <- median(immune.profiles.tsne$tSNE.1[II])
  loc.18to24.y <- median(immune.profiles.tsne$tSNE.2[II])
  
  ii <- floor((loc.12to18.x+60)/2)+1
  jj <- floor((loc.12to18.y+60)/2)+1
  if(sum(I)*sum(II) > 0){
    from12to18_to_18to24_n[ii,jj] <- from12to18_to_18to24_n[ii,jj] + 1
    from12to18_to_18to24_x[ii,jj] <- from12to18_to_18to24_n[ii,jj] + loc.18to24.x - loc.12to18.x
    from12to18_to_18to24_y[ii,jj] <- from12to18_to_18to24_n[ii,jj] + loc.18to24.y - loc.12to18.y
  }
}

from12to18_to_18to24_x <-  from12to18_to_18to24_x / from12to18_to_18to24_n
from12to18_to_18to24_y <-  from12to18_to_18to24_y / from12to18_to_18to24_n
from12to18_to_18to24_x[is.na(from12to18_to_18to24_x) | from12to18_to_18to24_n < 3] <- NA
from12to18_to_18to24_y[is.na(from12to18_to_18to24_y) | from12to18_to_18to24_n < 3] <- NA


## average vector field

fromstart_to_end_x <- array(NA, dim=c(61,61))
fromstart_to_end_y <- array(NA, dim=c(61,61))
for (ii in 1:61) {
  for (jj in 1:61) {
    temp <- c(  from0to1_to_1to3_x[ii,jj],
                from1to3_to_3to6_x[ii,jj],
                from3to6_to_6to12_x[ii,jj],
                from6to12_to_12to18_x[ii,jj],
                from12to18_to_18to24_x[ii,jj])
    fromstart_to_end_x[ii,jj] <- mean(temp,na.rm = T)
  }
}

for (ii in 1:61) {
  for (jj in 1:61) {
    temp <- c(  from0to1_to_1to3_y[ii,jj],
                from1to3_to_3to6_y[ii,jj],
                from3to6_to_6to12_y[ii,jj],
                from6to12_to_12to18_y[ii,jj],
                from12to18_to_18to24_y[ii,jj])
    fromstart_to_end_y[ii,jj] <- mean(temp,na.rm = T)
  }
}


## average vector field (re-scale)

xx <- c()
yy <- c()

for (i in 1 : nrow(fromstart_to_end_x)) {
  for (j in 1: ncol(fromstart_to_end_x)) {
    x1 <- x[i] + 1
    y1 <- y[j] + 1
    xx <- rbind(xx, c(x1, y1, fromstart_to_end_x[i,j]))
    yy <- rbind(yy, c(x1, y1, fromstart_to_end_y[i,j]))
  }
}

xx <- xx[!is.nan(xx[,3]),]
yy <- yy[!is.nan(yy[,3]),]

x <- -12:12*5  #re-scale
y <- -12:12*5  #re-scale
fromstart_to_end_xx <- array(NA,dim = c(length(x),length(y)))
fromstart_to_end_yy <- array(NA,dim = c(length(x),length(y)))


for(i in 1:length(x)){
  for (j in 1:length(y)) {
    distance <- (x[i] - xx[,1]) ^ 2 + (y[j] - xx[,2]) ^ 2
    fromstart_to_end_xx[i,j] <- mean(xx[distance < quantile(distance, 0.05), 3])
    fromstart_to_end_yy[i,j] <- mean(yy[distance < quantile(distance, 0.05), 3])
    
  }
}


for(i in -60:60){
  for (j in -60:60) {
    if (!is.na(fromstart_to_end_xx[i/5+13,j/5+13])){
      if (point.in.pol(profile[,1],profile[,2],i,j)==0){
        fromstart_to_end_xx[i/5+13, j/5+13] <- NA
      }
    }
  }
} 

for(i in -60:60){
  for (j in -60:60) {
    if (!is.na(fromstart_to_end_yy[i/5+13,j/5+13])){
      if (point.in.pol(profile[,1],profile[,2],i,j)==0){
        fromstart_to_end_yy[i/5+13, j/5+13] <- NA
      }
    }
  }
} 


## plot 'vector field'
plot(-100, -100, xlim=c(-60,60), ylim=c(-60,60), main="", xlab = "", ylab = "", axes = F)
lines(-profile[,1], -profile[,2], type = 'l') 

equilibrium <- c()
for (ii in 1:length(x)) {
  loc.start.x <- x[ii] + 2.5
  for (jj in 1:length(y)) {
    loc.start.y <- y[jj] + 2.5
    
    if (!is.na(loc.start.x + loc.start.x + fromstart_to_end_xx[ii,jj]) & 
        (fromstart_to_end_xx[ii,jj]) != 0) {
      leng <- sqrt(fromstart_to_end_xx[ii,jj]^2 + fromstart_to_end_yy[ii,jj]^2)
      points(-loc.start.x, -loc.start.y, cex=0.5, pch = 16)
      lines(c(-loc.start.x,-(loc.start.x + fromstart_to_end_xx[ii,jj]/(leng^0))),
            c(-loc.start.y, -(loc.start.y + fromstart_to_end_yy[ii,jj]/(leng^0))),
            lwd = 1)
      if(leng < 1.44){
        e <- c(loc.start.x, loc.start.y)
        equilibrium <- rbind(equilibrium, e)
        points(-loc.start.x, -loc.start.y, cex=0.5, pch = 16, col = "green")
        lines(c(-loc.start.x, -(loc.start.x + fromstart_to_end_xx[ii,jj]/(leng^0))),
              c(-loc.start.y, -(loc.start.y + fromstart_to_end_yy[ii,jj]/(leng^0))),
              lwd = 1, col = "green")
      }
    }
    
  }
}

points(-equilibrium[c(1, 2, 3, 5),], cex=0.5, pch = 16) # remove the uncontiguous points


## calculate equilibrium points (????)

equilibrium <- equilibrium[-c(1, 2 , 3, 5),] #remove the outlier
region <- as.data.frame(rbind(cbind(equilibrium[,1] - 2.5, equilibrium[,2] - 2.5),
                              cbind(equilibrium[,1] - 2.5, equilibrium[,2] + 2.5),
                              cbind(equilibrium[,1] + 2.5, equilibrium[,2] - 2.5),
                              cbind(equilibrium[,1] + 2.5, equilibrium[,2] + 2.5)))

region1 <- sqldf("select V1, V2
      from region 
      group by V1, V2
      having count(*) < 4")

points(region1,col="orange")
e_x <- c( 10, 10, 15, 15, 15, 20, 25, 25, 30, 30, 30, 30, 25, 25, 20, 20, 20, 15, 10)
e_y <- c(-25,-20,-20,-15,-10,-10,-10,-15,-15,-20,-25,-30,-30,-35,-35,-30,-25,-25,-25)
border <- cbind(e_x, e_y)
lines(border)

for(i in 1:nrow(immune.profiles.tsne)){
  x <- immune.profiles.tsne$tSNE.1[i]
  y <- immune.profiles.tsne$tSNE.2[i]
  immune.profiles.tsne$in_border[i] <- as.numeric(point.in.pol(e_x, e_y, x, y))
}

table(immune.profiles.tsne$in_border)

points(immune.profiles.tsne$tSNE.1[immune.profiles.tsne$in_border == 1],
       immune.profiles.tsne$tSNE.2[immune.profiles.tsne$in_border == 1])

data_in_border <- immune.profiles.tsne[immune.profiles.tsne$in_border == 1,]

summary(data_in_border$`lymphocyte percentage in nucleated cells`)
summary(data_in_border$`lymphocyte count`)
summary(data_in_border$`B cell count`)
summary(data_in_border$`T cell count`)
summary(data_in_border$`Treg cell count`)
summary(data_in_border$`CD4+ T cell count`)
summary(data_in_border$`CD8+ T cell percentage in lymphocytes`)
summary(data_in_border$`NK cell count`)



