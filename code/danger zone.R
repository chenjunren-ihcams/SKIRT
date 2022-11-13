library(readxl)
library(survival)
library(survminer)
library(pheatmap)
library(lubridate)
library(sqldf)


## Data preparation 

setwd("C:/Users/hp/Desktop/免疫重建-联合投影/upload files/")

SKIRT <- read_excel("./data/SKIRT.xlsx", guess_max = 1945)
type <- array('CV', nrow(SKIRT))
type[as.Date(SKIRT$`transplantation date`) >= '2019/01/01'] <- 'HO'
SKIRT <- cbind(SKIRT, type)
SKIRT <- SKIRT[SKIRT$center %in% c("IHCAMS-adult", "IHCAMS-child")]
patients <- SKIRT$ID
immune.profiles.tsne <- read_excel("./data/immune.profiles.tsne.xlsx", guess_max = 11150)


## Define functions

outline.func <- function(x, y, prob = 0.55, char = T){
  # input：coordinate of original points
  # output：coordinate of profile line
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

rotate <- function(x) t(apply(x, 2, rev))

point.in.pol <- function(pol.x, pol.y, point.x, point.y){
  #reference: https://blog.csdn.net/lynon/article/details/82015834
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



## Finding 'danger zone' of immune profiles
## by ten-fold cross-validation (CV) and hold-out validation (HO)


# Parameter setting

grid_x <- 121
grid_y <- 121
radius <- 6
case <- 10
alpha <- 1
p_threshold <- 0.01

validation <- c("CV", "HO")
pvalues <- c()

for(v in 1 : length(validation)){
  
  for (s in seq(0, 9, 3)){
    
    start <- s
    end <- s + 3
    
    x <- array(NA,nrow(SKIRT))
    y <- array(NA,nrow(SKIRT))
    patient_data <- c()
    
    
    # Data of specified time interval
    
    for (i in 1: length(patients)) {
      patient <- patients[i]
      
      II <- immune.profiles.tsne$ID == patient &
        immune.profiles.tsne$`days post-transplant` > start * 30 & immune.profiles.tsne$`days post-transplant` <= end * 30
      if (sum(II) > 1) {
        patient.data <- colMeans(immune.profiles.tsne[II, c("days post-transplant", "tSNE.1", "tSNE.2")])
        patient.data <- c(patient, as.numeric(patient.data))
      } else {
        if (sum(II) == 1) {
          patient.data <- immune.profiles.tsne[II, c("days post-transplant", "tSNE.1", "tSNE.2")]
          patient.data <- c(patient, as.numeric(patient.data))
        }
      }
      
      if (sum(II) > 0) {
        
        x[i] <- as.numeric(patient.data[3])
        y[i] <- as.numeric(patient.data[4])
        patient_data <- rbind(patient_data, patient.data)
        
      }
    }
    
    patient_data <- as.data.frame(patient_data,row.names = NULL)
    row.names(patient_data) <- NULL
    names(patient_data) <- names(immune.profiles.tsne)
    
    
    # Scatter plot of immune profiles
    
    J <- !is.na(SKIRT$`HLA-matched`) & 
      is.na(x) == FALSE 
    
    D <- J & SKIRT$death==1 
    
    
    JJ <- J & SKIRT$type == validation[v]
    DD <- JJ & SKIRT$death==1
    
    plot(x[JJ], y[JJ],
         cex=0.5, lwd=2, col=4, xlim=c(-60,60), ylim=c(-60,60), axes = FALSE,
         main=paste0("Immune profiles from days ",start*30 + 1,' to ',end*30,
                     '\n(cases = ',sum(JJ),', deaths (red) = ',sum(DD),")"),
         xlab='', ylab='')
    
    points(x[DD], y[DD], cex=0.5, lwd=2, col="red")
    lines(profile[,1], profile[,2], type = 'l')
    
    # Distance calculation
    
    x_seq <- seq(-60,60,length.out = grid_x)
    y_seq <- seq(-60,60,length.out = grid_y)
    delta_x_sq <- matrix(NA, nrow = nrow(SKIRT), ncol = grid_x)
    delta_y_sq <- matrix(NA, nrow = nrow(SKIRT), ncol = grid_y)
    dist_sq <- array(NA,c(grid_x,grid_y,nrow(SKIRT)))
    
    for (i in 1:grid_x) {
      delta_x_sq[,i] <- (x-x_seq[i])^2
    }
    
    for (j in 1:grid_y) {
      delta_y_sq[,j] <- (y-y_seq[j])^2
    }
    
    for (k in 1:nrow(SKIRT)) {
      for (i in 1:grid_x) {
        for(j in 1:grid_y) {
          dist_sq[i,j,k]<- delta_x_sq[k,i]+delta_y_sq[k,j]
        }
      }
    }
    
    if(dir.exists(paste0("./results/"))==0)
    {dir.create(paste0("./results/"))}
    
    
    # 10-fold cross-validation
    
    if(validation[v]=='CV'){
      
      set.seed(1)
      randomization <- sample(1:10, length(patients), replace=TRUE)
      beta_p_collect <- list()  
      
      for (kk in 1:max(randomization)) {
        
        beta <- matrix(NA, nrow= grid_y, ncol = grid_x)
        p <-  matrix(NA, nrow= grid_y, ncol = grid_x)
        beta_p <- matrix(NA, nrow= grid_y, ncol = grid_x)
        beta_p1_inflated <- matrix(NA, nrow= grid_y, ncol = grid_x)
        
        at.risk <- array(NA, length(patients))
        at.risk[J] <- 0
        
        for (i in 1:grid_x) {
          for(j in 1:grid_y) {
            
            group <- as.numeric(dist_sq[i,j,] <= (radius/1)^2)  ###debug：
            SKIRT.i.j <- data.frame(SKIRT[J & randomization!=kk & type=='CV',], group=group[J & randomization!=kk & type=='CV']) ###debug：
            
            
            if(sum(SKIRT.i.j$group==1) >= case){
              
              coxresult <- coxph(Surv(time,SKIRT.i.j$death) ~ group,
                                 data = SKIRT.i.j)
              
              beta[grid_y-j+1,i] <- coxresult$coefficients[names(coxresult$coefficients)=='group']
              p[grid_y-j+1,i] <- as.numeric(summary(coxresult)$coefficients["group", "Pr(>|z|)"] < p_threshold)
              
              
              if(p[grid_y-j+1,i] > 0){
                beta_p[grid_y-j+1,i] <- sign(beta[grid_y-j+1,i])
              }
              if(p[grid_y-j+1,i] == 0){
                beta_p[grid_y-j+1,i] <- NA   #将0改为NA就可以去掉
              }
              
              beta_p_collect[[kk]] <- beta_p
            }
          }
          cat('.')
        }
        
        
        if(TRUE){ # heatmap for each fold
          
          beta_p1 <- rotate(beta_p)
          
          for(i in 1 : nrow(beta_p1)){
            for (j in 1 : ncol(beta_p1)) {
              if( !is.na(beta_p1[i,j]) & beta_p1[i,j] == 1){
                beta_p1_inflated [i,j] <- 1
                for (m in 1 : nrow(beta_p1_inflated )) {
                  for (n in 1 : ncol(beta_p1_inflated )) {
                    if (((m-i)^2 + (n-j)^2) <= (radius/alpha)^2 & is.na(beta_p1_inflated [m,n]))
                    {beta_p1_inflated [m,n] <- 1}
                    
                  }
                }
              }
            }
          }
          
          for(i in -60:60){
            for (j in -60:60) {
              if (!is.na(beta_p1_inflated[i+61,j+61])){
                if (point.in.pol(profile[,1],profile[,2],i,j)==0){
                  beta_p1_inflated[i+61,j+61] <- NA
                }
              }
            }
          } 
          
          beta_p1_inflated[1,121] <- -1
          beta_p1_inflated[1,120] <- 1
          
          opar <- par(usr = c(-60, 60, -60, 60), xaxs = "i", yaxs = "i")
          image(x = seq(-60, 60, 1), 
                y = seq(-60, 60, 1),
                z = beta_p1_inflated, 
                col = hcl.colors(12, "Blue-Red 2", rev = F),
                xlim = c(-60, 60), 
                ylim=c(-60, 60), 
                ann = F, 
                xaxt = "n",
                yaxt = "n",
                bty = "n")
          lines(profile[,1], profile[,2], type = 'l', xlim = c(-60, 60), ylim=c(-60, 60))
          title(main = paste0('Danger zone during days ',start*30 + 1,' to ',end*30, '\n', 'fold = ',kk), cex.main = 1)
          par(opar)
          
        }
        
        
        cat('\n')
      }
      
      
      at.risk <- array(NA, length(patients))
      at.risk[J] <- 0
      for (kk in 1:max(randomization)) {
        for (i in 1:grid_x) {
          for(j in 1:grid_y) {
            if (!is.na(beta_p_collect[[kk]][grid_y-j+1,i])) {
              if (beta_p_collect[[kk]][grid_y-j+1,i] > 0) {
                at.risk[J][dist_sq[i,j,J] <= (radius/alpha)^2 & randomization[J]==kk & type[J]=='CV'] <- 1
              } 
              if (beta_p_collect[[kk]][grid_y-j+1,i] < 0) {
                at.risk[J][dist_sq[i,j,J] <= (radius/alpha)^2 & randomization[J]==kk & type[J]=='CV'] <- -1
              } 
            }
          }
        }
      }
    }
    
    
    
    
    # HO validation
    
    if(validation[v]=='HO'){
      
      beta <- matrix(NA, nrow= grid_y, ncol = grid_x)
      p <-  matrix(NA, nrow= grid_y, ncol = grid_x)
      beta_p <- matrix(NA, nrow= grid_y, ncol = grid_x)
      beta_p1_inflated <- matrix(NA, nrow= grid_y, ncol = grid_x)
      
      at.risk <- array(NA, length(patients))
      at.risk[J] <- 0
      
      for (i in 1:grid_x) {
        for(j in 1:grid_y) {
          
          group <- as.numeric(dist_sq[i,j,] <= (radius/1)^2)  
          SKIRT.i.j <- data.frame(SKIRT[J & type=='CV',], group=group[J & type=='CV'])   
          
          if(sum(SKIRT.i.j$group==1) >= case){
            
            coxresult <- coxph(Surv(time,SKIRT.i.j$death) ~ group,
                               data = SKIRT.i.j)
            
            beta[grid_y-j+1,i] <- coxresult$coefficients[names(coxresult$coefficients)=='group']
            p[grid_y-j+1,i] <- as.numeric(summary(coxresult)$coefficients["group", "Pr(>|z|)"] < p_threshold)
            
            
            if(p[grid_y-j+1,i]>0){
              beta_p[grid_y-j+1,i] <- sign(beta[grid_y-j+1,i])
            }
            if(p[grid_y-j+1,i]==0){
              beta_p[grid_y-j+1,i] <- NA   
            }
            
          }
        }
        cat('.')
      } 
      
      
      if(TRUE){ 
        
        beta_p1 <- rotate(beta_p)
        
        for(i in 1 : nrow(beta_p1)){
          for (j in 1 : ncol(beta_p1)) {
            if( !is.na(beta_p1[i,j]) & beta_p1[i,j] == 1){
              beta_p1_inflated [i,j] <- 1
              for (m in 1 : nrow(beta_p1_inflated )) {
                for (n in 1 : ncol(beta_p1_inflated )) {
                  if (((m-i)^2 + (n-j)^2) <= (radius/alpha)^2 & is.na(beta_p1_inflated [m,n]))
                  {beta_p1_inflated [m,n] <- 1}
                  
                }
              }
            }
          }
        }
        
        
        for(i in -60:60){
          for (j in -60:60) {
            if (!is.na(beta_p1_inflated[i+61,j+61])){
              if (point.in.pol(profile[,1],profile[,2],i,j)==0){
                beta_p1_inflated[i+61,j+61] <- NA
              }
            }
          }
        } 
        
        beta_p1_inflated[1,121] <- -1
        beta_p1_inflated[1,120] <- 1
        
        opar <- par(usr = c(-60, 60, -60, 60), xaxs = "i", yaxs = "i")
        
        image(x = seq(-60, 60, 1), 
              y = seq(-60, 60, 1),
              z = beta_p1_inflated, 
              col = hcl.colors(12, "Blue-Red 2", rev = F),
              xlim = c(-60, 60), 
              ylim =c(-60, 60), 
              ann = F, 
              xaxt = "n",
              yaxt = "n",
              xlab ='t-SNE 1',
              ylab ='t-SNE 2',
              bty = "n")
        lines(profile[,1], profile[,2], type = 'l', xlim = c(-60, 60), ylim=c(-60, 60))
        title(main = paste0('Danger zone during days ',start*30 + 1,' to ',end*30), cex.main = 1)
        
        par(opar)
      }
      
      
      write.csv(beta_p1_inflated, paste0("./results/Danger zone during days ",start*30 + 1, ' to ', end*30,'.csv'))
      
      for (i in 1:grid_x) {
        for(j in 1:grid_y) {
          if (!is.na(beta_p[grid_y-j+1,i])) {
            if (beta_p[grid_y-j+1,i]>0) {
              at.risk[J][dist_sq[i,j,J] <= (radius/alpha)^2 & type[J]=='HO'] <- 1
            } 
            if (beta_p[grid_y-j+1,i]<0) {
              at.risk[J][dist_sq[i,j,J] <= (radius/alpha)^2 & type[J]=='HO'] <- -1
            } 
          }
        }
      }
      
    }
    
    
    
    
    # Survival analysis
    
    datax <- cbind(SKIRT[,c("ID","center","death","time","type")], at.risk)[!is.na(at.risk),]
    datax$at.risk[datax$at.risk == -1] <- 0
    datax <- merge(datax, patient_data, by.x = 'ID', by.y = 'ID', all.x = TRUE, all.y = FALSE)
    
    if(validation[v] == 'HO'){
      datax <- datax[datax$type == 'HO',]
    }else{ 
      datax <- datax[datax$type == 'CV',]
    }
    
    surv_diff <- survdiff(Surv(time, death) ~ datax$at.risk, data = datax)
    p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
    pvalue <- c(validation[v], start*30 + 1, end*30, p.value)
    pvalues <- rbind(pvalues, pvalue)
    fit1 <- survfit(Surv(time, death) ~ datax$at.risk, data = datax)
    figure <- ggsurvplot(fit1,
                         title = paste0("Danger zone during days ",start*30 + 1," - ",end*30),
                         ggtheme = theme_bw(),
                         xlab="Time since transplantation (years) ", 
                         ylab="Overall survival",
                         xlim=c(0,4),
                         ylim=c(0.5,1),
                         break.time.by = 1,
                         risk.table = TRUE, 
                         tables.y.text=TRUE,
                         risk.table.fontsize=4,
                         censor.size=4,
                         palette=c('blue','red') ,
                         pval = T,
                         pval.method=T,
                         pval.coord=c(0,0.55),
                         pval.method.coord=c(0,0.6))
    print(figure)
    
    write.csv(datax,
              file = paste0("./results/Days from ",start*30 + 1, ' to ', end*30, " (", validation[v], ")", '.csv'),
              row.names = FALSE)
    
    cat('\n')
    
  }

}

pvalues <- as.data.frame(pvalues)
names(pvalues) <- c("validation", "start (days)", "end (days)", "p-value")
write.csv(pvalues, "./results/pvalues.csv", row.names = FALSE)















