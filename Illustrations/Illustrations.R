#### Load Libraries ####

library(FDboost)
library(fda)
library(refund)
source('wilcox.fda.R')

#### Resin Viscosity ####
data("viscosity")

Xv    <- matrix(viscosity$visAll, nrow = nrow(viscosity$visAll), ncol = ncol(viscosity$visAll))
fXv   <- fpca.face(Xv)
Yvis  <- fXv$Yhat

##### temperature of resin #####
TR  <- viscosity$T_A

wilcox.fda(Yvis ~ TR)

pdf('resin_temp.pdf')
col_TR <- ifelse(TR == 'low', rgb(t(col2rgb('royalblue3'))/255, alpha = 0.75),
                 rgb(t(col2rgb('navajowhite2'))/255, alpha = 0.75))
matplot(viscosity$timeAll, t(Yvis), type = 'n', xlab = 'Seconds',
        ylab = 'Viscosity', main = 'Viscosity of Resin\nby Temperature of Resion', 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25, lwd = 2)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
matplot(viscosity$timeAll, t(Yvis), add = TRUE, col = col_TR, type = 'l', lwd = 2)
dev.off()

##### temperature of curing agent #####
TC  <- viscosity$T_B

wilcox.fda(Yvis ~ TC)

pdf('cure_temp.pdf')
col_TC <- ifelse(TC == 'low', rgb(t(col2rgb('royalblue3'))/255, alpha = 0.75),
                 rgb(t(col2rgb('navajowhite2'))/255, alpha = 0.75))
matplot(viscosity$timeAll, t(Yvis), type = 'n', xlab = 'Seconds',
        ylab = 'Viscosity', main = 'Viscosity of Resin\nby Curing Agent', 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25, lwd = 2)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
matplot(viscosity$timeAll, t(Yvis), add = TRUE, col = col_TC, type = 'l', lwd = 2)
dev.off()

##### temperature of tools #####
TT  <- viscosity$T_C

wilcox.fda(Yvis ~ TT)

pdf('tool_temp.pdf')
col_TT <- ifelse(TT == 'low', rgb(t(col2rgb('royalblue3'))/255, alpha = 0.75),
                 rgb(t(col2rgb('navajowhite2'))/255, alpha = 0.75))
matplot(viscosity$timeAll, t(Yvis), type = 'n', xlab = 'Seconds',
        ylab = 'Viscosity', main = 'Viscosity of Resin\nby Temperature of Tools', 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25, lwd = 2)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
matplot(viscosity$timeAll, t(Yvis), add = TRUE, col = col_TT, type = 'l', lwd = 2)
dev.off()

##### rotational speed #####
RS  <- viscosity$rspeed

wilcox.fda(Yvis ~ RS)

pdf('rot_speed.pdf')
col_RS <- ifelse(RS == 'low', rgb(t(col2rgb('royalblue3'))/255, alpha = 0.75),
                 rgb(t(col2rgb('navajowhite2'))/255, alpha = 0.75))
matplot(viscosity$timeAll, t(Yvis), type = 'n', xlab = 'Seconds',
        ylab = 'Viscosity', main = 'Viscosity of Resin\nby Rotational Speed', 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25, lwd = 2)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
matplot(viscosity$timeAll, t(Yvis), add = TRUE, col = col_RS, type = 'l', lwd = 2)
dev.off()

##### mass flow #####
MF   <- viscosity$mflow

wilcox.fda(Yvis ~ MF)

pdf('mass_flow.pdf')
col_MF <- ifelse(MF == 'low', rgb(t(col2rgb('royalblue3'))/255, alpha = 0.75),
                 rgb(t(col2rgb('navajowhite2'))/255, alpha = 0.75))
matplot(viscosity$timeAll, t(Yvis), type = 'n', xlab = 'Seconds',
        ylab = 'Viscosity', main = 'Viscosity of Resin\nby Mass Flow', 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25, lwd = 2)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
matplot(viscosity$timeAll, t(Yvis), add = TRUE, col = col_MF, type = 'l', lwd = 2)
dev.off()



#### Canadian Weather ####

R   <- CanadianWeather$region

col_R   <- apply(as.matrix(R), 1, 
                 function(x){ switch(x,
                                     Arctic = rgb(t(col2rgb('royalblue3'))/255, alpha = 0.75),
                                     Atlantic = rgb(t(col2rgb('navajowhite2'))/255, alpha = 0.75),
                                     Continental = rgb(t(col2rgb('indianred'))/255, alpha = 0.75),
                                     Pacific = rgb(t(col2rgb('forestgreen'))/255, alpha = 0.75))})

##### temperature #####
XT    <- t(CanadianWeather$dailyAv[,,'Temperature.C'])
fXT   <- fpca.face(XT)
YT    <- fXT$Yhat

kruskal.fda(YT ~ R)

pdf('cw_temp.pdf')
matplot(1:365, t(YT), type = 'n', xlab = 'Days',
        ylab = 'Temperature (C)', main = 'Canadian Weather by Region', 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
matplot(1:365, t(YT), add = TRUE, col = col_R, type = 'l', lwd = 2)
dev.off()

##### precipitation #####
XP    <- t(CanadianWeather$dailyAv[,,'Precipitation.mm'])
fXP   <- fpca.face(XP)
YP    <- fXP$Yhat

kruskal.fda(YP ~ R)

pdf('cw_precip.pdf')
matplot(1:365, t(YP), type = 'n', xlab = 'Days',
        ylab = 'Precipitation (mm)', main = 'Canadian Weather by Region', 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
matplot(1:365, t(YP), add = TRUE, col = col_R, type = 'l', lwd = 2)
dev.off()


#### COVID Mobility ####

driveRequests   <- read.table('driving_requests.txt', header = TRUE)

driveRIA   <- driveRequests[which(driveRequests$state == 'Iowa'), -62]
driveRMN   <- driveRequests[which(driveRequests$state == 'Minnesota'), -62]

driveRMD   <- driveRequests[which(driveRequests$state == 'Maryland'), -62]
driveRVA   <- driveRequests[which(driveRequests$state == 'Virginia'), -62]
driveRWV   <- driveRequests[which(driveRequests$state == 'West Virginia'), -62]

driveRCO   <- driveRequests[which(driveRequests$state == 'Colorado'), -62]
driveRUT   <- driveRequests[which(driveRequests$state == 'Utah'), -62]


driveR_IAMN         <- rbind(driveRIA, driveRMN)
driveR_IAMN$state   <- c(rep('IA', nrow(driveRIA)), rep('MN', nrow(driveRMN)))

driveR_MVW          <- rbind(driveRMD, driveRVA, driveRWV)
driveR_MVW$state    <- c(rep('MD', nrow(driveRMD)), rep('VA', nrow(driveRVA)), rep('WV', nrow(driveRWV)))

driveR_COUT         <- rbind(driveRCO, driveRUT)
driveR_COUT$state   <- c(rep('CO', nrow(driveRCO)), rep('UT', nrow(driveRUT)))

ex_columns  <- paste('drives.', 1:61, sep = '')

IAMN      <- as.matrix(driveR_IAMN[,ex_columns])
fXIM      <- fpca.face(IAMN)
YIM       <- fXIM$Yhat
state_IM  <- driveR_IAMN$state

wilcox.fda(YIM ~ state_IM)

pdf('IAvMN.pdf')
col_IAMN <- ifelse(state_IM == 'IA', rgb(t(col2rgb('royalblue3'))/255, alpha = 0.75),
                   rgb(t(col2rgb('navajowhite2'))/255, alpha = 0.75))
matplot(1:61 - 31, t(YIM), type = 'n', xlab = 'Days Since Policy Implementation',
        ylab = 'Percent Change', main = 'Driving Requests\nIowa vs Minnesota',
        ylim = c(0, 250), 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25, lwd = 2)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
matplot(1:61 - 31, t(YIM), add = TRUE, col = col_IAMN, type = 'l', lwd = 2)
dev.off()

COUT      <- as.matrix(driveR_COUT[,ex_columns])
fXCU      <- fpca.face(COUT)
YCU       <- fXCU$Yhat
state_CU  <- driveR_COUT$state

wilcox.fda(YCU ~ state_CU)

pdf('COvUT.pdf')
col_COUT <- ifelse(state_CU == 'UT', rgb(t(col2rgb('royalblue3'))/255, alpha = 0.75),
                   rgb(t(col2rgb('navajowhite2'))/255, alpha = 0.75))
matplot(1:61 - 31, t(YCU), type = 'n', xlab = 'Days Since Policy Implementation',
        ylab = 'Percent Change', main = 'Driving Requests\nColordao vs Utah',
        ylim = c(0, 500), 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25, lwd = 2)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
matplot(1:61 - 31, t(YCU), add = TRUE, col = col_COUT, type = 'l', lwd = 2)
dev.off()


MVW         <- as.matrix(driveR_MVW[,ex_columns])
fXMVW       <- fpca.face(MVW)
YMVW        <- fXMVW$Yhat
state_MVW   <- driveR_MVW$state

kruskal.fda(YMVW ~ state_MVW)

pdf('MDvVAvWV.pdf')
col_MVW <- apply(as.matrix(state_MVW), 1, 
                 function(x){ switch(x,
                                     MD = rgb(t(col2rgb('royalblue3'))/255, alpha = 0.75),
                                     VA = rgb(t(col2rgb('indianred'))/255, alpha = 0.75),
                                     WV = rgb(t(col2rgb('navajowhite2'))/255, alpha = 0.75))})
matplot(1:61 - 31, t(YMVW), type = 'n', xlab = 'Days Since Policy Implementation',
        ylab = 'Percent Change', main = 'Driving Requests\nMaryland, Virgina, & West Virginia',
        ylim = c(0, 350), 
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25, lwd = 2)
abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
matplot(1:61 - 31, t(YMVW), add = TRUE, col = col_MVW, type = 'l', lwd = 2)
dev.off()
