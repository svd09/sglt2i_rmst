library (survRM2)
library(IPDfromKM)
library (metaRMST)

db=file.choose()
active <- read.table(db,
  sep="\t", header=FALSE)

db=file.choose()
comparator <- read.table(db,
                     sep="\t", header=FALSE)

db=file.choose()
no_at_risk<- read.table(db,header=FALSE)

active$V2 = 100-active$V2
comparator$V2 = 100-comparator$V2

### Get data from the sample dataset=======================

trisk <- no_at_risk$V1
nrisk_comparator <- no_at_risk$V2
nrisk_active <- no_at_risk$V3

### Estimate the IPD for the comparator therapy treatment group ====================
pre_comparator <- preprocess(dat=comparator, trisk=trisk,nrisk=nrisk_comparator,maxy=100)
est_comparator <- getIPD(prep=pre_comparator,armID=0,tot.events=NULL)

### Estimate the IPD for the active treatment group ====================
pre_active <- preprocess(dat=active, trisk=trisk,nrisk=nrisk_active,maxy=100)
est_active <- getIPD(prep=pre_active,armID=1,tot.events=NULL)

### survival report for two arms ===================
surv2 <- survreport(ipd1=est_comparator$IPD,ipd2=est_active$IPD,arms=2,
                    interval=8,s=c(0.75,0.5,0.25),showplots=TRUE)
print(surv2)

#combine both arms into one dataset
total_ipd <- rbind(est_comparator$IPD,est_active$IPD)

### COMBINE TRIALS
canvas <- total_ipd # then run load trial line 5->
credence <- total_ipd # then run load trial line 5->
empareg <- total_ipd # then run load trial line 5->
vertis <- total_ipd # then run load trial line 5->
declare <- total_ipd # then run load trial line 5->

## LABEL TRIALS
canvas$trialID =1
credence$trialID=2
empareg$trialID=3
vertis$trialID=4
declare$trialID=5

# Write IPD data
db2=file.choose()
write.csv (total_ipd, db2)

# combine trial ipd
overall <- rbind (canvas, empareg, vertis)

## RENAME variables for the function
names(overall) [1] <- "Time"
names(overall) [2] <- "Event"
names(overall) [3] <- "Arm"
names(overall) [4] <- "trialID"

# Write IPD data
db2=file.choose()
write.csv (overall, db2)

### metaRMST from saved ipd file
ff=file.choose()
CVDeath <- read.csv(ff, header=TRUE)
obj <- RMSTcurves(CVDeath, time_horizons=c(1:84), tmax=84, nboot=500,  MA_mvma = FALSE, MA_mvma_boot = FALSE,
                  MA_uni = FALSE, MA_uni_flex = TRUE)
RMSTplot(obj, xlim=c(0,84), ylim=c(-0.1,3), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="Cardiovascular Death", col = c("red", "blue", "green","orange"), trial_legend=TRUE, MA_legend = FALSE)
                                                                                                                                 

ff=file.choose()
MACE <- read.csv(ff, header=TRUE)
obj <- RMSTcurves(MACE, time_horizons=c(1:84), tmax=84, nboot=500,  MA_mvma = FALSE, MA_mvma_boot = FALSE,
                  MA_uni = FALSE, MA_uni_flex = TRUE)
RMSTplot(obj, xlim=c(0,84), ylim=c(-0.1,2), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="Major Adverse Cardiovascular Events (MI, Stroke, CV Death)", col = c("red", "green", "orange","purple"), trial_legend=TRUE, MA_legend = FALSE)

ff=file.choose()
HFH <- read.csv(ff, header=TRUE)
obj <- RMSTcurves(HFH, time_horizons=c(1:84), tmax=84, nboot=500,  MA_mvma = FALSE, MA_mvma_boot = FALSE,
                  MA_uni = FALSE, MA_uni_flex = TRUE)
RMSTplot(obj, xlim=c(0,84), ylim=c(-0.1,2), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="Heart Failure Hospitalizations", col = c("red", "green","orange"), trial_legend=TRUE, MA_legend = FALSE)

ff=file.choose()
ACM <- read.csv(ff, header=TRUE)
obj <- RMSTcurves(ACM, time_horizons=c(1:84), tmax=84, nboot=500,  MA_mvma = FALSE, MA_mvma_boot = FALSE,
                  MA_uni = FALSE, MA_uni_flex = TRUE)
RMSTplot(obj, xlim=c(0,84), ylim=c(-0.1,2), yby=0.5, ylab="RMSTD (mos)", xlab="Time (mos)", main="All-Cause Mortality", col = c("red", "blue", "green","orange"), trial_legend=TRUE, MA_legend = FALSE)
