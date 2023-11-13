#MM_RSA_TKA

#Fitting Michaelis-Menten curves to radiostereometric analysis (RSA) mean maximum 
#total point motion (MTPM) migration data for total knee arthroplasty (TKA)

#Version 1. 2023-11-13 by Elise Laende

#Two packages for fitting non-linear models are used 
  #"drc" Analysis of Dose-Response Curves (https://cran.r-project.org/web/packages/drc/drc.pdf)
  #"nls" Non-linear regression

#Revision History
#Version 1. 2019-07-17 Script created for analysis of PTKA dataset

#Inputs

#Dataframe of time (follow-up visits in months) and group mean MTPM data in mm

#Example data format for single group:
  # FU_month        MTPM 
  # 0               0.00
  # 1.5             0.37
  # 3               0.40
  # 6               0.48
  # 12              0.49
  # 24              0.50

#Example data format for multiple groups identified by group ID:
#  ID    FU_month  MTPM
#  1      0.0      0.00
#  1      1.5      0.43
#  1      6.0      0.55
#  1     12.0      0.57
#  1     24.0      0.47
#  2      0.0      0.00
#  2      1.5      0.53
#  2      6.0      0.51
#  2     12.0      0.53
#  2     24.0      0.45
#  3      0.0      0.00
#  3      6.0      0.28
#  3     12.0      0.30
#  3     24.0      0.34
#  4      0.0      0.00
#  4      3.0      0.93
#  4     12.0      1.00
#  4     24.0      0.95


#We recommend time in months, but time can be in any consistent unit (i.e. weeks or months or years) 
    #*Be sure to not mix units* (use all months or all years)
    #*Include time 0 with 0 mm of MTPM migration
    #*Need to use a time variable, not follow-up number
    
  #Examples of acceptable times:
      # (0, 1.5, 3, 6, 12, 24)  in months  *recommended*
      # (0, 6, 12, 26, 52, 104)  in weeks
      # (0, 0.125, 0.25, 0.5, 1, 2) in years

  #Examples of INCORRECT time: (0, 6 weeks, 3 months, 6 months, 1 year, 2 years) *All visits must use the same measure of time* 
                            #  (X1, X2, X3, X4) *Visits must be in a unit of time*

#LIMITATIONS
  #*Only intended for group MEAN MTPM values (not individual cases)
  #*It may not be possible to fit curves on all data (models will not converge)

#Note: data need not be ordered

#Outputs

#MTPMemax: estimated maximum MTPM (mm) 
#K: k-value is half the time to the estimated maximum


#Install packages if necessary ---------------------------
#install.packages("drc")
#install.packages("plyr")
#install.packages("reshape2")

# Load libraries ---------------------------
library(drc)
library(plyr)


################################################################################
## PART 1 - Fit a single curve to a single group (mean MPTM) ###################
################################################################################


#Example data
FU_month<-c(0, 1.5, 3, 6, 12, 24) #time in months
MTPM<-c(0.00, 0.37,0.39,0.44,0.48, 0.49) #group mean MTPM values in mm

#generate data frame of time and mean MTPM migration data
df<-data.frame(FU_month, MTPM)



#OPTION 1: Use drc package ----------------------------------------------------
mm.drc <- drm(MTPM ~ FU_month, data = df, fct = MM.2())

#view model fitting results
summary(mm.drc)

#view model fit coefficients ("e" = MTPMemax, "d" = K)
mm.drc

#Plot the data points and the fitted curve
#(y limits for the plot are set to ensure the MTPMemax will be visible)
plot(mm.drc, log = '', pch = 17, main = "MM curve fitted with DRC", xlim=c(0,24), ylab="MTPM (mm)", xlab="Follow-up (months)")
abline(h=coef(mm.drc)[1], col="orange", lty=2)  #Add an orange dashed line at the MTPMemax
abline(v=coef(mm.drc)[2], col="green", lty=2)



#OPTION 2: use nls package ----------------------------------------------------

#fit model using nls (non-linear least squares) package
mm.nls<-nls(MTPM ~ MTPMemax*FU_month/(K+FU_month), data = df, start=list(K=max(df$FU_month/2),MTPMemax=max(df$MTPM)))

#view model fitting results
mm.nls

#view model fit coefficients (MTPMemax, K)
coef(mm.nls)
confint(mm.nls) #Can also calculate confidence intervals on coefficients

#plot raw data and fitted curve

attach(df)
plot(FU_month, MTPM, las=1, pch=16, main="MM curve fitted with NLS", xlim=c(0,24), ylab="MTPM (mm)", xlab="Follow-up (months)") #plot raw data as black dots
x <- seq(min(FU_month), max(FU_month), length=100) #calculate predicted model results
y <- predict(mm.nls, list(FU_month=x))
points(x, y, type='l', col='blue') #added fitted model (in blue)
#ylim(0, max(c(MTPM, coef(mm.nls)[2]))*1.1)
abline(h=coef(mm.nls)[2], col="orange", lty=2) #add horizontal line at MTPMemax
abline(v=coef(mm.nls)[1], col="green", lty=2) #add vertical line at K
detach()



################################################################################
## PART 2 - Use a function to fit curves to multiple groups  ###################
################################################################################

# #Function for fitting MM curves with error handling for non-convergence cases
# #Note: error message will display, but will not halt execution

# #Tries fitting first with DRC, then NLS
# #DRC works better with fewer visits
# #NLS works for some that don't converge (allows negative k values)
# Output includes MTPMemax, K, method used (drc or nls), number of datapoints including (0,0), number of follow-ups (excludes (0,0))
# Input: Dataframe with minimum of these 3 columns:
        #Group identifier (variable name can vary)
        #"FU_month" follow-up time in months
        #"MTPM" group mean MTPM in mm at associated follow-up month

MMRSA_DRC_NLS<-function(DataIn) {
  #check if fitting will generate an error
  possibleError <- tryCatch(
    drm(MTPM ~ FU_month, data = DataIn, fct = MM.2()),
    error=function(e) e
  )
  #if no error, fit curve, find coefficients
  if(!inherits(possibleError, "error")){
    mm<-drm(MTPM ~ FU_month, data = DataIn, fct = MM.2())
    mmout<-data.frame("MTPMemax"=mm$coefficients[[1]],"K"=mm$coefficients[[2]],"method"="drc", "n_datapoints"=length(DataIn$MTPM), "n_followups"=(length(DataIn$MTPM)-1))
    #if there is an error, try nls approach
  }else{
    #Try nls
    possibleError2 <- tryCatch(
      nls(MTPM ~ MTPMmax*FU_month/(KM+FU_month), data = DataIn, start=list(KM=max(DataIn$FU_month/2),MTPMmax=max(DataIn$MTPM))),
      error=function(e) e
    )
    #if no error, fit curve, find coefficients
    if(!inherits(possibleError2, "error")){
      mm.nls<-nls(MTPM ~ MTPMmax*FU_month/(KM+FU_month), data = DataIn, start=list(KM=max(DataIn$FU_month/2),MTPMmax=max(DataIn$MTPM)))
      mmout<-data.frame("MTPMemax"=coef(mm.nls)[[2]],"K"=coef(mm.nls)[[1]], "method"="nls", "n_datapoints"=length(DataIn$MTPM),"n_followups"=(length(DataIn$MTPM)-1))
      #if there is an error, output NA NA instead of coefficient
    }else{
      mmout<-data.frame("MTPMemax"=NA, "K"=NA, "method"=NA, "n_datapoints"=length(DataIn$MTPM), "n_followups"=(length(DataIn$MTPM)-1)) #if data cannot be fit, output "NA" values
    }
  }
}



#Sample data set with multiple groups with mean MTPM migrations
ID<-c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4) #group ID
FU_month<-c(0.0, 1.5 , 6.0 ,12.0 ,24.0,  0.0,  1.5 , 6.0, 12.0, 24.0,  0.0,  6.0, 12.0, 24.0 , 0.0 , 3.0 ,12.0 ,24.0) #follow-up times in months
MTPM<-c(0.00, 0.43, 0.55, 0.57, 0.47, 0.00, 0.53, 0.51, 0.53, 0.45, 0.00, 0.28, 0.30, 0.34, 0.00, 0.93, 1.00, 0.95) #group mean MTPM values in mm


#generate data frame of time and mean MTPM migration data
SampleData <-data.frame(ID, FU_month, MTPM)

#Use ddply to call MMRSA_DRC_NLS function and fit curve to each group (by "ID")
mm.curve.results<-ddply(SampleData, ~ID, MMRSA_DRC_NLS) 
#Note: error message indicates that for group 2, the drc method did not converge and model was fit with nls instead


#view results
mm.curve.results


 
# # #Note: Data need not be sorted
# 
# #time zero and MTPM zero at end....
# FU_month<-c(1.5, 3, 6, 12, 24,0) #time in months
# MTPM<-c(0.37,0.39,0.44,0.48, 0.49, 0) #group mean MTPM values in mm
# 
# #...give the same result as time zero and MTPM zero at the start.
# FU_month<-c(0, 1.5, 3, 6, 12, 24) #time in months
# MTPM<-c(0.00, 0.37,0.39,0.44,0.48, 0.49) #group mean MTPM values in mm


#References: https://davetang.org/muse/2013/05/17/fitting-a-michaelis-mentens-curve-using/
