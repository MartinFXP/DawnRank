#' patspeccuttoff
#' 
#' patspeccuttoff determines the patient specific cutoffs give
#' 
#' 
#' @param patient, the patient of interest
#' @param ms, a dataframe containing dawn rank output (specifically, the third
#' element of DawnRank output)
#' @param default, the default cutoff in case there are not enough common
#' drivers to make a clear conclusion
#' @return a data frame containing information regarding the genetic
#' information of a specific patient, and the cutoff itself
#' @export patspeccutoff
patspeccutoff <-
function(patient,ms,default=.95){
  require(maxstat)
  msp<-ms[ms[,2]==patient,]
  msp$significant<-0
  msp<-msp[order(msp[,3],decreasing=TRUE),]
  if((sum(msp$isCGC==1)>1)&(nrow(msp)>10)){
    maxstat.test(isCGC~PercentRank,data=msp,smethod="Wilcox",pmethod="HL")->trialcutoff
    msp$significant[msp$PercentRank>=trialcutoff$estimate]<-1
    cutoff<-trialcutoff$estimate
  }else{
    msp$significant[msp$PercentRank>=default]<-1
    cutoff<-default
  }
  if(sum(ms$significant)<1){
    msp$significant[1]<-1
    cutoff<-msp[1,3]
  }
  return(list(msp,cutoff))
}
