rm(list=ls(all=TRUE))  
graphics.off()

###############
#reading the data
###############

data<- read.csv(file = "Data/GrowthRespPhotoData.csv")
#making a data frame of the desired data out of the read data
data<- data.frame(data$FinalID, data$OriginalTraitName ,data$OriginalTraitValue, data$ConTemp, data$ConSpecies)




##################
#organizing data
##################

#getting rid of the data that has less that data points for each thermal response
test_data <- data[data$data.FinalID %in% names(which(table(data$data.FinalID)>4)), ]
#getting rid of the NAs in the data
test_data<-na.omit(test_data)
#getting the smallest value
min_value<-min(test_data$data.OriginalTraitValue)
# make a 0.001 difference so that there wont be any 0 in the data set
min_value_plus<-( abs(min_value)+0.001)
# adding the min_value .By doing this , we make sure the there will not be any negative values or zeros in the dataframe
test_data$data.OriginalTraitValue<-test_data$data.OriginalTraitValue + min_value_plus
#making index use for pandas in python data analysis
test_data$panda_index <- seq(1:length(test_data$data.FinalID))


#new column of b0 eh el e? th tl 
#lmfit , you have the value for the slope and the intercept
  #subset
  #unique

#making all the temperature form celcius to kelvin to a new column
test_data$ConTemp_kelvin<-test_data$data.ConTemp+273.15


#arrange a new column of unique id for each data point for data manipulation
test_data$unique_id <- cumsum(!duplicated(test_data$data.FinalID))

#new column of log trait value
test_data$log_OriginalTraitValue<-log(test_data$data.OriginalTraitValue)

#new column of 1 over kt
k<-8.617 * 10 ^ -5
test_data$one_over_kt<- (1/k*(test_data$ConTemp_kelvin))

#bo and e , others constant
#subset each id
#log b0 should be slope


E <- c()
logbo<-c()
T_h<-c()
T_l<-c()
E_h<-c()
E_l<-c()
B0 <- c()
#x<- list()

#E l is the enzyme’s low-temperarure de-activation energy (eV) which controls the behavior of the enzyme (and the curve) at very low temperatures
#T l is the at which the enzyme is   50% low- temperature deactivated
#E h is the enzyme’s high-temperature de-activation energy (eV) which controls the behavior of the enzyme (and the curve) at very high temperatures,
#T h is the at which the enzyme is 50% high-temperature deactivated.
#E is the activation energy (eV) which controls the rise of the curve up to the peak in the “normal operating range” for the enzyme (below the peak of the curve and above T h ).

for (i in 1:length(unique(test_data$unique_id))){
  info<-subset(test_data,unique_id == i, select = c(log_OriginalTraitValue,one_over_kt))
  model<-lm(log_OriginalTraitValue~one_over_kt , data = info)
  logbo<-append(logbo,rep(abs(model$coefficients[[1]]),nrow(info)))
  E<-append(E, rep(abs(model$coefficients[[2]]), nrow(info)))
  E_h<-append(E_h, rep(5*model$coefficients[[2]], nrow(info)))
  E_l<-append(E_l, rep(model$coefficients[[2]]*0.5, nrow(info)))
  T_h<-append(T_h,rep(max(subset(test_data, unique_id == i, select = c(ConTemp_kelvin))),nrow(info)))
  T_l<-append(T_l,rep(min(subset(test_data, unique_id == i, select = c(ConTemp_kelvin))),nrow(info)))
  B0<-append(B0, rep(min(subset(test_data, unique_id == i, select = c(data.OriginalTraitValue))),nrow(info)))
 }

  #dllm<-min(test_data$data.OriginalTraitValue)

#for (i in test_data$data.OriginalTraitValue){
 # if (i < -200){
   # print (i)
  #}
#}

model<-summary(lm(log_OriginalTraitValue~one_over_kt , data = info))

test_data$logbo<-logbo
test_data$E<-E
test_data$E_h<-E_h
test_data$E_l<-E_l
test_data$T_h<-T_h
test_data$T_l<-T_l
test_data$B0 <- B0

#plot(info$one_over_kt,info$log_OriginalTraitValue)
#save the to a list, and then

FinalID <- test_data$data.FinalID
OriginalTraitValue <- test_data$data.OriginalTraitValue
unique_id <- test_data$unique_id
Temps <- test_data$ConTemp_kelvin
logB0 <- test_data$logbo
E <- test_data$E
E_h <- test_data$E_h
E_l <- test_data$E_l
T_h <-test_data$T_h
T_l <- test_data$T_l
one_over_kt<-test_data$one_over_kt
index_panda<-test_data$panda_index
B0 <- test_data$B0 



output<- data.frame(FinalID ,OriginalTraitValue,unique_id ,Temps, logB0,E,E_h , E_l,T_h,T_l,one_over_kt,index_panda,B0)





write.csv(output, file = "Data/output.csv", quote = FALSE, row.names = FALSE)


#iddd<-levels(thermal$data.FinalID)

