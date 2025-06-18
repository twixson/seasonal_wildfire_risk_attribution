library(tidyverse)
library(evd)
library(knitr)

setwd("C:/Users/troyw/Documents/CSU/MSPHD/Research/Time_Series_Extension_Spring2022/")


source("innovations_functions.R")

file.name <- "ERA5_data/ERA5_CA_FWI_05.csv"

# The following are the values that we must choose for forward and back transformation
past <- 1959:1978
present <- 2002:2021
gpd_quant <- 0.975 # quantile for gpd in original marginal estimation
window_size <- 14 # number of days on each side of the moving window for marginal estimation 
maxlag <- 20 # How many lags to estimate tpdf 
r0_quant <- 0.975 # high radial quantile for tpdf estimation 
set.seed(92432) # noise sequence seed
seq_length <- 50000 # length of noise sequence to generate MA(q)'s for fitting
max_theta <- maxlag # maximum reasonable order of for MA(q)
chosen_order <- 10 # what order MA(q) do we want to use 


#load the data 
data.temp <- read_csv(file.name, show_col_types = FALSE)

#restrict to months we are considering
data.temp <- filter(data.temp, MON %in% 6:10)

#restrict to years we are considering 
data_early <- filter(data.temp, YR %in% past)
data_late <- filter(data.temp, YR %in% present)

#Get rid of "bad" years (bad if it has fewer than 50 days or fewer than 5 observations above some value)
keep <- c()

#past climate
#0.7 quantile of past ISI used to determine if a year has enough moderately high observations 
# to considered a full year   
quant_good <- quantile(data_early$ISI, probs = 0.7) 
for(yr in unique(data_early$YR)){
  data_temp <- filter(data_early, YR == yr)
  if((length(data_temp$ISI) > 50) && (sum(data_temp$ISI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_early <- filter(data_early, YR %in% keep)

#current climate
keep <- c()

for(yr in unique(data_late$YR)){
  data_temp <- filter(data_late, YR == yr)
  if((length(data_temp$ISI) > 50) && (sum(data_temp$ISI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_late <- filter(data_late, YR %in% keep)

# Add a variable "days_since" to represent days since June 1
data_early$days_since <- as.numeric(data_early$DATE - as.Date(paste(data_early$YR, "-05-31", sep = "")))
data_late$days_since <- as.numeric(data_late$DATE - as.Date(paste(data_late$YR, "-05-31", sep = "")))

# save quantiles (i.e. sort the (2k+1) day moving window of observations for each day)
ecdf_early <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)
ecdf_late <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)

for(i in 1:153){
  # consider only days between day i-k and i+k and sort the values on those (2k+1) days
  temp_data_early <- sort(data_early$ISI[which(data_early$days_since %in% 
                                                 max(1, i-window_size):min(153, i+window_size))])
  temp_data_late <- sort(data_late$ISI[which(data_late$days_since %in% 
                                               max(1, i-window_size):min(153, i+window_size))])
  
  # save the sorted days within the moving window as the ecdf quantiles
  ecdf_early[1:length(temp_data_early), i] <- temp_data_early
  ecdf_late[1:length(temp_data_late), i] <- temp_data_late
}


# save the threshold values
mu_early <- mu_late <- rep(NA, 153)

for (i in 1:153){
  q_temp_early <- floor(gpd_quant*(1+length(na.omit(ecdf_early[,i]))))
  mu_early[i] <- ecdf_early[q_temp_early, i]
  
  q_temp_late <- floor(gpd_quant*(1+length(na.omit(ecdf_late[,i]))))
  mu_late[i] <- ecdf_late[q_temp_late, i]
}

# Compute and save the amount that observations exceed the threshold. 
exceedance_early <- exceedance_late <- c()
yrs_early <- yrs_late <- c()
days_early <- days_late <- c()

for (i in 1:153){
  indices_temp <- which(data_early$days_since == i) #take all of the ith days
  temp_data_early <- data_early$ISI[indices_temp] 
  temp_ex_early <- temp_data_early - mu_early[i] #compute exceedances
  exceedance_early <- c(exceedance_early, temp_ex_early) #append
  yrs_early <- c(yrs_early, data_early$YR[indices_temp])
  days_early <- c(days_early, rep(i, length(indices_temp)))
  
  indices_temp <- which(data_late$days_since == i)
  temp_data_late <- data_late$ISI[indices_temp]
  temp_ex_late <- temp_data_late - mu_late[i]
  exceedance_late <- c(exceedance_late, temp_ex_late)
  yrs_late <- c(yrs_late, data_late$YR[indices_temp])
  days_late <- c(days_late, rep(i, length(indices_temp)))
}

# Rearrange exceedances so that they are in the same order as the data
exceedance_early <- data.frame(exceedance_early, yrs_early, days_early) %>% arrange(yrs_early)
exceedance_late <- data.frame(exceedance_late, yrs_late, days_late) %>% arrange(yrs_late)

# fit the GPD to the exceedances and save the likelihood value
single_model_early <- fpot(exceedance_early$exceedance_early, threshold = 0)
single_model_late <- fpot(exceedance_late$exceedance_late, threshold = 0)

#save parameter values
scale_early <- unname(single_model_early$estimate[1])
shape_early <- unname(single_model_early$estimate[2])
scale_late <- unname(single_model_late$estimate[1])
shape_late <- unname(single_model_late$estimate[2])

### Transform to uniform
#column is year
ts_uniform_early <- ts_uniform_late <- matrix(NA, nrow = 153, ncol = 20) 

#early
for(i in 1:length(data_early$ISI)){
  j <- data_early$days_since[i]
  yr <- which(unique(data_early$YR) == data_early$YR[i])
  if(data_early$ISI[i] < mu_early[j]) {
    ts_uniform_early[j, yr] <- (which(ecdf_early[,j] == data_early$ISI[i]) /
      (length(na.omit(ecdf_early[,j]))+1))[1]
  } else {
    ts_uniform_early[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_early$ISI[i], 
                                                                loc = mu_early[j],
                                                                scale = single_model_early$estimate[1], 
                                                                shape = single_model_early$estimate[2])
  }
}

#late
for(i in 1:length(data_late$ISI)){
  j <- data_late$days_since[i]
  yr <- which(unique(data_late$YR) == data_late$YR[i])
  if(data_late$ISI[i] < mu_late[j]) {
    ts_uniform_late[j, yr] <- (which(ecdf_late[,j] == data_late$ISI[i]) /
      (length(na.omit(ecdf_late[,j]))+1))[1]
  } else {
    ts_uniform_late[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_late$ISI[i], 
                                                               loc = mu_late[j],
                                                               scale = single_model_late$estimate[1], 
                                                               shape = single_model_late$estimate[2])
  }
}


############
#Frechet with alpha=2
###########
RV_ts_early <- 1/sqrt(log(1/ts_uniform_early))
RV_ts_late <- 1/sqrt(log(1/ts_uniform_late))


# Subtract mean (Bias correction)
me <- mean(na.omit(c(RV_ts_early)))
ml <- mean(na.omit(c(RV_ts_late)))

RV_ts_early <- ifelse(RV_ts_early - me > 0,
                      RV_ts_early - me, 0)
RV_ts_late <- ifelse(RV_ts_late - ml > 0,
                     RV_ts_late - ml, 0)



# first we make a plot of the daily 0.95 quantile when considering all years
quantiles_95_per_day <- rep(NA, 153)
seasonality_quantile <- 0.975
num_days <- c(30, 31, 31, 30, 31)
num_days_in_prev_months <- 0
for (i in 6:10){
  num_days_this_month <- num_days[(i-5)]
  for (j in 1:num_days_this_month){
    quantiles_95_per_day[j + num_days_in_prev_months] <- quantile(filter(data.temp, MON == i,
                                                                         DAY == j)$ISI,
                                                                  seasonality_quantile, type = 1)
  }
  num_days_in_prev_months <- num_days_in_prev_months + num_days_this_month
}

k <- 15
q <- c(0.25, 0.5, 0.7, 0.85, 0.975)

data.temp$days_since <- rep(seq(0,152,1), length(unique(data.temp$YR)))
quant <- matrix(NA, nrow = (154-k), ncol=5)
means <- rep(NA, 154-k)

for (i in 0:(153-k-1)) {
  temp <- filter(data.temp, days_since %in% c(i:(i+k)))$ISI
  quant[(i+1),] <- quantile(temp, q, type = 1, na.rm = T)
  means[(i+1)] <- mean(temp)
}

plot(quantiles_95_per_day, type = "l",
     xlab = "Days Since June 1", ylab = "ISI")
lines(8:146, quant[,5], col = "blue")


quants_early <- quants_late <- matrix(NA, nrow = 153, ncol = 6)

for(i in 1:153){
  quants_early[i,1:5] <- quantile(na.omit(ecdf_early[,i]), probs = q)
  quants_early[i,6] <- mean(na.omit(ecdf_early[,i]))
  
  quants_late[i, 1:5] <- quantile(na.omit(ecdf_late[,i]), probs = q)
  quants_late[i, 6] <- mean(na.omit(ecdf_late[,i]))
}



# Compare early to late 0.95 quantile 
dates <- format(as.Date(as.Date('2020-06-01'):as.Date('2020-10-31'), 
                        origin = '2020-01-01'), format = '%m-%d')

plot(quants_late[,5], type = "l", 
     ylim = c(min(quants_early[,5]) -5, max(quants_late[,5]) +5), xlim = c(0,153),
     xlab = "", ylab = "ISI", xaxt="n")
axis(1, at = c(seq(1, 153, 14),153), labels = dates[c(seq(1, 153, 14),153)], las=2)
lines(quants_early[,5], lty = 2)
legend(x = 117, y = 72, c("Present", "", "", "", "Past"), bty = "n")


# gpd parameters
parameter <- c("scale", "shape")
past_paras <- c(scale_early, shape_early)
present_paras <- c(scale_late, shape_late)

gpd_paras <- data.frame("." = parameter, "past" = round(past_paras, 4), "present" = round(present_paras, 4))

kable(gpd_paras)



#########
#Estimate TPDF
#########

# compute lags
lags_early <- lags_late <- matrix(NA, nrow = (153*20), ncol = maxlag+1)

for(i in 1:20){
  indices <- c(((i-1)*153+1):(i*153))
  lags_early[indices, (maxlag+1)] <- RV_ts_early[,i] #save transformed data into a vector
  lags_late[indices, (maxlag+1)] <- RV_ts_late[,i]
  for(j in 1:maxlag){
    lags_early[indices, j] <- lag(RV_ts_early[,i], j) 
    #append each years lagged vectors onto previous years
    lags_late[indices, j] <- lag(RV_ts_late[,i], j)
  }
}

# compute radial components
rt_early <- rt_late <- matrix(NA, nrow <- (153*20), ncol = maxlag)

for(j in 1:maxlag){
  rt_early[,j] <- sqrt(lags_early[,j]^2 + lags_early[,(maxlag+1)]^2)
  
  rt_late[,j] <- sqrt(lags_late[,j]^2 + lags_late[,(maxlag+1)]^2)
}

# compute "high threshold of radial component" (0.95 quantile of radial components)
r0_early <- r0_late <- rep(NA, maxlag)

for (i in 1:maxlag){
  len <- length(na.omit(rt_early[,i]))
  ind <- floor(r0_quant*len)
  r0_early[i] <- sort(na.omit(rt_early[,i]))[ind]
  
  len <- length(na.omit(rt_late[,i]))
  ind <- floor(r0_quant*len)
  r0_late[i] <- sort(na.omit(rt_late[,i]))[ind]
}


# Compute tpdf
tpdf_early <- tpdf_late <- rep(NA, maxlag)

for (i in 1:maxlag) {
  wi_temp <- lags_early[,i]*lags_early[,(maxlag+1)]/(rt_early[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_early[,i]), 0, rt_early[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_early[i]) # how many above radial threshold
  tpdf_early[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_early[i])])
  
  wi_temp <- lags_late[,i]*lags_late[,(maxlag+1)]/(rt_late[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_late[,i]), 0, rt_late[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_late[i]) # how many above radial threshold
  tpdf_late[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_late[i])])
}

#add lag 0 to tpdf
tpdf_early <- c(1, tpdf_early)
tpdf_late <- c(1, tpdf_late)  


#plot
plot(0:maxlag, tpdf_early, type = "h", ylim = c(0,1), main = "ISI tpdf CA")
legend("topright", legend = c("Past", "Present"), bty = "n", text.col = c("black", "blue"))
lines(0:maxlag + 0.1, tpdf_late, type = "h", ylim = c(0,1), col = "blue")










###########33
###############3
###########3
rm(list = ls())






source("innovations_functions.R")

file.name <- "ERA5_data/ERA5_CA_FWI_05.csv"

# The following are the values that we must choose for forward and back transformation
past <- 1959:1978
present <- 2002:2021
gpd_quant <- 0.975 # quantile for gpd in original marginal estimation
window_size <- 14 # number of days on each side of the moving window for marginal estimation 
maxlag <- 50 # How many lags to estimate tpdf 
r0_quant <- 0.975 # high radial quantile for tpdf estimation 
set.seed(92432) # noise sequence seed
seq_length <- 50000 # length of noise sequence to generate MA(q)'s for fitting
max_theta <- maxlag # maximum reasonable order of for MA(q)
chosen_order <- 10 # what order MA(q) do we want to use 




#load the data 
data.temp <- read_csv(file.name, show_col_types = FALSE)

#restrict to months we are considering
data.temp <- filter(data.temp, MON %in% 6:10)

#restrict to years we are considering 
data_early <- filter(data.temp, YR %in% past)
data_late <- filter(data.temp, YR %in% present)

#Get rid of "bad" years (bad if it has fewer than 50 days or fewer than 5 observations above some value)
keep <- c()

#past climate
#0.7 quantile of past BUI used to determine if a year has enough moderately high observations 
# to considered a full year   
quant_good <- quantile(data_early$BUI, probs = 0.7) 
for(yr in unique(data_early$YR)){
  data_temp <- filter(data_early, YR == yr)
  if((length(data_temp$BUI) > 50) && (sum(data_temp$BUI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_early <- filter(data_early, YR %in% keep)

#current climate
keep <- c()

for(yr in unique(data_late$YR)){
  data_temp <- filter(data_late, YR == yr)
  if((length(data_temp$BUI) > 50) && (sum(data_temp$BUI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_late <- filter(data_late, YR %in% keep)

# Add a variable "days_since" to represent days since June 1
data_early$days_since <- as.numeric(data_early$DATE - as.Date(paste(data_early$YR, "-05-31", sep = "")))
data_late$days_since <- as.numeric(data_late$DATE - as.Date(paste(data_late$YR, "-05-31", sep = "")))

# save quantiles (i.e. sort the (2k+1) day moving window of observations for each day)
ecdf_early <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)
ecdf_late <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)

for(i in 1:153){
  # consider only days between day i-k and i+k and sort the values on those (2k+1) days
  temp_data_early <- sort(data_early$BUI[which(data_early$days_since %in% 
                                                 max(1, i-window_size):min(153, i+window_size))])
  temp_data_late <- sort(data_late$BUI[which(data_late$days_since %in% 
                                               max(1, i-window_size):min(153, i+window_size))])
  
  # save the sorted days within the moving window as the ecdf quantiles
  ecdf_early[1:length(temp_data_early), i] <- temp_data_early
  ecdf_late[1:length(temp_data_late), i] <- temp_data_late
}


# save the threshold values
mu_early <- mu_late <- rep(NA, 153)

for (i in 1:153){
  q_temp_early <- floor(gpd_quant*(1+length(na.omit(ecdf_early[,i]))))
  mu_early[i] <- ecdf_early[q_temp_early, i]
  
  q_temp_late <- floor(gpd_quant*(1+length(na.omit(ecdf_late[,i]))))
  mu_late[i] <- ecdf_late[q_temp_late, i]
}

# Compute and save the amount that observations exceed the threshold. 
exceedance_early <- exceedance_late <- c()
yrs_early <- yrs_late <- c()
days_early <- days_late <- c()

for (i in 1:153){
  indices_temp <- which(data_early$days_since == i) #take all of the ith days
  temp_data_early <- data_early$BUI[indices_temp] 
  temp_ex_early <- temp_data_early - mu_early[i] #compute exceedances
  exceedance_early <- c(exceedance_early, temp_ex_early) #append
  yrs_early <- c(yrs_early, data_early$YR[indices_temp])
  days_early <- c(days_early, rep(i, length(indices_temp)))
  
  indices_temp <- which(data_late$days_since == i)
  temp_data_late <- data_late$BUI[indices_temp]
  temp_ex_late <- temp_data_late - mu_late[i]
  exceedance_late <- c(exceedance_late, temp_ex_late)
  yrs_late <- c(yrs_late, data_late$YR[indices_temp])
  days_late <- c(days_late, rep(i, length(indices_temp)))
}

# Rearrange exceedances so that they are in the same order as the data
exceedance_early <- data.frame(exceedance_early, yrs_early, days_early) %>% arrange(yrs_early)
exceedance_late <- data.frame(exceedance_late, yrs_late, days_late) %>% arrange(yrs_late)

# fit the GPD to the exceedances and save the likelihood value
single_model_early <- fpot(exceedance_early$exceedance_early, threshold = 0)
single_model_late <- fpot(exceedance_late$exceedance_late, threshold = 0)

#save parameter values
scale_early <- unname(single_model_early$estimate[1])
shape_early <- unname(single_model_early$estimate[2])
scale_late <- unname(single_model_late$estimate[1])
shape_late <- unname(single_model_late$estimate[2])

### Transform to uniform
#column is year
ts_uniform_early <- ts_uniform_late <- matrix(NA, nrow = 153, ncol = 20) 

#early
for(i in 1:length(data_early$BUI)){
  j <- data_early$days_since[i]
  yr <- which(unique(data_early$YR) == data_early$YR[i])
  if(data_early$BUI[i] < mu_early[j]) {
    ts_uniform_early[j, yr] <- (which(ecdf_early[,j] == data_early$BUI[i]) /
      (length(na.omit(ecdf_early[,j]))+1))[1]
  } else {
    ts_uniform_early[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_early$BUI[i], 
                                                                loc = mu_early[j],
                                                                scale = single_model_early$estimate[1], 
                                                                shape = single_model_early$estimate[2])
  }
}

#late
for(i in 1:length(data_late$BUI)){
  j <- data_late$days_since[i]
  yr <- which(unique(data_late$YR) == data_late$YR[i])
  if(data_late$BUI[i] < mu_late[j]) {
    ts_uniform_late[j, yr] <- (which(ecdf_late[,j] == data_late$BUI[i]) /
      (length(na.omit(ecdf_late[,j]))+1))[1]
  } else {
    ts_uniform_late[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_late$BUI[i], 
                                                               loc = mu_late[j],
                                                               scale = single_model_late$estimate[1], 
                                                               shape = single_model_late$estimate[2])
  }
}


############
#Frechet with alpha=2
###########
RV_ts_early <- 1/sqrt(log(1/ts_uniform_early))
RV_ts_late <- 1/sqrt(log(1/ts_uniform_late))


# Subtract mean (Bias correction)
me <- mean(na.omit(c(RV_ts_early)))
ml <- mean(na.omit(c(RV_ts_late)))

RV_ts_early <- ifelse(RV_ts_early - me > 0,
                      RV_ts_early - me, 0)
RV_ts_late <- ifelse(RV_ts_late - ml > 0,
                     RV_ts_late - ml, 0)



# first we make a plot of the daily 0.95 quantile when considering all years
quantiles_95_per_day <- rep(NA, 153)
seasonality_quantile <- 0.975
num_days <- c(30, 31, 31, 30, 31)
num_days_in_prev_months <- 0
for (i in 6:10){
  num_days_this_month <- num_days[(i-5)]
  for (j in 1:num_days_this_month){
    quantiles_95_per_day[j + num_days_in_prev_months] <- quantile(filter(data.temp, MON == i,
                                                                         DAY == j)$BUI,
                                                                  seasonality_quantile, type = 1)
  }
  num_days_in_prev_months <- num_days_in_prev_months + num_days_this_month
}

k <- 15
q <- c(0.25, 0.5, 0.7, 0.85, 0.975)

data.temp$days_since <- rep(seq(0,152,1), length(unique(data.temp$YR)))
quant <- matrix(NA, nrow = (154-k), ncol=5)
means <- rep(NA, 154-k)

for (i in 0:(153-k-1)) {
  temp <- filter(data.temp, days_since %in% c(i:(i+k)))$BUI
  quant[(i+1),] <- quantile(temp, q, type = 1, na.rm = T)
  means[(i+1)] <- mean(temp)
}

plot(quantiles_95_per_day, type = "l",
     xlab = "Days Since June 1", ylab = "BUI")
lines(8:146, quant[,5], col = "blue")


quants_early <- quants_late <- matrix(NA, nrow = 153, ncol = 6)

for(i in 1:153){
  quants_early[i,1:5] <- quantile(na.omit(ecdf_early[,i]), probs = q)
  quants_early[i,6] <- mean(na.omit(ecdf_early[,i]))
  
  quants_late[i, 1:5] <- quantile(na.omit(ecdf_late[,i]), probs = q)
  quants_late[i, 6] <- mean(na.omit(ecdf_late[,i]))
}


# Compare early to late 0.95 quantile 
dates <- format(as.Date(as.Date('2020-06-01'):as.Date('2020-10-31'), 
                        origin = '2020-01-01'), format = '%m-%d')

plot(quants_late[,5], type = "l", 
     ylim = c(min(quants_early[,5]) -5, max(quants_late[,5]) +5), xlim = c(0,153),
     xlab = "", ylab = "BUI", xaxt="n")
axis(1, at = c(seq(1, 153, 14),153), labels = dates[c(seq(1, 153, 14),153)], las=2)
lines(quants_early[,5], lty = 2)
legend(x = 117, y = 72, c("Present", "", "", "", "Past"), bty = "n")


# gpd parameters
parameter <- c("scale", "shape")
past_paras <- c(scale_early, shape_early)
present_paras <- c(scale_late, shape_late)

gpd_paras <- data.frame("." = parameter, "past" = round(past_paras, 4), "present" = round(present_paras, 4))

kable(gpd_paras)



#########
#Estimate TPDF
#########

# compute lags
lags_early <- lags_late <- matrix(NA, nrow = (153*20), ncol = maxlag+1)

for(i in 1:20){
  indices <- c(((i-1)*153+1):(i*153))
  lags_early[indices, (maxlag+1)] <- RV_ts_early[,i] #save transformed data into a vector
  lags_late[indices, (maxlag+1)] <- RV_ts_late[,i]
  for(j in 1:maxlag){
    lags_early[indices, j] <- lag(RV_ts_early[,i], j) 
    #append each years lagged vectors onto previous years
    lags_late[indices, j] <- lag(RV_ts_late[,i], j)
  }
}

# compute radial components
rt_early <- rt_late <- matrix(NA, nrow <- (153*20), ncol = maxlag)

for(j in 1:maxlag){
  rt_early[,j] <- sqrt(lags_early[,j]^2 + lags_early[,(maxlag+1)]^2)
  
  rt_late[,j] <- sqrt(lags_late[,j]^2 + lags_late[,(maxlag+1)]^2)
}

# compute "high threshold of radial component" (0.95 quantile of radial components)
r0_early <- r0_late <- rep(NA, maxlag)

for (i in 1:maxlag){
  len <- length(na.omit(rt_early[,i]))
  ind <- floor(r0_quant*len)
  r0_early[i] <- sort(na.omit(rt_early[,i]))[ind]
  
  len <- length(na.omit(rt_late[,i]))
  ind <- floor(r0_quant*len)
  r0_late[i] <- sort(na.omit(rt_late[,i]))[ind]
}


# Compute tpdf
tpdf_early <- tpdf_late <- rep(NA, maxlag)

for (i in 1:maxlag) {
  wi_temp <- lags_early[,i]*lags_early[,(maxlag+1)]/(rt_early[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_early[,i]), 0, rt_early[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_early[i]) # how many above radial threshold
  tpdf_early[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_early[i])])
  
  wi_temp <- lags_late[,i]*lags_late[,(maxlag+1)]/(rt_late[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_late[,i]), 0, rt_late[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_late[i]) # how many above radial threshold
  tpdf_late[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_late[i])])
}

#add lag 0 to tpdf
tpdf_early <- c(1, tpdf_early)
tpdf_late <- c(1, tpdf_late)  


#plot
plot(0:maxlag, tpdf_early, type = "h", ylim = c(0,1), main = "BUI tpdf CA")
legend("topright", legend = c("Past", "Present"), bty = "n", text.col = c("black", "blue"))
lines(0:maxlag + 0.1, tpdf_late, type = "h", ylim = c(0,1), col = "blue")























rm(list = ls())
###############
########################################3
# COLORADO 
#######################################
##############

















file.name <- "ERA5_data/ERA5_CO_FWI_07_15.csv"

# The following are the values that we must choose for forward and back transformation
past <- 1959:1978
present <- 2002:2021
gpd_quant <- 0.975 # quantile for gpd in original marginal estimation
window_size <- 14 # number of days on each side of the moving window for marginal estimation 
maxlag <- 20 # How many lags to estimate tpdf 
r0_quant <- 0.975 # high radial quantile for tpdf estimation 
set.seed(92432) # noise sequence seed
seq_length <- 50000 # length of noise sequence to generate MA(q)'s for fitting
max_theta <- maxlag # maximum reasonable order of for MA(q)
chosen_order <- 10 # what order MA(q) do we want to use 


#load the data 
data.temp <- read_csv(file.name, show_col_types = FALSE)

#restrict to months we are considering
data.temp <- filter(data.temp, MON %in% 6:10)

#restrict to years we are considering 
data_early <- filter(data.temp, YR %in% past)
data_late <- filter(data.temp, YR %in% present)

#Get rid of "bad" years (bad if it has fewer than 50 days or fewer than 5 observations above some value)
keep <- c()

#past climate
#0.7 quantile of past ISI used to determine if a year has enough moderately high observations 
# to considered a full year   
quant_good <- quantile(data_early$ISI, probs = 0.7) 
for(yr in unique(data_early$YR)){
  data_temp <- filter(data_early, YR == yr)
  if((length(data_temp$ISI) > 50) && (sum(data_temp$ISI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_early <- filter(data_early, YR %in% keep)

#current climate
keep <- c()

for(yr in unique(data_late$YR)){
  data_temp <- filter(data_late, YR == yr)
  if((length(data_temp$ISI) > 50) && (sum(data_temp$ISI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_late <- filter(data_late, YR %in% keep)

# Add a variable "days_since" to represent days since June 1
data_early$days_since <- as.numeric(data_early$DATE - as.Date(paste(data_early$YR, "-05-31", sep = "")))
data_late$days_since <- as.numeric(data_late$DATE - as.Date(paste(data_late$YR, "-05-31", sep = "")))

# save quantiles (i.e. sort the (2k+1) day moving window of observations for each day)
ecdf_early <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)
ecdf_late <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)

for(i in 1:153){
  # consider only days between day i-k and i+k and sort the values on those (2k+1) days
  temp_data_early <- sort(data_early$ISI[which(data_early$days_since %in% 
                                                 max(1, i-window_size):min(153, i+window_size))])
  temp_data_late <- sort(data_late$ISI[which(data_late$days_since %in% 
                                               max(1, i-window_size):min(153, i+window_size))])
  
  # save the sorted days within the moving window as the ecdf quantiles
  ecdf_early[1:length(temp_data_early), i] <- temp_data_early
  ecdf_late[1:length(temp_data_late), i] <- temp_data_late
}


# save the threshold values
mu_early <- mu_late <- rep(NA, 153)

for (i in 1:153){
  q_temp_early <- floor(gpd_quant*(1+length(na.omit(ecdf_early[,i]))))
  mu_early[i] <- ecdf_early[q_temp_early, i]
  
  q_temp_late <- floor(gpd_quant*(1+length(na.omit(ecdf_late[,i]))))
  mu_late[i] <- ecdf_late[q_temp_late, i]
}

# Compute and save the amount that observations exceed the threshold. 
exceedance_early <- exceedance_late <- c()
yrs_early <- yrs_late <- c()
days_early <- days_late <- c()

for (i in 1:153){
  indices_temp <- which(data_early$days_since == i) #take all of the ith days
  temp_data_early <- data_early$ISI[indices_temp] 
  temp_ex_early <- temp_data_early - mu_early[i] #compute exceedances
  exceedance_early <- c(exceedance_early, temp_ex_early) #append
  yrs_early <- c(yrs_early, data_early$YR[indices_temp])
  days_early <- c(days_early, rep(i, length(indices_temp)))
  
  indices_temp <- which(data_late$days_since == i)
  temp_data_late <- data_late$ISI[indices_temp]
  temp_ex_late <- temp_data_late - mu_late[i]
  exceedance_late <- c(exceedance_late, temp_ex_late)
  yrs_late <- c(yrs_late, data_late$YR[indices_temp])
  days_late <- c(days_late, rep(i, length(indices_temp)))
}

# Rearrange exceedances so that they are in the same order as the data
exceedance_early <- data.frame(exceedance_early, yrs_early, days_early) %>% arrange(yrs_early)
exceedance_late <- data.frame(exceedance_late, yrs_late, days_late) %>% arrange(yrs_late)

# fit the GPD to the exceedances and save the likelihood value
single_model_early <- fpot(exceedance_early$exceedance_early, threshold = 0)
single_model_late <- fpot(exceedance_late$exceedance_late, threshold = 0)

#save parameter values
scale_early <- unname(single_model_early$estimate[1])
shape_early <- unname(single_model_early$estimate[2])
scale_late <- unname(single_model_late$estimate[1])
shape_late <- unname(single_model_late$estimate[2])

### Transform to uniform
#column is year
ts_uniform_early <- ts_uniform_late <- matrix(NA, nrow = 153, ncol = 20) 

#early
for(i in 1:length(data_early$ISI)){
  j <- data_early$days_since[i]
  yr <- which(unique(data_early$YR) == data_early$YR[i])
  if(data_early$ISI[i] < mu_early[j]) {
    ts_uniform_early[j, yr] <- (which(ecdf_early[,j] == data_early$ISI[i]) /
      (length(na.omit(ecdf_early[,j]))+1))[1]
  } else {
    ts_uniform_early[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_early$ISI[i], 
                                                                loc = mu_early[j],
                                                                scale = single_model_early$estimate[1], 
                                                                shape = single_model_early$estimate[2])
  }
}

#late
for(i in 1:length(data_late$ISI)){
  j <- data_late$days_since[i]
  yr <- which(unique(data_late$YR) == data_late$YR[i])
  if(data_late$ISI[i] < mu_late[j]) {
    ts_uniform_late[j, yr] <- (which(ecdf_late[,j] == data_late$ISI[i]) /
      (length(na.omit(ecdf_late[,j]))+1))[1]
  } else {
    ts_uniform_late[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_late$ISI[i], 
                                                               loc = mu_late[j],
                                                               scale = single_model_late$estimate[1], 
                                                               shape = single_model_late$estimate[2])
  }
}


############
#Frechet with alpha=2
###########
RV_ts_early <- 1/sqrt(log(1/ts_uniform_early))
RV_ts_late <- 1/sqrt(log(1/ts_uniform_late))


# Subtract mean (Bias correction)
me <- mean(na.omit(c(RV_ts_early)))
ml <- mean(na.omit(c(RV_ts_late)))

RV_ts_early <- ifelse(RV_ts_early - me > 0,
                      RV_ts_early - me, 0)
RV_ts_late <- ifelse(RV_ts_late - ml > 0,
                     RV_ts_late - ml, 0)



# first we make a plot of the daily 0.95 quantile when considering all years
quantiles_95_per_day <- rep(NA, 153)
seasonality_quantile <- 0.975
num_days <- c(30, 31, 31, 30, 31)
num_days_in_prev_months <- 0
for (i in 6:10){
  num_days_this_month <- num_days[(i-5)]
  for (j in 1:num_days_this_month){
    quantiles_95_per_day[j + num_days_in_prev_months] <- quantile(filter(data.temp, MON == i,
                                                                         DAY == j)$ISI,
                                                                  seasonality_quantile, type = 1)
  }
  num_days_in_prev_months <- num_days_in_prev_months + num_days_this_month
}

k <- 15
q <- c(0.25, 0.5, 0.7, 0.85, 0.975)

data.temp$days_since <- rep(seq(0,152,1), length(unique(data.temp$YR)))
quant <- matrix(NA, nrow = (154-k), ncol=5)
means <- rep(NA, 154-k)

for (i in 0:(153-k-1)) {
  temp <- filter(data.temp, days_since %in% c(i:(i+k)))$ISI
  quant[(i+1),] <- quantile(temp, q, type = 1, na.rm = T)
  means[(i+1)] <- mean(temp)
}

plot(quantiles_95_per_day, type = "l",
     xlab = "Days Since June 1", ylab = "ISI")
lines(8:146, quant[,5], col = "blue")


quants_early <- quants_late <- matrix(NA, nrow = 153, ncol = 6)

for(i in 1:153){
  quants_early[i,1:5] <- quantile(na.omit(ecdf_early[,i]), probs = q)
  quants_early[i,6] <- mean(na.omit(ecdf_early[,i]))
  
  quants_late[i, 1:5] <- quantile(na.omit(ecdf_late[,i]), probs = q)
  quants_late[i, 6] <- mean(na.omit(ecdf_late[,i]))
}



# Compare early to late 0.95 quantile 
dates <- format(as.Date(as.Date('2020-06-01'):as.Date('2020-10-31'), 
                        origin = '2020-01-01'), format = '%m-%d')

plot(quants_late[,5], type = "l", 
     ylim = c(min(quants_early[,5]) -5, max(quants_late[,5]) +5), xlim = c(0,153),
     xlab = "", ylab = "ISI", xaxt="n")
axis(1, at = c(seq(1, 153, 14),153), labels = dates[c(seq(1, 153, 14),153)], las=2)
lines(quants_early[,5], lty = 2)
legend(x = 117, y = 72, c("Present", "", "", "", "Past"), bty = "n")


# gpd parameters
parameter <- c("scale", "shape")
past_paras <- c(scale_early, shape_early)
present_paras <- c(scale_late, shape_late)

gpd_paras <- data.frame("." = parameter, "past" = round(past_paras, 4), "present" = round(present_paras, 4))

kable(gpd_paras)



#########
#Estimate TPDF
#########

# compute lags
lags_early <- lags_late <- matrix(NA, nrow = (153*20), ncol = maxlag+1)

for(i in 1:20){
  indices <- c(((i-1)*153+1):(i*153))
  lags_early[indices, (maxlag+1)] <- RV_ts_early[,i] #save transformed data into a vector
  lags_late[indices, (maxlag+1)] <- RV_ts_late[,i]
  for(j in 1:maxlag){
    lags_early[indices, j] <- lag(RV_ts_early[,i], j) 
    #append each years lagged vectors onto previous years
    lags_late[indices, j] <- lag(RV_ts_late[,i], j)
  }
}

# compute radial components
rt_early <- rt_late <- matrix(NA, nrow <- (153*20), ncol = maxlag)

for(j in 1:maxlag){
  rt_early[,j] <- sqrt(lags_early[,j]^2 + lags_early[,(maxlag+1)]^2)
  
  rt_late[,j] <- sqrt(lags_late[,j]^2 + lags_late[,(maxlag+1)]^2)
}

# compute "high threshold of radial component" (0.95 quantile of radial components)
r0_early <- r0_late <- rep(NA, maxlag)

for (i in 1:maxlag){
  len <- length(na.omit(rt_early[,i]))
  ind <- floor(r0_quant*len)
  r0_early[i] <- sort(na.omit(rt_early[,i]))[ind]
  
  len <- length(na.omit(rt_late[,i]))
  ind <- floor(r0_quant*len)
  r0_late[i] <- sort(na.omit(rt_late[,i]))[ind]
}


# Compute tpdf
tpdf_early <- tpdf_late <- rep(NA, maxlag)

for (i in 1:maxlag) {
  wi_temp <- lags_early[,i]*lags_early[,(maxlag+1)]/(rt_early[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_early[,i]), 0, rt_early[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_early[i]) # how many above radial threshold
  tpdf_early[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_early[i])])
  
  wi_temp <- lags_late[,i]*lags_late[,(maxlag+1)]/(rt_late[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_late[,i]), 0, rt_late[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_late[i]) # how many above radial threshold
  tpdf_late[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_late[i])])
}

#add lag 0 to tpdf
tpdf_early <- c(1, tpdf_early)
tpdf_late <- c(1, tpdf_late)  


#plot
plot(0:maxlag, tpdf_early, type = "h", ylim = c(0,1), main = "ISI tpdf CO")
legend("topright", legend = c("Past", "Present"), bty = "n", text.col = c("black", "blue"))
lines(0:maxlag + 0.1, tpdf_late, type = "h", ylim = c(0,1), col = "blue")










###########33
###############3
###########3
rm(list = ls())



file.name <- "ERA5_data/ERA5_CO_FWI_07_15.csv"

# The following are the values that we must choose for forward and back transformation
past <- 1959:1978
present <- 2002:2021
gpd_quant <- 0.975 # quantile for gpd in original marginal estimation
window_size <- 14 # number of days on each side of the moving window for marginal estimation 
maxlag <- 50 # How many lags to estimate tpdf 
r0_quant <- 0.975 # high radial quantile for tpdf estimation 
set.seed(92432) # noise sequence seed
seq_length <- 50000 # length of noise sequence to generate MA(q)'s for fitting
max_theta <- maxlag # maximum reasonable order of for MA(q)
chosen_order <- 10 # what order MA(q) do we want to use 




#load the data 
data.temp <- read_csv(file.name, show_col_types = FALSE)

#restrict to months we are considering
data.temp <- filter(data.temp, MON %in% 6:10)

#restrict to years we are considering 
data_early <- filter(data.temp, YR %in% past)
data_late <- filter(data.temp, YR %in% present)

#Get rid of "bad" years (bad if it has fewer than 50 days or fewer than 5 observations above some value)
keep <- c()

#past climate
#0.7 quantile of past BUI used to determine if a year has enough moderately high observations 
# to considered a full year   
quant_good <- quantile(data_early$BUI, probs = 0.7) 
for(yr in unique(data_early$YR)){
  data_temp <- filter(data_early, YR == yr)
  if((length(data_temp$BUI) > 50) && (sum(data_temp$BUI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_early <- filter(data_early, YR %in% keep)

#current climate
keep <- c()

for(yr in unique(data_late$YR)){
  data_temp <- filter(data_late, YR == yr)
  if((length(data_temp$BUI) > 50) && (sum(data_temp$BUI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_late <- filter(data_late, YR %in% keep)

# Add a variable "days_since" to represent days since June 1
data_early$days_since <- as.numeric(data_early$DATE - as.Date(paste(data_early$YR, "-05-31", sep = "")))
data_late$days_since <- as.numeric(data_late$DATE - as.Date(paste(data_late$YR, "-05-31", sep = "")))

# save quantiles (i.e. sort the (2k+1) day moving window of observations for each day)
ecdf_early <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)
ecdf_late <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)

for(i in 1:153){
  # consider only days between day i-k and i+k and sort the values on those (2k+1) days
  temp_data_early <- sort(data_early$BUI[which(data_early$days_since %in% 
                                                 max(1, i-window_size):min(153, i+window_size))])
  temp_data_late <- sort(data_late$BUI[which(data_late$days_since %in% 
                                               max(1, i-window_size):min(153, i+window_size))])
  
  # save the sorted days within the moving window as the ecdf quantiles
  ecdf_early[1:length(temp_data_early), i] <- temp_data_early
  ecdf_late[1:length(temp_data_late), i] <- temp_data_late
}


# save the threshold values
mu_early <- mu_late <- rep(NA, 153)

for (i in 1:153){
  q_temp_early <- floor(gpd_quant*(1+length(na.omit(ecdf_early[,i]))))
  mu_early[i] <- ecdf_early[q_temp_early, i]
  
  q_temp_late <- floor(gpd_quant*(1+length(na.omit(ecdf_late[,i]))))
  mu_late[i] <- ecdf_late[q_temp_late, i]
}

# Compute and save the amount that observations exceed the threshold. 
exceedance_early <- exceedance_late <- c()
yrs_early <- yrs_late <- c()
days_early <- days_late <- c()

for (i in 1:153){
  indices_temp <- which(data_early$days_since == i) #take all of the ith days
  temp_data_early <- data_early$BUI[indices_temp] 
  temp_ex_early <- temp_data_early - mu_early[i] #compute exceedances
  exceedance_early <- c(exceedance_early, temp_ex_early) #append
  yrs_early <- c(yrs_early, data_early$YR[indices_temp])
  days_early <- c(days_early, rep(i, length(indices_temp)))
  
  indices_temp <- which(data_late$days_since == i)
  temp_data_late <- data_late$BUI[indices_temp]
  temp_ex_late <- temp_data_late - mu_late[i]
  exceedance_late <- c(exceedance_late, temp_ex_late)
  yrs_late <- c(yrs_late, data_late$YR[indices_temp])
  days_late <- c(days_late, rep(i, length(indices_temp)))
}

# Rearrange exceedances so that they are in the same order as the data
exceedance_early <- data.frame(exceedance_early, yrs_early, days_early) %>% arrange(yrs_early)
exceedance_late <- data.frame(exceedance_late, yrs_late, days_late) %>% arrange(yrs_late)

# fit the GPD to the exceedances and save the likelihood value
single_model_early <- fpot(exceedance_early$exceedance_early, threshold = 0)
single_model_late <- fpot(exceedance_late$exceedance_late, threshold = 0)

#save parameter values
scale_early <- unname(single_model_early$estimate[1])
shape_early <- unname(single_model_early$estimate[2])
scale_late <- unname(single_model_late$estimate[1])
shape_late <- unname(single_model_late$estimate[2])

### Transform to uniform
#column is year
ts_uniform_early <- ts_uniform_late <- matrix(NA, nrow = 153, ncol = 20) 

#early
for(i in 1:length(data_early$BUI)){
  j <- data_early$days_since[i]
  yr <- which(unique(data_early$YR) == data_early$YR[i])
  if(data_early$BUI[i] < mu_early[j]) {
    ts_uniform_early[j, yr] <- (which(ecdf_early[,j] == data_early$BUI[i]) /
      (length(na.omit(ecdf_early[,j]))+1))[1]
  } else {
    ts_uniform_early[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_early$BUI[i], 
                                                                loc = mu_early[j],
                                                                scale = single_model_early$estimate[1], 
                                                                shape = single_model_early$estimate[2])
  }
}

#late
for(i in 1:length(data_late$BUI)){
  j <- data_late$days_since[i]
  yr <- which(unique(data_late$YR) == data_late$YR[i])
  if(data_late$BUI[i] < mu_late[j]) {
    ts_uniform_late[j, yr] <- (which(ecdf_late[,j] == data_late$BUI[i]) /
      (length(na.omit(ecdf_late[,j]))+1))[1]
  } else {
    ts_uniform_late[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_late$BUI[i], 
                                                               loc = mu_late[j],
                                                               scale = single_model_late$estimate[1], 
                                                               shape = single_model_late$estimate[2])
  }
}


############
#Frechet with alpha=2
###########
RV_ts_early <- 1/sqrt(log(1/ts_uniform_early))
RV_ts_late <- 1/sqrt(log(1/ts_uniform_late))


# Subtract mean (Bias correction)
me <- mean(na.omit(c(RV_ts_early)))
ml <- mean(na.omit(c(RV_ts_late)))

RV_ts_early <- ifelse(RV_ts_early - me > 0,
                      RV_ts_early - me, 0)
RV_ts_late <- ifelse(RV_ts_late - ml > 0,
                     RV_ts_late - ml, 0)



# first we make a plot of the daily 0.95 quantile when considering all years
quantiles_95_per_day <- rep(NA, 153)
seasonality_quantile <- 0.975
num_days <- c(30, 31, 31, 30, 31)
num_days_in_prev_months <- 0
for (i in 6:10){
  num_days_this_month <- num_days[(i-5)]
  for (j in 1:num_days_this_month){
    quantiles_95_per_day[j + num_days_in_prev_months] <- quantile(filter(data.temp, MON == i,
                                                                         DAY == j)$BUI,
                                                                  seasonality_quantile, type = 1)
  }
  num_days_in_prev_months <- num_days_in_prev_months + num_days_this_month
}

k <- 15
q <- c(0.25, 0.5, 0.7, 0.85, 0.975)

data.temp$days_since <- rep(seq(0,152,1), length(unique(data.temp$YR)))
quant <- matrix(NA, nrow = (154-k), ncol=5)
means <- rep(NA, 154-k)

for (i in 0:(153-k-1)) {
  temp <- filter(data.temp, days_since %in% c(i:(i+k)))$BUI
  quant[(i+1),] <- quantile(temp, q, type = 1, na.rm = T)
  means[(i+1)] <- mean(temp)
}

plot(quantiles_95_per_day, type = "l",
     xlab = "Days Since June 1", ylab = "BUI")
lines(8:146, quant[,5], col = "blue")


quants_early <- quants_late <- matrix(NA, nrow = 153, ncol = 6)

for(i in 1:153){
  quants_early[i,1:5] <- quantile(na.omit(ecdf_early[,i]), probs = q)
  quants_early[i,6] <- mean(na.omit(ecdf_early[,i]))
  
  quants_late[i, 1:5] <- quantile(na.omit(ecdf_late[,i]), probs = q)
  quants_late[i, 6] <- mean(na.omit(ecdf_late[,i]))
}


# Compare early to late 0.95 quantile 
dates <- format(as.Date(as.Date('2020-06-01'):as.Date('2020-10-31'), 
                        origin = '2020-01-01'), format = '%m-%d')

plot(quants_late[,5], type = "l", 
     ylim = c(min(quants_early[,5]) -5, max(quants_late[,5]) +5), xlim = c(0,153),
     xlab = "", ylab = "BUI", xaxt="n")
axis(1, at = c(seq(1, 153, 14),153), labels = dates[c(seq(1, 153, 14),153)], las=2)
lines(quants_early[,5], lty = 2)
legend(x = 117, y = 72, c("Present", "", "", "", "Past"), bty = "n")


# gpd parameters
parameter <- c("scale", "shape")
past_paras <- c(scale_early, shape_early)
present_paras <- c(scale_late, shape_late)

gpd_paras <- data.frame("." = parameter, "past" = round(past_paras, 4), "present" = round(present_paras, 4))

kable(gpd_paras)



#########
#Estimate TPDF
#########

# compute lags
lags_early <- lags_late <- matrix(NA, nrow = (153*20), ncol = maxlag+1)

for(i in 1:20){
  indices <- c(((i-1)*153+1):(i*153))
  lags_early[indices, (maxlag+1)] <- RV_ts_early[,i] #save transformed data into a vector
  lags_late[indices, (maxlag+1)] <- RV_ts_late[,i]
  for(j in 1:maxlag){
    lags_early[indices, j] <- lag(RV_ts_early[,i], j) 
    #append each years lagged vectors onto previous years
    lags_late[indices, j] <- lag(RV_ts_late[,i], j)
  }
}

# compute radial components
rt_early <- rt_late <- matrix(NA, nrow <- (153*20), ncol = maxlag)

for(j in 1:maxlag){
  rt_early[,j] <- sqrt(lags_early[,j]^2 + lags_early[,(maxlag+1)]^2)
  
  rt_late[,j] <- sqrt(lags_late[,j]^2 + lags_late[,(maxlag+1)]^2)
}

# compute "high threshold of radial component" (0.95 quantile of radial components)
r0_early <- r0_late <- rep(NA, maxlag)

for (i in 1:maxlag){
  len <- length(na.omit(rt_early[,i]))
  ind <- floor(r0_quant*len)
  r0_early[i] <- sort(na.omit(rt_early[,i]))[ind]
  
  len <- length(na.omit(rt_late[,i]))
  ind <- floor(r0_quant*len)
  r0_late[i] <- sort(na.omit(rt_late[,i]))[ind]
}


# Compute tpdf
tpdf_early <- tpdf_late <- rep(NA, maxlag)

for (i in 1:maxlag) {
  wi_temp <- lags_early[,i]*lags_early[,(maxlag+1)]/(rt_early[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_early[,i]), 0, rt_early[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_early[i]) # how many above radial threshold
  tpdf_early[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_early[i])])
  
  wi_temp <- lags_late[,i]*lags_late[,(maxlag+1)]/(rt_late[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_late[,i]), 0, rt_late[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_late[i]) # how many above radial threshold
  tpdf_late[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_late[i])])
}

#add lag 0 to tpdf
tpdf_early <- c(1, tpdf_early)
tpdf_late <- c(1, tpdf_late)  


#plot
plot(0:maxlag, tpdf_early, type = "h", ylim = c(0,1), main = "BUI tpdf CO")
legend("topright", legend = c("Past", "Present"), bty = "n", text.col = c("black", "blue"))
lines(0:maxlag + 0.1, tpdf_late, type = "h", ylim = c(0,1), col = "blue")






















#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rm(list = ls())
###############
########################################3
# COLORADO RAWS
#######################################
##############
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
















file.name <- "RAWS/Harbison_FWI_all.csv"

# The following are the values that we must choose for forward and back transformation
past <- 1959:1978
present <- 2002:2021
gpd_quant <- 0.975 # quantile for gpd in original marginal estimation
window_size <- 14 # number of days on each side of the moving window for marginal estimation 
maxlag <- 20 # How many lags to estimate tpdf 
r0_quant <- 0.975 # high radial quantile for tpdf estimation 
set.seed(92432) # noise sequence seed
seq_length <- 50000 # length of noise sequence to generate MA(q)'s for fitting
max_theta <- maxlag # maximum reasonable order of for MA(q)
chosen_order <- 10 # what order MA(q) do we want to use 


#load the data 
data.temp <- read_csv(file.name, show_col_types = FALSE)

#restrict to months we are considering
data.temp <- filter(data.temp, MON %in% 6:10)

#restrict to years we are considering 
data_early <- filter(data.temp, YR %in% past)
data_late <- filter(data.temp, YR %in% present)

#Get rid of "bad" years (bad if it has fewer than 50 days or fewer than 5 observations above some value)
keep <- c()

#past climate
#0.7 quantile of past ISI used to determine if a year has enough moderately high observations 
# to considered a full year   
quant_good <- quantile(data_early$ISI, probs = 0.7) 
for(yr in unique(data_early$YR)){
  data_temp <- filter(data_early, YR == yr)
  if((length(data_temp$ISI) > 50) && (sum(data_temp$ISI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_early <- filter(data_early, YR %in% keep)

#current climate
keep <- c()

for(yr in unique(data_late$YR)){
  data_temp <- filter(data_late, YR == yr)
  if((length(data_temp$ISI) > 50) && (sum(data_temp$ISI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_late <- filter(data_late, YR %in% keep)

# Add a variable "days_since" to represent days since June 1
data_early$days_since <- as.numeric(data_early$DATE - as.Date(paste(data_early$YR, "-05-31", sep = "")))
data_late$days_since <- as.numeric(data_late$DATE - as.Date(paste(data_late$YR, "-05-31", sep = "")))

# save quantiles (i.e. sort the (2k+1) day moving window of observations for each day)
ecdf_early <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)
ecdf_late <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)

for(i in 1:153){
  # consider only days between day i-k and i+k and sort the values on those (2k+1) days
  temp_data_early <- sort(data_early$ISI[which(data_early$days_since %in% 
                                                 max(1, i-window_size):min(153, i+window_size))])
  temp_data_late <- sort(data_late$ISI[which(data_late$days_since %in% 
                                               max(1, i-window_size):min(153, i+window_size))])
  
  # save the sorted days within the moving window as the ecdf quantiles
  ecdf_early[1:length(temp_data_early), i] <- temp_data_early
  ecdf_late[1:length(temp_data_late), i] <- temp_data_late
}


# save the threshold values
mu_early <- mu_late <- rep(NA, 153)

for (i in 1:153){
  q_temp_early <- floor(gpd_quant*(1+length(na.omit(ecdf_early[,i]))))
  mu_early[i] <- ecdf_early[q_temp_early, i]
  
  q_temp_late <- floor(gpd_quant*(1+length(na.omit(ecdf_late[,i]))))
  mu_late[i] <- ecdf_late[q_temp_late, i]
}

# Compute and save the amount that observations exceed the threshold. 
exceedance_early <- exceedance_late <- c()
yrs_early <- yrs_late <- c()
days_early <- days_late <- c()

for (i in 1:153){
  indices_temp <- which(data_early$days_since == i) #take all of the ith days
  temp_data_early <- data_early$ISI[indices_temp] 
  temp_ex_early <- temp_data_early - mu_early[i] #compute exceedances
  exceedance_early <- c(exceedance_early, temp_ex_early) #append
  yrs_early <- c(yrs_early, data_early$YR[indices_temp])
  days_early <- c(days_early, rep(i, length(indices_temp)))
  
  indices_temp <- which(data_late$days_since == i)
  temp_data_late <- data_late$ISI[indices_temp]
  temp_ex_late <- temp_data_late - mu_late[i]
  exceedance_late <- c(exceedance_late, temp_ex_late)
  yrs_late <- c(yrs_late, data_late$YR[indices_temp])
  days_late <- c(days_late, rep(i, length(indices_temp)))
}

# Rearrange exceedances so that they are in the same order as the data
exceedance_early <- data.frame(exceedance_early, yrs_early, days_early) %>% arrange(yrs_early)
exceedance_late <- data.frame(exceedance_late, yrs_late, days_late) %>% arrange(yrs_late)

# fit the GPD to the exceedances and save the likelihood value
single_model_early <- fpot(exceedance_early$exceedance_early, threshold = 0)
single_model_late <- fpot(exceedance_late$exceedance_late, threshold = 0)

#save parameter values
scale_early <- unname(single_model_early$estimate[1])
shape_early <- unname(single_model_early$estimate[2])
scale_late <- unname(single_model_late$estimate[1])
shape_late <- unname(single_model_late$estimate[2])

### Transform to uniform
#column is year
ts_uniform_early <- ts_uniform_late <- matrix(NA, nrow = 153, ncol = 20) 

#early
for(i in 1:length(data_early$ISI)){
  j <- data_early$days_since[i]
  yr <- which(unique(data_early$YR) == data_early$YR[i])
  if(data_early$ISI[i] < mu_early[j]) {
    ts_uniform_early[j, yr] <- min(which(ecdf_early[,j] == data_early$ISI[i]) /
                                     (length(na.omit(ecdf_early[,j]))+1))
  } else {
    ts_uniform_early[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_early$ISI[i], 
                                                                loc = mu_early[j],
                                                                scale = single_model_early$estimate[1], 
                                                                shape = single_model_early$estimate[2])
  }
}

#late
for(i in 1:length(data_late$ISI)){
  j <- data_late$days_since[i]
  yr <- which(unique(data_late$YR) == data_late$YR[i])
  if(data_late$ISI[i] < mu_late[j]) {
    ts_uniform_late[j, yr] <- which(ecdf_late[,j] == data_late$ISI[i]) /
      (length(na.omit(ecdf_late[,j]))+1)
  } else {
    ts_uniform_late[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_late$ISI[i], 
                                                               loc = mu_late[j],
                                                               scale = single_model_late$estimate[1], 
                                                               shape = single_model_late$estimate[2])
  }
}


############
#Frechet with alpha=2
###########
RV_ts_early <- 1/sqrt(log(1/ts_uniform_early))
RV_ts_late <- 1/sqrt(log(1/ts_uniform_late))


# Subtract mean (Bias correction)
me <- mean(na.omit(c(RV_ts_early)))
ml <- mean(na.omit(c(RV_ts_late)))

RV_ts_early <- ifelse(RV_ts_early - me > 0,
                      RV_ts_early - me, 0)
RV_ts_late <- ifelse(RV_ts_late - ml > 0,
                     RV_ts_late - ml, 0)


# Subtract mean (Bias correction)
me <- mean(na.omit(c(RV_ts_early)))
ml <- mean(na.omit(c(RV_ts_late)))

RV_ts_early <- ifelse(RV_ts_early - me > 0,
                      RV_ts_early - me, 0)
RV_ts_late <- ifelse(RV_ts_late - ml > 0,
                     RV_ts_late - ml, 0)



# first we make a plot of the daily 0.95 quantile when considering all years
quantiles_95_per_day <- rep(NA, 153)
seasonality_quantile <- 0.975
num_days <- c(30, 31, 31, 30, 31)
num_days_in_prev_months <- 0
for (i in 6:10){
  num_days_this_month <- num_days[(i-5)]
  for (j in 1:num_days_this_month){
    quantiles_95_per_day[j + num_days_in_prev_months] <- quantile(filter(data.temp, MON == i,
                                                                         DAY == j)$ISI,
                                                                  seasonality_quantile, type = 1)
  }
  num_days_in_prev_months <- num_days_in_prev_months + num_days_this_month
}

k <- 15
q <- c(0.25, 0.5, 0.7, 0.85, 0.975)

data.temp$days_since <- as.numeric(data.temp$DATE - 
                                     as.Date(paste(data.temp$YR, "-05-31", sep = "")))
quant <- matrix(NA, nrow = (154-k), ncol=5)
means <- rep(NA, 154-k)

for (i in 0:(153-k-1)) {
  temp <- filter(data.temp, days_since %in% c(i:(i+k)))$ISI
  quant[(i+1),] <- quantile(temp, q, type = 1, na.rm = T)
  means[(i+1)] <- mean(temp)
}

plot(quantiles_95_per_day, type = "l",
     xlab = "Days Since June 1", ylab = "ISI")
lines(8:146, quant[,5], col = "blue")


quants_early <- quants_late <- matrix(NA, nrow = 153, ncol = 6)

for(i in 1:153){
  quants_early[i,1:5] <- quantile(na.omit(ecdf_early[,i]), probs = q)
  quants_early[i,6] <- mean(na.omit(ecdf_early[,i]))
  
  quants_late[i, 1:5] <- quantile(na.omit(ecdf_late[,i]), probs = q)
  quants_late[i, 6] <- mean(na.omit(ecdf_late[,i]))
}



# Compare early to late 0.95 quantile 
dates <- format(as.Date(as.Date('2020-06-01'):as.Date('2020-10-31'), 
                        origin = '2020-01-01'), format = '%m-%d')

plot(quants_late[,5], type = "l", 
     ylim = c(min(quants_early[,5]) -5, max(quants_late[,5]) +5), xlim = c(0,153),
     xlab = "", ylab = "ISI", xaxt="n")
axis(1, at = c(seq(1, 153, 14),153), labels = dates[c(seq(1, 153, 14),153)], las=2)
lines(quants_early[,5], lty = 2)
legend(x = 117, y = 72, c("Present", "", "", "", "Past"), bty = "n")


# gpd parameters
parameter <- c("scale", "shape")
past_paras <- c(scale_early, shape_early)
present_paras <- c(scale_late, shape_late)

gpd_paras <- data.frame("." = parameter, "past" = round(past_paras, 4), "present" = round(present_paras, 4))

kable(gpd_paras)



#########
#Estimate TPDF
#########

# compute lags
lags_early <- lags_late <- matrix(NA, nrow = (153*20), ncol = maxlag+1)

for(i in 1:20){
  indices <- c(((i-1)*153+1):(i*153))
  lags_early[indices, (maxlag+1)] <- RV_ts_early[,i] #save transformed data into a vector
  lags_late[indices, (maxlag+1)] <- RV_ts_late[,i]
  for(j in 1:maxlag){
    lags_early[indices, j] <- lag(RV_ts_early[,i], j) 
    #append each years lagged vectors onto previous years
    lags_late[indices, j] <- lag(RV_ts_late[,i], j)
  }
}

# compute radial components
rt_early <- rt_late <- matrix(NA, nrow <- (153*20), ncol = maxlag)

for(j in 1:maxlag){
  rt_early[,j] <- sqrt(lags_early[,j]^2 + lags_early[,(maxlag+1)]^2)
  
  rt_late[,j] <- sqrt(lags_late[,j]^2 + lags_late[,(maxlag+1)]^2)
}

# compute "high threshold of radial component" (0.95 quantile of radial components)
r0_early <- r0_late <- rep(NA, maxlag)

for (i in 1:maxlag){
  len <- length(na.omit(rt_early[,i]))
  ind <- floor(r0_quant*len)
  r0_early[i] <- sort(na.omit(rt_early[,i]))[ind]
  
  len <- length(na.omit(rt_late[,i]))
  ind <- floor(r0_quant*len)
  r0_late[i] <- sort(na.omit(rt_late[,i]))[ind]
}


# Compute tpdf
tpdf_early <- tpdf_late <- rep(NA, maxlag)

for (i in 1:maxlag) {
  wi_temp <- lags_early[,i]*lags_early[,(maxlag+1)]/(rt_early[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_early[,i]), 0, rt_early[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_early[i]) # how many above radial threshold
  tpdf_early[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_early[i])])
  
  wi_temp <- lags_late[,i]*lags_late[,(maxlag+1)]/(rt_late[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_late[,i]), 0, rt_late[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_late[i]) # how many above radial threshold
  tpdf_late[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_late[i])])
}

#add lag 0 to tpdf
tpdf_early <- c(1, tpdf_early)
tpdf_late <- c(1, tpdf_late)  


#plot
plot(0:maxlag, tpdf_early, type = "h", ylim = c(0,1), main = "ISI tpdf CO")
legend("topright", legend = c("Past", "Present"), bty = "n", text.col = c("black", "blue"))
lines(0:maxlag + 0.1, tpdf_late, type = "h", ylim = c(0,1), col = "blue")










###########33
###############3
###########3
rm(list = ls())



file.name <- "RAWS/Harbison_FWI_all.csv"

# The following are the values that we must choose for forward and back transformation
past <- 1959:1978
present <- 2002:2021
gpd_quant <- 0.975 # quantile for gpd in original marginal estimation
window_size <- 14 # number of days on each side of the moving window for marginal estimation 
maxlag <- 50 # How many lags to estimate tpdf 
r0_quant <- 0.975 # high radial quantile for tpdf estimation 
set.seed(92432) # noise sequence seed
seq_length <- 50000 # length of noise sequence to generate MA(q)'s for fitting
max_theta <- maxlag # maximum reasonable order of for MA(q)
chosen_order <- 10 # what order MA(q) do we want to use 




#load the data 
data.temp <- read_csv(file.name, show_col_types = FALSE)

#restrict to months we are considering
data.temp <- filter(data.temp, MON %in% 6:10)

#restrict to years we are considering 
data_early <- filter(data.temp, YR %in% past)
data_late <- filter(data.temp, YR %in% present)

#Get rid of "bad" years (bad if it has fewer than 50 days or fewer than 5 observations above some value)
keep <- c()

#past climate
#0.7 quantile of past BUI used to determine if a year has enough moderately high observations 
# to considered a full year   
quant_good <- quantile(data_early$BUI, probs = 0.7) 
for(yr in unique(data_early$YR)){
  data_temp <- filter(data_early, YR == yr)
  if((length(data_temp$BUI) > 50) && (sum(data_temp$BUI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_early <- filter(data_early, YR %in% keep)

#current climate
keep <- c()

for(yr in unique(data_late$YR)){
  data_temp <- filter(data_late, YR == yr)
  if((length(data_temp$BUI) > 50) && (sum(data_temp$BUI > quant_good) >= 5)){
    keep <- c(keep, yr)
  } 
}

data_late <- filter(data_late, YR %in% keep)

# Add a variable "days_since" to represent days since June 1
data_early$days_since <- as.numeric(data_early$DATE - as.Date(paste(data_early$YR, "-05-31", sep = "")))
data_late$days_since <- as.numeric(data_late$DATE - as.Date(paste(data_late$YR, "-05-31", sep = "")))

# save quantiles (i.e. sort the (2k+1) day moving window of observations for each day)
ecdf_early <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)
ecdf_late <- matrix(NA, nrow = (2*window_size+1)*20, ncol = 153)

for(i in 1:153){
  # consider only days between day i-k and i+k and sort the values on those (2k+1) days
  temp_data_early <- sort(data_early$BUI[which(data_early$days_since %in% 
                                                 max(1, i-window_size):min(153, i+window_size))])
  temp_data_late <- sort(data_late$BUI[which(data_late$days_since %in% 
                                               max(1, i-window_size):min(153, i+window_size))])
  
  # save the sorted days within the moving window as the ecdf quantiles
  ecdf_early[1:length(temp_data_early), i] <- temp_data_early
  ecdf_late[1:length(temp_data_late), i] <- temp_data_late
}


# save the threshold values
mu_early <- mu_late <- rep(NA, 153)

for (i in 1:153){
  q_temp_early <- floor(gpd_quant*(1+length(na.omit(ecdf_early[,i]))))
  mu_early[i] <- ecdf_early[q_temp_early, i]
  
  q_temp_late <- floor(gpd_quant*(1+length(na.omit(ecdf_late[,i]))))
  mu_late[i] <- ecdf_late[q_temp_late, i]
}

# Compute and save the amount that observations exceed the threshold. 
exceedance_early <- exceedance_late <- c()
yrs_early <- yrs_late <- c()
days_early <- days_late <- c()

for (i in 1:153){
  indices_temp <- which(data_early$days_since == i) #take all of the ith days
  temp_data_early <- data_early$BUI[indices_temp] 
  temp_ex_early <- temp_data_early - mu_early[i] #compute exceedances
  exceedance_early <- c(exceedance_early, temp_ex_early) #append
  yrs_early <- c(yrs_early, data_early$YR[indices_temp])
  days_early <- c(days_early, rep(i, length(indices_temp)))
  
  indices_temp <- which(data_late$days_since == i)
  temp_data_late <- data_late$BUI[indices_temp]
  temp_ex_late <- temp_data_late - mu_late[i]
  exceedance_late <- c(exceedance_late, temp_ex_late)
  yrs_late <- c(yrs_late, data_late$YR[indices_temp])
  days_late <- c(days_late, rep(i, length(indices_temp)))
}

# Rearrange exceedances so that they are in the same order as the data
exceedance_early <- data.frame(exceedance_early, yrs_early, days_early) %>% arrange(yrs_early)
exceedance_late <- data.frame(exceedance_late, yrs_late, days_late) %>% arrange(yrs_late)

# fit the GPD to the exceedances and save the likelihood value
single_model_early <- fpot(exceedance_early$exceedance_early, threshold = 0)
single_model_late <- fpot(exceedance_late$exceedance_late, threshold = 0)

#save parameter values
scale_early <- unname(single_model_early$estimate[1])
shape_early <- unname(single_model_early$estimate[2])
scale_late <- unname(single_model_late$estimate[1])
shape_late <- unname(single_model_late$estimate[2])

### Transform to uniform
#column is year
ts_uniform_early <- ts_uniform_late <- matrix(NA, nrow = 153, ncol = 20) 

#early
for(i in 1:length(data_early$BUI)){
  j <- data_early$days_since[i]
  yr <- which(unique(data_early$YR) == data_early$YR[i])
  if(data_early$BUI[i] < mu_early[j]) {
    ts_uniform_early[j, yr] <- (which(ecdf_early[,j] == data_early$BUI[i]) /
                                  (length(na.omit(ecdf_early[,j]))+1))[1]
  } else {
    ts_uniform_early[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_early$BUI[i], 
                                                                loc = mu_early[j],
                                                                scale = single_model_early$estimate[1], 
                                                                shape = single_model_early$estimate[2])
  }
}

#late
for(i in 1:length(data_late$BUI)){
  j <- data_late$days_since[i]
  yr <- which(unique(data_late$YR) == data_late$YR[i])
  if(data_late$BUI[i] < mu_late[j]) {
    ts_uniform_late[j, yr] <- (which(ecdf_late[,j] == data_late$BUI[i]) /
                                 (length(na.omit(ecdf_late[,j]))+1))[1]
  } else {
    ts_uniform_late[j, yr] <- gpd_quant + (1-gpd_quant) * pgpd(data_late$BUI[i], 
                                                               loc = mu_late[j],
                                                               scale = single_model_late$estimate[1], 
                                                               shape = single_model_late$estimate[2])
  }
}


############
#Frechet with alpha=2
###########
RV_ts_early <- 1/sqrt(log(1/ts_uniform_early))
RV_ts_late <- 1/sqrt(log(1/ts_uniform_late))


# Subtract mean (Bias correction)
me <- mean(na.omit(c(RV_ts_early)))
ml <- mean(na.omit(c(RV_ts_late)))

RV_ts_early <- ifelse(RV_ts_early - me > 0,
                      RV_ts_early - me, 0)
RV_ts_late <- ifelse(RV_ts_late - ml > 0,
                     RV_ts_late - ml, 0)



# first we make a plot of the daily 0.95 quantile when considering all years
quantiles_95_per_day <- rep(NA, 153)
seasonality_quantile <- 0.975
num_days <- c(30, 31, 31, 30, 31)
num_days_in_prev_months <- 0
for (i in 6:10){
  num_days_this_month <- num_days[(i-5)]
  for (j in 1:num_days_this_month){
    quantiles_95_per_day[j + num_days_in_prev_months] <- quantile(filter(data.temp, MON == i,
                                                                         DAY == j)$ISI,
                                                                  seasonality_quantile, type = 1)
  }
  num_days_in_prev_months <- num_days_in_prev_months + num_days_this_month
}

k <- 15
q <- c(0.25, 0.5, 0.7, 0.85, 0.975)

data.temp$days_since <- as.numeric(data.temp$DATE - 
                                     as.Date(paste(data.temp$YR, "-05-31", sep = "")))
quant <- matrix(NA, nrow = (154-k), ncol=5)
means <- rep(NA, 154-k)

for (i in 0:(153-k-1)) {
  temp <- filter(data.temp, days_since %in% c(i:(i+k)))$ISI
  quant[(i+1),] <- quantile(temp, q, type = 1, na.rm = T)
  means[(i+1)] <- mean(temp)
}

plot(quantiles_95_per_day, type = "l",
     xlab = "Days Since June 1", ylab = "ISI")
lines(8:146, quant[,5], col = "blue")


quants_early <- quants_late <- matrix(NA, nrow = 153, ncol = 6)

for(i in 1:153){
  quants_early[i,1:5] <- quantile(na.omit(ecdf_early[,i]), probs = q)
  quants_early[i,6] <- mean(na.omit(ecdf_early[,i]))
  
  quants_late[i, 1:5] <- quantile(na.omit(ecdf_late[,i]), probs = q)
  quants_late[i, 6] <- mean(na.omit(ecdf_late[,i]))
}


# Compare early to late 0.95 quantile 
dates <- format(as.Date(as.Date('2020-06-01'):as.Date('2020-10-31'), 
                        origin = '2020-01-01'), format = '%m-%d')

plot(quants_late[,5], type = "l", 
     ylim = c(min(quants_early[,5]) -5, max(quants_late[,5]) +5), xlim = c(0,153),
     xlab = "", ylab = "BUI", xaxt="n")
axis(1, at = c(seq(1, 153, 14),153), labels = dates[c(seq(1, 153, 14),153)], las=2)
lines(quants_early[,5], lty = 2)
legend(x = 117, y = 72, c("Present", "", "", "", "Past"), bty = "n")


# gpd parameters
parameter <- c("scale", "shape")
past_paras <- c(scale_early, shape_early)
present_paras <- c(scale_late, shape_late)

gpd_paras <- data.frame("." = parameter, "past" = round(past_paras, 4), "present" = round(present_paras, 4))

kable(gpd_paras)



#########
#Estimate TPDF
#########

# compute lags
lags_early <- lags_late <- matrix(NA, nrow = (153*20), ncol = maxlag+1)

for(i in 1:20){
  indices <- c(((i-1)*153+1):(i*153))
  lags_early[indices, (maxlag+1)] <- RV_ts_early[,i] #save transformed data into a vector
  lags_late[indices, (maxlag+1)] <- RV_ts_late[,i]
  for(j in 1:maxlag){
    lags_early[indices, j] <- lag(RV_ts_early[,i], j) 
    #append each years lagged vectors onto previous years
    lags_late[indices, j] <- lag(RV_ts_late[,i], j)
  }
}

# compute radial components
rt_early <- rt_late <- matrix(NA, nrow <- (153*20), ncol = maxlag)

for(j in 1:maxlag){
  rt_early[,j] <- sqrt(lags_early[,j]^2 + lags_early[,(maxlag+1)]^2)
  
  rt_late[,j] <- sqrt(lags_late[,j]^2 + lags_late[,(maxlag+1)]^2)
}

# compute "high threshold of radial component" (0.95 quantile of radial components)
r0_early <- r0_late <- rep(NA, maxlag)

for (i in 1:maxlag){
  len <- length(na.omit(rt_early[,i]))
  ind <- floor(r0_quant*len)
  r0_early[i] <- sort(na.omit(rt_early[,i]))[ind]
  
  len <- length(na.omit(rt_late[,i]))
  ind <- floor(r0_quant*len)
  r0_late[i] <- sort(na.omit(rt_late[,i]))[ind]
}


# Compute tpdf
tpdf_early <- tpdf_late <- rep(NA, maxlag)

for (i in 1:maxlag) {
  wi_temp <- lags_early[,i]*lags_early[,(maxlag+1)]/(rt_early[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_early[,i]), 0, rt_early[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_early[i]) # how many above radial threshold
  tpdf_early[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_early[i])])
  
  wi_temp <- lags_late[,i]*lags_late[,(maxlag+1)]/(rt_late[,i]^2) #angular components
  rt_temp <- ifelse(is.na(rt_late[,i]), 0, rt_late[,i]) #set NA's to zer
  k <- sum(rt_temp > r0_late[i]) # how many above radial threshold
  tpdf_late[i] <- 2/k * sum(wi_temp[which(rt_temp>r0_late[i])])
}

#add lag 0 to tpdf
tpdf_early <- c(1, tpdf_early)
tpdf_late <- c(1, tpdf_late)  


#plot
plot(0:maxlag, tpdf_early, type = "h", ylim = c(0,1), main = "BUI tpdf CO")
legend("topright", legend = c("Past", "Present"), bty = "n", text.col = c("black", "blue"))
lines(0:maxlag + 0.1, tpdf_late, type = "h", ylim = c(0,1), col = "blue")












# for plotting gpd's
#curve(dgpd(x, loc = 0, scale =3.3964, shape = -0.2708), from = 0, to = 15, ylim = c(0, 0.35))
#curve(dgpd(x, loc = 0, scale = 2.8267, shape = -0.0799), from = 0, to = 15, add=T, col = "blue")
