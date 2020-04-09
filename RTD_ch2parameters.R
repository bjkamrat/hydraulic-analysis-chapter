# Author: Brock Kamrath
# Date Created: 4/8/2020
# Built using the priniciples from Bodin et al (2013) and Black and Martinez (2003)

#Objective: Hydraulic analysis of common dataset format for tracer tests at walnut Cove

#load packages
rm(list=ls(all=TRUE))
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(readxl)
library(formattable)


data <- read_excel("hydraulic_analysis/tt_o24.xlsx", sheet = "main_ext")
head(data)
  
# Wetland Description
# plate # at outlet weir
weir_num <- 4

# Volume of tracer used (L)
trace_vol <- 0.15

if(weir_num == 4){
  #nominal wetland volume (L)
  V <- 1750033.722
} else {
  #nominal wetland volume (L)
  V <- 1031524.711
}


#Mass of tracer added (ug)
mass_added <- 238000000*trace_vol
m10 <- 0.1*mass_added

#find the time weighted average flow during the tracer test
# Time weighted average flow (Lpd)
Q_avg <- sum(data$Q_Lpd*data$deltat_d)/sum(data$deltat_d)
# Time weighted average flow (gpm)
Q_avg_gpm <- Q_avg/(3.785*24*60)

#Nominal retention time (days)
t_nom <- V/Q_avg


#Calculate the theta value
data <- data%>%
  mutate(time_prev = lag(time_d),
         Q_Lpd_prev = lag(Q_Lpd),
         C_ugL_prev = lag(C_ugL))

#set a vector for f(t) to NA values
data$f_t <- as.numeric(rep(NA,nrow(data)))
# set the first value to 0
data$f_t[1] <- 0

data <- data%>%
  mutate(M0 = (((Q_Lpd_prev*C_ugL_prev)+(Q_Lpd*C_ugL))/2)*(deltat_d),
         f_t = (Q_Lpd*C_ugL)/sum(M0,na.rm=TRUE),
         M1 = (((time_prev*Q_Lpd_prev*C_ugL_prev)+(time_d*Q_Lpd*C_ugL))/2)*(deltat_d),
         M2 = ((((time_prev^2)*Q_Lpd_prev*C_ugL_prev)+((time_d^2)*Q_Lpd*C_ugL))/2)*(deltat_d)
         )

#mass recovered (%)
mass_out <- sum(data$M0,na.rm=TRUE)
mass_rec <- (mass_out/mass_added)*100

# Check on values for unity, should be 1
unity <- sum(data$f_t*data$deltat_d)
# Actual retention time, First normalized moment (d)
tau <- sum(data$M1,na.rm=TRUE)/sum(data$M0,na.rm=TRUE)

# Effective volume as decimal
e <- tau/t_nom

# The normalized second moment of the RTD, varience (d^2)
sigma2 <- (sum(data$M2,na.rm=TRUE)/sum(data$M0,na.rm=TRUE))-(tau^2)

#N value for method of moments
N <- (tau^2)/sigma2

# λe is similar to λm due to the numerically equal value between λm and e, 
# but is animprovement on λm. λe reflects the RTD shape due to addition of N in
# the formula
lamda_e <- e*(1-(1/N))

#sigma theta is a dimensionless measure of dispersive processes
sigma_theta2 <- sigma2/(tau^2)


# Calucualte the when 10% of the mass added reaches the outlet relative to t_nom
M0 <- na.omit(data$M0)
mass_cum <- cumsum(M0)
time_d <- data$time_d[-1]

for(i in 1:length(mass_cum)){
  if(m10 > mass_cum[i]){
    t10 <- ((m10-mass_cum[i])/(mass_cum[i+1]-mass_cum[i])*(time_d[i+1]-time_d[i]))+time_d[i]
  }
}

#normalize t10
t10n <- t10/t_nom


obs <-  data.frame(data$time_d, data$f_t)
colnames(obs)

obs1 <- obs %>% 
  setNames(c('date','F_t'))


ggplot(data)+
  geom_point(aes(x=time_d,y = f_t))+
  geom_vline(aes(xintercept=tau))

#Moment analysis using Gamma distribution
gamma <- function(data,N,tau){
  dat <- data%>%
  mutate(f_t = dgamma(time_d, shape = N, scale=tau/N))
  
  outputs <- data.frame(date=dat$time_d,f_t=dat$f_t)
  return(outputs)
}

  
#1.par: Define initial values - N, tau----
#use nominal values provided
initial_values <- c(N,tau)

#2.fn: function that returns the SSE ----
#for model, which will be minimized
error <- function(parameters,observed_repsonses,explanatory_variables){
  #Assign initial values to parameters
  #N value from data
  N <- parameters[1]
  #tau
  tau <- parameters[2]
  
  #Run model
  simulated_responses <- gamma(explanatory_variables,N,tau)
  
  #Merge with observations
  sim_obs <- inner_join(observed_repsonses,simulated_responses,by="date")
  
  
  #Calculate SSE
  sse <- sum((sim_obs$F_t-sim_obs$f_t)^2)
  
  #Return
  return(sse)
}

#read in weather file
dat <- data
observations <- obs1

# Calibrate
cal <- optim(par=initial_values,
             fn=error,
             gr=NULL,
             observed_repsonses=observations,
             explanatory_variables=dat,
             method="Nelder-Mead")

cal

#apply calibrated maize model
simulation <- gamma(data,N=N,tau=tau)
simulation_calibrated <- gamma(data,N=cal$par[1],tau=cal$par[2])

e_gamma <- cal$par[2]/t_nom
lamda_e_gamma <- (cal$par[2]/t_nom)*(1-(1/cal$par[1]))
sigma_theta2_gamma <- 1/cal$par[1]

results <- c(mass_rec,Q_avg,e,N,tau,lamda_e,sigma_theta2,e_gamma,cal$par[1],cal$par[2],lamda_e_gamma,sigma_theta2_gamma)
results <- formattable(results, digits = 2, format = "f")
results
#########################################################################

ggplot()+
  geom_point(data=data,aes(time_d,f_t))+
  geom_line(data=simulation,aes(date,f_t))+
  geom_line(data=simulation_calibrated,aes(date,f_t),color="red")+
  geom_vline(aes(xintercept=t_nom),color="blue")+
  labs(x="Time, day",y="f(t)")+
  theme_classic()
