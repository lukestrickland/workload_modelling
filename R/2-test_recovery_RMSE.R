library(FIPS)
library(tidyverse)

#First try recovery of three process model with 8 hours of sleep once per day at the same times

load("img/data_list_8h.RData")

free_ps <- c(d=-0.1, g = -0.5, Ca= 2.5, Ua=0.5, S0=3)
fixed_ps <- pvec.threeprocess[!(names(pvec.threeprocess)
                                %in% names(free_ps))]
#bounds
# decay rates (g, d) must be -ve
# g can't be exactly 0
# S0 cannot be higher than upper asymptote or lower than lower asymptote


lower <- c(-Inf,-Inf,-Inf,-Inf, 2.4)
upper <- c(-0.001,-0.001,Inf,Inf, 14.3)

#bad_datetimes <- full_df$datetime[1:2]

#full_df <- full_df %>% filter(!datetime %in% bad_datetimes)

FIPS_errsquared <- function(free, fixed, dat, modeltype = "TPM", free_names = NULL) {
  if (!is.null(free_names)) 
    names(free) <- free_names
  ps <- c(free, fixed)
  simmed_data = FIPS_simulation(FIPS_df = dat, modeltype = modeltype, 
                                pvec = ps)$FIPS_df
  (simmed_data$KSS[simmed_data$event=="observation"] - 
    dat$KSS[simmed_data$event=="observation"]
    )^2
}

FIPS_GROUP_MSE <- function(free, fixed, dat, modeltype = "TPM", free_names = NULL){
  sqerr <- do.call("c", dat %>% group_by(subject) %>% 
                     group_map(~ FIPS_errsquared(
                       free, fixed, .x, free_names=names(free), modeltype=modeltype)))
  mean(sqerr)
}

FIPS_errsquared(free_ps, fixed_ps, full_df %>% filter(subject==7))

result_nlminb <- nlminb(free_ps, FIPS_GROUP_MSE, fixed=fixed_ps, dat=full_df,
                        free_names = names(free_ps), upper=upper,
                        lower=lower)

#Basically, it's bad, even with just a few free parameters



load("img/data_list_sleepdep.RData")

#Add some noise
full_df$lapses <- full_df$lapses + rnorm(length(full_df$lapses),
                       0, 1)

#Remove decimals (couldn't exist in observations)
full_df$lapses <- round(full_df$lapses)

#try recovery of the unified model with ramarishkan's data 
FIPS_errsquared <- function(free, fixed, dat, modeltype = "TPM", free_names = NULL) {
  if (!is.null(free_names)) 
    names(free) <- free_names
  ps <- c(free, fixed)
  simmed_data = FIPS_simulation(FIPS_df = dat, modeltype = modeltype, 
                                pvec = ps)$FIPS_df
  (simmed_data$lapses[simmed_data$event=="observation"] - 
      dat$lapses[simmed_data$event=="observation"]
  )^2
}

FIPS_GROUP_MSE <- function(free, fixed, dat, modeltype = "TPM", free_names = NULL){
  sqerr <- do.call("c", dat %>% group_by(subject) %>% 
                     group_map(~ FIPS_errsquared(
                       free, fixed, .x, free_names=names(free), modeltype=modeltype)))
  mean(sqerr)
}


#LUKE: not totally sure if these bounds are good
lower <- c(0,0,0,0, -0.11, -0.11, -Inf)
upper <- c(Inf,Inf, Inf,Inf, Inf, Inf, Inf)

#LUKE: Decent recovery from both the below sets of start points
#start points 1
free_ps<- c(kappa=2, tau_la=75, tau_s=0.5, tau_w=25, S0=7, U0=50, phi=5)
#start points 2
#free_ps<- c(kappa=8, tau_la=150, tau_s=2, tau_w=60, S0=0, U0=85, phi=0.5)

rbind(lower, free_ps, upper)

fixed_ps <- unified_make_pvec()[!(names(unified_make_pvec())
                                %in% names(free_ps))]

result_nlminb <- nlminb(free_ps, FIPS_GROUP_MSE, fixed=fixed_ps, dat=full_df,
                        free_names = names(free_ps), upper=upper,
                        lower=lower, modeltype="unified")

recovery_df <- cbind(free_ps, result_nlminb$par, unified_make_pvec()[names(free_ps)])

colnames(recovery_df) <- c("start_points", "recovered", "true")

recovery_df


