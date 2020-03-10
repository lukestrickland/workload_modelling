#clear workspace
rm(list=ls())

#load packages
library(rstan)
library(FIPS)
library(truncnorm)

#load data
load("img/data_list_sleepdep.RData")

#Randomly sample parameter values for each participant

pvec <- unified_make_pvec()
names(pvec)[names(pvec)=="S0"] <- "S0_raw"
names(pvec)[names(pvec)=="L0"] <- "L0_raw"
names(pvec)[names(pvec)=="tau_s"] <- "tau_d"
names(pvec)[names(pvec)=="tau_w"] <- "tau_r"

#Transformations to put the L0 and S0 in a better format for constraints
pvec['S0_raw'] <- (pvec['S0_raw'] - pvec['L0_raw'])/(pvec['U0'] - pvec['L0_raw'])
pvec['L0_raw'] <- pvec['L0_raw']/pvec['U0']

pvec["sigma"] <- 0.5


p_list <- as.list(pvec)

fit_list <- c(data_list, p_list)
fit_list$fatigue <- rep(0, fit_list$Ntotal)



fit_generate <- stan( file="Stan/3PM_unified_group_gen_generating.stan" , 
                      data = fit_list , 
                      pars = c("pp"),
                      include = TRUE,
                      cores=1,
                      iter=1,
                      chains = 1,
                      seed=1234 )

fit_list$fatigue = round(as.vector(rstan::extract(fit_generate,"pp")$pp))

fit_recovery_fixed_tau <- stan( file="Stan/3PM_unified_group_gen_now.stan" , 
                      data = fit_list , 
                      #pars = c("pp"),
                      #include = TRUE,
                      cores=4,
                      #iter=1,
                      chains = 4,
                      seed=1234 )

save(fit_recovery_fixed_tau, file = 
    "C:/Users/282952E/Dropbox/SAMPLES/workload_modelling/recovery_unified_fixed_tau.RData"
    )


#Analyse recovery
smry = summary(fit_recovery_fixed_tau)
smry$summary[1:10,]



generating = tibble(L0 = fit_list$L0_raw,
                    S0 = fit_list$S0_raw,
                    U0 = fit_list$U0,
                    tau_r = fit_list$tau_r,
                    tau_d = fit_list$tau_d,
                    sigma = fit_list$sigma,
                    kappa = fit_list$kappa,
                    phi = fit_list$phi,
                    wc = fit_list$wc,
                    wd = fit_list$wd
)

pairs(generating)

recovered = tibble( L0_raw = smry[[1]][grep("L0_raw", rownames(smry[[1]])),"mean"],
                    S0_raw = smry[[1]][grep("S0_raw", rownames(smry[[1]])),"mean"],
                    U0 = smry[[1]][grep("U0", rownames(smry[[1]])),"mean"],
                    kappa = smry[[1]][grep("kappa", rownames(smry[[1]])),"mean"],
                    phi = smry[[1]][grep("phi", rownames(smry[[1]])),"mean"]
)

recovered

#Craziness
pairs(fit_recovery_fixed_tau,pars = c("tau_d", "tau_r", "U0", "S0_raw", "sigma"))

#Try just the parameters from the paper - c(kappa=2, tau_la=75, tau_s=0.5, tau_w=25, S0=7, U0=50, phi=5)



