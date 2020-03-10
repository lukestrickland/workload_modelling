
functions{
  //Equations from Borbely & Acherman (1999). Sleep Homeostasis and Models of Sleep Regulation. Journal of Biological Rythems
  
  //Note also that Rajdev (2013) Theoretical Biology has
  //some good info on estimating the model. They estimate   //S0, U (upper assymptote), L (lower assymptote), tau_w, tau_s, phi, kappa. Phi represents the phase (t0 I think), and kappa represents the relative influence of C over S.
  
  //Higher S and C indicate greater performance impairment
  
  //calculates S during wake
  real Sfun(real sw, //sw = S upon waking
             real taw, //taw = time awake
             real tau_r, //tau_d = controls rate of rise in S during wake
             real U0
             ){ 
    
    real r=exp(-taw/tau_r);        
    real S=U0-r*(U0-sw);
    return S;
  }
  
  //calculates S during sleep
  real Spfun(real ss, //ss = S upon falling asleep
             real tas, //tas = time asleep
             real tau_s, //tau_d = controls rate of decay in S during sleep
             real U0,     //upper assymptote
             real tau_la, //rate of change in lower assymptote
             real ls     //lower assymptote at sleep onset
             ){ 
    
    real term1 = ss*exp(-tas/tau_s);
    real term2 = -2*U0*(1-exp(-tas/tau_s));
    real term3 = (((ls + 2*U0)*tau_la)/(tau_la-tau_s)) * (exp(-tas/tau_la)-exp(-tas/tau_s));            
    real Sp=term1+term2+term3;
    return Sp;
  }
  
 //calculates L during wake
  real Lfun(real lw, //lower assymptote upon falling asleep
           real taw,    //time asleep
           real tau_la, //rate of change in lower assymptote
           real U0      //upper assymptote
          ){
    real term1 = lw*exp(-taw/tau_la);
    real term2 = U0*(1-exp(-taw/tau_la));
    real L=term1+term2;
    return L;
  }
  
  
  
  //calculates L during sleep
  real Lpfun(real ls, //lower assymptote upon falling asleep
           real tas,    //time asleep
           real tau_la, //rate of change in lower assymptote
           real U0      //upper assymptote
          ){
    real term1 = ls*exp(-tas/tau_la);
    real term2 = -2*U0*(1-exp(-tas/tau_la));
    real L=term1+term2;
    return L;
  }
  
  
  //calculates C (circadian process)
  real Cfun(real tod, //tod = time of day (in decimal hours)
            real phi,  //phi = phase at beginning of the simulation (I think this should be 0 if t = tod)
            real tau, //tau = period of C process
            real A //amplitute of process
            ){
    
    real omega = 2*pi()/tau;
    real term1 = 0.97*sin(omega*(tod+phi));
    real term2 = 0.22*sin(2*omega*(tod+phi));
    real term3 = 0.07*sin(3*omega*(tod+phi));
    real term4 = 0.03*sin(4*omega*(tod+phi));
    real term5 = 0.0001*sin(5*omega*(tod+phi));
    real C = A*(term1+term2+term3+term4+term5);
    return C;
  }
  
  //sleep intertia function. From 3PM. Calculates
  //effect of sleep intertia on alterness
  real Wfun(real taw,
            real wc, //Wc = extent of alterness reduction at time of waking (typically = -5.72, but sign is reversed for 2PM)
            real wd  //Wd = exponential recovery of alterness (typically = -1.51)
            ){
              
    real W=wc*exp(wd*taw);
    return W;
  }
}

data {
  int<lower=0> Nsubj;
  int<lower=0> Ntotal;
  int<lower=0> subject[Ntotal];
  int<lower=0> event_number[Ntotal];
  int<lower=0> previous_episode_type[Ntotal];
  real<lower=0> time_since_previous[Ntotal];
  vector<lower=0>[Ntotal] timeofday;
  vector<lower=0>[Ntotal] fatigue;
  int<lower=0> Nvalid;
  int<lower=0> valid[Nvalid];
  real<lower=0> U0;
  real<lower= -0.11, upper=1> S0_raw;
  real<lower=-0.11,upper=1> L0_raw;
  real phi;
  real<lower=0> kappa;
  real<lower=0> tau_d; //typically fixed to 4.2  h
  real<lower=0> tau_r; //typically fixed to 18.2 h
  real<lower=0> tau_la;
 // real<lower=0> A;     //typically fixed to 0.12
  real<lower=0> sigma; //typically not modeled, but necessary for parameter estimation
  real<lower=0> wc;
  real<upper=0> wd;
}  
 
parameters {

  //parameters
  real<lower=0,upper=1> dummy;
  
}

transformed parameters {

  //Fixed parameters
  //real t0 = -0.6;
  real tau = 24;
  real A = 1;
  real timeawake;
  //real U0 = L0+U0_raw;
  //real S0 = L0+U0_raw*S0_raw;
  
  //Calculate level of processes for each observation
  real s_prev;  //level of homeostatic process at previous event
  real l_prev;  //level of lower assymptote
  vector[Ntotal] S; //level of homeostatic process
  vector[Ntotal] C; //level of 24-hour circadian process
  vector[Ntotal] L; //lower assumptote of homeostatic process
  vector[Ntotal] W;
  

    //Centre parameters
  real L0 = L0_raw*U0;
  real S0 = L0 + (U0-L0)*S0_raw; 

  //Centre parameters
  // real L0 = (L0_raw*3 - 2)*U0;
  // real S0 = L0 + (U0-L0)*S0_raw; 

  //Calculate S and C for each event
  for(i in 1:Ntotal){
    //if it is subject's first event, assign S to be S0;
    if(event_number[i] == 1){
      S[i] = S0;
      L[i] = L0;
      s_prev = S[i];
      l_prev = L[i];
      W[i] = 0;
      timeawake = 0;
    }
    //if it is not subject's first event, update S based on S at previous event
    if(event_number[i] > 1){
      //if most recent episode was sleep (1)
      if(previous_episode_type[i] == 1){
        S[i] = Spfun(s_prev,time_since_previous[i],tau_d,U0,tau_la,l_prev);
        L[i] = Lpfun(l_prev,time_since_previous[i],tau_la,U0);
        s_prev = S[i];
        l_prev = L[i]; //L process only updates when sleeps
        W[i] = 0;
        timeawake = 0;
      }
      //if most recent episode was wake (2)
      if(previous_episode_type[i] == 2){
        S[i] = Sfun(s_prev,time_since_previous[i],tau_r,U0);
         L[i] = Lfun(l_prev,time_since_previous[i],tau_la,U0);
        s_prev = S[i];
        l_prev = L[i];
        timeawake += time_since_previous[i];
        W[i] = Wfun(timeawake,wc,wd);
      }
      //if most recent episode was work (3)
     // if(previous_episode_type[i] == 3){
    //    S[i] = Sfun(s_prev,time_since_previous[i],tau_r_w,U0);
    //    s_prev = S[i];
    //  }
    }  
    
    C[i]=Cfun(timeofday[i],phi,tau,A);
    
  }
}

model {
  
   //hyperparameters
  // S0_a ~ normal(0,4);
  // S0_b ~ normal(0,4);
  // 
  // U0_mean ~ normal(0,5);
  // U0_sd ~ normal(0,2);
  // 
  // L0_a ~ normal(0,4);
  // L0_b ~ normal(0,4);
  // 
  // phi_mean ~ normal(0,5);
  // phi_sd ~ normal(0,2);
  // 
  // kappa_mean ~ normal(0,5);
  // kappa_sd ~ normal(0,2);
  // 
  // tau_d_mean ~ normal(0,10); 
  // tau_d_sd ~ normal(0,2); 
  // 
  // tau_r_mean ~ normal(0,10); 
  // tau_r_sd ~ normal(0,2);
  // 
  // tau_la_mean ~ normal(0,10); 
  // tau_la_sd ~ normal(0,2);
  // 
  // sigma_mean ~ normal(0,3);
  // sigma_sd ~ normal(0,2);
  // 
  // wc_mean ~ normal(0,5);
  // wc_sd ~ normal(0,2);
  // 
  // wd_mean ~ normal(0,5);
  // wd_sd ~ normal(0,2);
  // 
  // 
  // //parameters
  // S0_raw ~ beta(S0_a,S0_b);
  // U0_raw ~ normal(0,1);
  // L0_raw ~ beta(L0_a,L0_b);
  // phi_raw ~ normal(0,1);
  // kappa_raw ~ normal(0,1);
  // tau_d_raw ~ normal(0,1); //typically fixed to 4.2  h
  // tau_r_raw ~ normal(0,1); //typically fixed to 18.2 h
  // tau_la_raw ~ normal(0,1); 
  // sigma_raw ~ normal(0,1);
  // wc_raw ~ normal(0,1);
  // wd_raw ~ normal(0,1);
  // 
  // //S0 ~ normal(0,10);
  // // U0_raw ~ normal(0,10);
  // // L0 ~ normal(0,10);
  // // phi ~ normal(0,5);
  // // kappa ~ normal(0,5);
  // // tau_d ~ normal(0,10); //typically fixed to 4.2  h
  // // tau_r ~ normal(0,10); //typically fixed to 18.2 h
  // // //tau_r_nw ~ normal(0,10); //typically fixed to 18.2 h
  // // //A ~ normal(0,1);     //typically fixed to 0.12
  // // sigma ~ normal(0,3);
  // 
  //  for(i in 1:Nvalid){
  //   fatigue[valid[i]] ~ normal(S[valid[i]]+W[valid[i]]+kappa*C[valid[i]],sigma);
  // } 
}

generated quantities{
  real pp[Ntotal];
  vector[Nvalid] log_lik;
  for(i in 1:Ntotal){
    pp[i] = normal_rng(S[i]+kappa*C[i],sigma);
    // Test equivalence with FIPS by removing noise
    //pp[i] = S[i]+kappa*C[i];
  }
  for(j in 1:Nvalid){
    log_lik[j] = normal_lpdf(fatigue[valid[j]] | S[valid[j]]+kappa*C[valid[j]],sigma);
  }
}
