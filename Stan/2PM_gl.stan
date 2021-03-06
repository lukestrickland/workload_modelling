
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
             real tau_d, //tau_d = controls rate of decay in S during sleep
             real L0
             ){ 
                
    real d=exp(-tas/tau_d);            
    real Sp=L0-d*(L0-ss);
    return Sp;
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
}

data {
  int<lower=0> Nsubj;
  int<lower=0> Ntotal;
  int<lower=0> subject[Ntotal];
  int<lower=0> observation[Ntotal];
  int<lower=0> day[Ntotal];
  vector<lower=0>[Ntotal] timeawake;
  vector<lower=0>[Ntotal] timeasleep;
  vector<lower=0>[Ntotal] totaltimeawake;
  vector<lower=0>[Ntotal] timeofday;
  vector<lower=0>[Ntotal] fatigue;
  int<lower=0> Nvalid;
  int<lower=0> valid[Nvalid];
}  
 
parameters {
  real<lower=0> S0;
  real<lower=0> U0;
  real<lower=0> L0;
  real phi;
  real<lower=0> kappa;
  real<lower=0> tau_d; //typically fixed to 4.2  h
  real<lower=0> tau_r; //typically fixed to 18.2 h
 // real<lower=0> A;     //typically fixed to 0.12
  real<lower=0> sigma; //typically not modeled, but necessary for parameter estimation
}

transformed parameters {

  //Fixed parameters
  //real t0 = -0.6;
  real tau = 24;
  real A = 1;
    
  //Calculate level of processes for each observation
  real sw;  //level of homeostatic process at time of waking
  real ss; //level of homeostatic process at time of sleep
  vector[Ntotal] S; //level of homeostatic process
  vector[Ntotal] C; //level of 24-hour circadian process
  
  for(i in 1:Ntotal){
    if(observation[i]==1){ //new wake episode: update S on wake
      if(day[i]==1){ //if first day, S on wake takes a value of 14 (as is convention)
        sw=S0;
      } 
      if(day[i]<1){//if after first day, S wake calculated based on timeasleep and S at sleep
        sw=Spfun(ss,timeasleep[i],tau_d,L0);
      }
      ss=Sfun(sw,totaltimeawake[i],tau_r,U0); //update Ssleep based on the time they went to bed that day;
    }
    
    S[i]=Sfun(sw,timeawake[i],tau_r,U0);
    C[i]=Cfun(timeofday[i],phi,tau,A);
  }
}

model {
  S0 ~ normal(0,10);
  U0 ~ normal(0,10);
  L0 ~ normal(0,10);
  phi ~ normal(0,5);
  kappa ~ normal(0,5);
  tau_d ~ normal(0,10); //typically fixed to 4.2  h
  tau_r ~ normal(0,10); //typically fixed to 18.2 h
  //A ~ normal(0,1);     //typically fixed to 0.12
  sigma ~ normal(0,3);

   for(i in 1:Nvalid){
    fatigue[valid[i]] ~ normal(S[valid[i]]+kappa*C[valid[i]],sigma);
  } 
}

generated quantities{
  real pp[Ntotal];
  vector[Nvalid] log_lik;
  for(i in 1:Ntotal){
    pp[i] = normal_rng(S[i]+kappa*C[i],sigma);
  }
  for(j in 1:Nvalid){
    log_lik[j] = normal_lpdf(fatigue[valid[j]] | S[valid[j]]+kappa*C[valid[j]],sigma);
  }
}
