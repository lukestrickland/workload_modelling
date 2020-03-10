
functions{
  //Equations from Rajdev (2013). A unified mathematical model to quantify performance impairment for both chronic sleep restriction and total sleep deprivati. Journal of Theoretical Biology
  
  //Note also that Rajdev (2013) Theoretical Biology has
  //some good info on estimating the model. They estimate   //S0, U (upper assymptote), L (lower assymptote), tau_w, tau_s, phi, kappa. Phi represents the phase (t0 I think), and kappa represents the relative influence of C over S.
  
  //Higher S and C indicate greater performance impairment
  
  //calculates S during wake
  real Sfun(real sw, //sw = S upon waking
             real taw, //taw = time awake
             real tau_w, //tau_r = controls rate of rise in S during wake
             real U0     //upper assymptote
             ){ 
    real r=exp(-taw/tau_w);        
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
  
  real Lfun(real ls, //lower assymptote upon falling asleep
           real tas,    //time asleep
           real tau_la, //rate of change in lower assymptote
           real U0      //upper assymptote
          ){
    real term1 = ls*exp(-tas/tau_la);
    real term2 = -2*U0*(1-exp(-tas/tau_la));
    real L=term1+term2;
    return(L);
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
  real<lower=-2,upper=1> L0;  //Sleep debt. Needs to be between -2 and 1 to start off.
  real<lower=1> tau_s; 
  real<lower=0,upper=40> tau_w; //typically fixed to 18.2 h
  real<lower=0> phi;
  real<lower=0> kappa; //influence of circadian process - represents A in this model
  real<lower=0> tau_la;
  real<lower=0> sigma; //typically not modeled, but necessary for parameter estimation
}

transformed parameters {

  //Fixed parameters
  real A = 1;
  real tau = 24;
    
  //Calculate level of processes for each observation
  real sw;  //level of homeostatic process at time of waking
  real ss; //level of homeostatic process at time of sleep
  real lw;  //level of lower assymptote at time of waking
  real ls; //level of lower assymptote at time of sleep (identical to Lwake for a given day according to UM)
  vector[Ntotal] S; //level of homeostatic process
  vector[Ntotal] C; //level of 24-hour circadian process
  vector[Ntotal] L; //lower assumptote of homeostatic process
  
  for(i in 1:Ntotal){
    if(observation[i]==1){ //new wake episode: update S on wake
      if(day[i]==1){ //if first day, S on wake takes a value of 14 (as is convention)
        sw=S0;
        lw=L0;
      } 
      if(day[i]>1){//if after first day, S wake calculated based on timeasleep and S at sleep
        sw=Spfun(ss,timeasleep[i],tau_s,U0,tau_la,ls);
        lw=Lfun(ls,timeasleep[i],tau_la,U0);
      }
      ss=Sfun(sw,totaltimeawake[i],tau_w,U0); //update Ssleep based on the time they went to bed that day;
      ls=lw;
    }
    
    L[i]=lw;
    S[i]=Sfun(sw,timeawake[i],tau_w,U0);
    C[i]=Cfun(timeofday[i],phi,tau,A);
  }
}

model {
  S0 ~ normal(0,10); //originally estimated value = 0
  U0 ~ normal(10,10); //originally estimated value = 24.12
  L0 ~ normal(-0.5,1); //originally estimated value = 0
  tau_s ~ normal(0,10); //originally estimated value = 1
  tau_w ~ normal(10,10); //originally estimated value = 40
  tau_la ~ normal(0,10);  //originally estimated value = 4.06
  phi ~ normal(0,5);     //originally estimated value = 2.02
  kappa ~ normal(0,5);   //originally estimated value = 4.13
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
