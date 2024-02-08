data{
  int K;//number of units
  int NT;//number of times
  int tACF;//time of ACF
  matrix[K,NT] Y;   //notifications
}
transformed data{
  row_vector[NT] times;
  matrix[K,NT] TZ; //matrix of 1:NT in each row
  for(i in 1:NT) times[i] = 1.0*i;
  TZ = rep_matrix(times,NT);
}
parameters{
  // MVN stuff
  matrix[K,3] z;
  cholesky_factor_corr[3] L_Omega;
  vector<lower=0,upper=pi()/2>[3] tau_unif;  // prior scale
  row_vector[3] m;                        // means
  vector<lower=0>[K] cvec;                 //intercept vector
  vector[K] svec;                 //slope vector
  real<lower=0> bigeps; //noise prior SD
  vector[K] eps;       //noises
}
transformed parameters{
  matrix[K,NT] M;                 //mean predicted response
  matrix[K,NT] cpts;              //intercept matrix
  matrix[K,NT] slopes;            //slope matrix
  matrix[K,3] ml = rep_matrix(m,K);     // global means of log effects duplicated over units
  //MVN stuff
  vector<lower=0>[3] tau = 2.5 * tan(tau_unif);
  matrix[K,3] B = ml +  z * diag_pre_multiply(tau, L_Omega)';
  matrix[K,3] RRs = exp(B);       //effects in real space
  // for outputting cors
  matrix[3, 3] Sigma_beta = diag_pre_multiply(tau, L_Omega) * diag_pre_multiply(tau, L_Omega)';
  //regression
  cpts = rep_matrix(cvec,NT);
  slopes = rep_matrix(svec,NT);   
  //multiply trajectories after tACF by a RR factor & capture slope steepening
  // (indices starting at 1)
  // want mint[j+tACF] = A*m[tACF+1] + B*s*(j-1)   //[j is time after tACF]
  // = A*(c+s*tACF) + B*s*(j-1) = A*c + A*s*tACF + B*s*(j-1)
  // k = j+tACF 
  // mint[k] = A*c + A*s*tACF + B*s*(k-1-tACF) = A*c + (A-B)*s*tACF + B*s*(k-1)
  for(j in (tACF+1):NT){ //times after ACF
    cpts[:,j] = cpts[:,j] .* RRs[:,2] + (RRs[:,2]-RRs[:,3]).*svec*tACF;
    slopes[:,j] = slopes[:,j] .* RRs[:,3];
  }
  M = cpts + slopes .* TZ; 
  M[:,tACF] = M[:,tACF] .* RRs[:,1]; //ACF peak effect
}
model{
  //prior for MVN: see https://mc-stan.org/docs/2_26/stan-users-guide/multivariate-hierarchical-priors-section.html
  to_vector(z) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  //rest 
  cvec ~ cauchy(0,10); //intercept prior
  svec ~ normal(0,5); //slopes prior
  m ~ normal(0,5);//global mean log effects
  bigeps ~ normal(0,100);//global prior for unit noise SD
  eps ~ normal(0,bigeps);//noise-level specific to unit
  // data likelihood
  for(i in 1:K){
    //beta[i] ~ multi_normal(m,s); //log effects in unit i ~ some global distribution with correlation
    Y[i] ~ normal(M[i],eps[i]);//time-series in unit i
  }
}

