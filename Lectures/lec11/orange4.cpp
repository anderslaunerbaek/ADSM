#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n);
  DATA_VECTOR(y);                                 // Data vector transmitted from R
  DATA_VECTOR(t);                                 // Data vector transmitted from R
  DATA_FACTOR(tree);                              // Data vector transmitted from R
  DATA_FACTOR(times);                              // Data vector transmitted from R


  PARAMETER_VECTOR(u1);                             // Random effects
  PARAMETER_VECTOR(u2);                             // Random effects
  // Parameters
  PARAMETER_VECTOR(beta);                          // Parameter value transmitted from R
  PARAMETER(sigma_u1);                              // Parameter value transmitted from R
  PARAMETER(sigma_u2);                              // Parameter value transmitted from R
  PARAMETER(sigma);                                // Parameter value transmitted from R
   PARAMETER(phi);                                // Parameter value transmitted from R

  using namespace density;
  int nobs = y.size();
  int ntrees = u1.size();
  int ntimes = u2.size();
  Type mean_ran = Type(0);
  vector<Type> pred(n);
  vector<Type> res(n);


  int j;
  int k;

  Type f = 0;                                      // Declare the "objective function" (neg. log. likelihood)

  for(int j=0; j < ntrees; j++){
    f -= dnorm(u1[j], mean_ran, sigma_u1, true);
  }

  for(int k=0; k < ntimes; k++){
    f -= dnorm(u2[k], mean_ran, sigma_u2, true);
  }

  for(int i =0; i < nobs; i++){
    j = tree[i];
    k = times[i];
    pred[i] = (beta[0]+u1[j]+u2[k])/(1+exp(-(t[i]-beta[1])/beta[2]));
    res[i] = y[i]-pred[i];
  }

  matrix<Type> cov(n,n);
  for (int i=0;i<ntrees;i++)
  {
    for (int j=0;j<ntimes;j++)
      {
	cov(i*ntimes+j,i*ntimes+j) = sigma * sigma;
	for (int k=0;k<j;k++)
	  {
	    cov(i*ntimes+k,i*ntimes+j) = sigma * sigma * exp(- phi * Type(2)/Type(365) * (t[j]-t[k]));
	    cov(i*ntimes+j,i*ntimes+k) = cov(i*ntimes+k,i*ntimes+j);
	  }
      }
  }  
  
  MVNORM_t<Type> neg_log_density(cov);

  f += neg_log_density(res);
  
  return f;
}

