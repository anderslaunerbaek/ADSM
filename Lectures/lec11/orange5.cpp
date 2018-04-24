#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(ntimes);
  DATA_VECTOR(y);                                 // Data vector transmitted from R
  DATA_VECTOR(t);                                 // Data vector transmitted from R
  DATA_VECTOR(season);                                 // Data vector transmitted from R
  DATA_FACTOR(tree);                              // Data vector transmitted from R

  PARAMETER_VECTOR(u);                             // Random effects
  // Parameters
  PARAMETER_VECTOR(beta);                          // Parameter value transmitted from R
  PARAMETER(sigma_u);                              // Parameter value transmitted from R
  PARAMETER(sigma);                                // Parameter value transmitted from R
  PARAMETER(phi);

   using namespace density;
  int nobs = y.size();
  int ntrees = u.size();
  Type mean_ran = Type(0);
  vector<Type> pred(nobs);
  vector<Type> res(nobs);

  int j;

  Type f = 0;                                      // Declare the "objective function" (neg. log. likelihood)

  for(int j=0; j < ntrees; j++){
    f -= dnorm(u[j], beta[0] , sigma_u, true);
  }

  for(int i =0; i < nobs; i++){
    j = tree[i];
    pred[i] = u[j]/(1+exp(-((t[i]-beta[1])/beta[2] + season[i]*beta[3])));
    res[i] = y[i]-pred[i];
  }

  matrix<Type> cov(nobs,nobs);
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
