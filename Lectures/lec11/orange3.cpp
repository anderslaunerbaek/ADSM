#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);                                 // Data vector transmitted from R
  DATA_VECTOR(t);                                 // Data vector transmitted from R
  DATA_VECTOR(season);                                 // Data vector transmitted from R
  DATA_FACTOR(tree);                              // Data vector transmitted from R

  PARAMETER_VECTOR(u);                             // Random effects
  // Parameters
  PARAMETER_VECTOR(beta);                          // Parameter value transmitted from R
  PARAMETER(sigma_u);                              // Parameter value transmitted from R
  PARAMETER(sigma);                                // Parameter value transmitted from R

  int nobs = y.size();
  int ntrees = u.size();
  Type mean_ran = Type(0);

  int j;

  Type f = 0;                                      // Declare the "objective function" (neg. log. likelihood)

  for(int j=0; j < ntrees; j++){
    f -= dnorm(u[j], beta[0] , sigma_u, true);
  }

  for(int i =0; i < nobs; i++){
    j = tree[i];
    f -= dnorm(y[i], u[j]/(1+exp(-((t[i]-beta[1])/beta[2] +
				   season[i]*beta[3]))), sigma, true);
  }

  return f;
}
