#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);                                 // Data vector transmitted from R
  DATA_VECTOR(n);                                 // Data vector transmitted from R
  DATA_VECTOR(logc);                              // Data vector transmitted from R

  PARAMETER_VECTOR(B);                             // Random effects
  // Parameters
  PARAMETER_VECTOR(beta);                          // Parameter value transmitted from R
  PARAMETER(sigma_b);                              // Parameter value transmitted from R
  
  int nobs = y.size();
  Type mean_ran = Type(0);
  Type p = Type(0);

  Type f = 0;                                      // Declare the "objective function" (neg. log. likelihood)

  for(int i =0; i < nobs; i++){
    f -= dnorm(B[i], mean_ran , sigma_b, true);
    p = 1/(1+exp(-beta[0]-beta[1]*logc[i]-B[i]));
    f -= dbinom(y[i],n[i],p, true);
  }

  return f;
}
