#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);                                 // Data vector transmitted from R
  PARAMETER(mu);                                  // Parameter value transmitted from R
  PARAMETER(sigma);                               //                 

  Type f;                                         // Declare the "objective function" (neg. log. likelihood)
  f = -sum(dnorm(x,mu,sigma,true));               // Use R-style call to normal density

  return f;
}
