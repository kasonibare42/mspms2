#ifndef _RANDOM_H_
#define _RANDOM_H_

/* Random number generator */
void rmarin(int ij, int kl);
void ranmar(double rvec[], int len);
/* Generate Gaussian random numbers
 From Numerical Recipes 
 */
double gaussran();

#endif
