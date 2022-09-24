/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : bsfunn.c                                       */
/*                                                                           */
/* Created:       2014/06/07 (MPu)                                           */
/* Last modified: 2016/03/01 (JLe)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Basis functions for homogeneous diffusion flux solver        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "BsfunN:"

/*****************************************************************************/

void BsfunN(long nG, double *v, double *r, complex **T, complex **U, 
            complex **E)
{
  /* Compute normal derivative of DE basis function at r */
  /* function: expm( v'*r sqrtm(A) )                     */
  /* A = U*T*U', T tirangular, U unitary                 */

  long i, j, dim; 
  double prod; 
  complex z; 

  complex *d, **As,**B; 


  /* Allocate auxiliary variables */

  d  = (complex * )Mem(MEM_ALLOC, nG, sizeof(complex  )); 
  As = (complex **)Mem(MEM_ALLOC, nG, sizeof(complex *)); 
  B  = (complex **)Mem(MEM_ALLOC, nG, sizeof(complex *)); 
  for (i=0; i<nG; i++){
    As[i] = (complex *)Mem(MEM_ALLOC, nG, sizeof(complex)); 
    B [i] = (complex *)Mem(MEM_ALLOC, nG, sizeof(complex)); 
  }

  /*-----------------------------------------------------------------------*/


  /* 3-D space */

  dim = 3; 

  /* Compute inner product v' * r:  */ 

  prod = 0.0; 
  for (i=0; i<dim; i++){
    prod = prod + v[i] * r[i]; 
  }

  /*-----------------------------------------------------------------------*/

  /* Use Parlett method to compute e^( prod * sqrt( A ) ) */


  /* Compute B = expm (prod * sqrtm (A)) */

  if (prod > 1E-20 || prod <-1E-20){ /* Compute normal derivative on the boundary */

    for (i=0; i<nG; i++){
      z    = T[i][i]; 
      z    = c_sqrt(z); 
      z.re = z.re * prod; 
      z.im = z.im * prod; 

      CheckValue(FUNCTION_NAME, "z.re","", z.re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "z.im","", z.im, -INFTY, INFTY);
      
      d[i] = c_exp(z);     
    }

  /* Parlett algorithm */

    parlett(nG, d, T, U, B); 

  }
  else{ /* expm of zero  matrix => identity matrix */
  
    for (i=0; i<nG; i++){
      B[i][i].re = 1.0; 
    }    
  }

  for (i=0; i<nG; i++){
    for (j=0; j<nG; j++){
      CheckValue(FUNCTION_NAME, "B[i][j].re","", B[i][j].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "B[i][j].im","", B[i][j].im, -INFTY, INFTY);
    }
  }


  /* Compute As = sqrtm( A ) */

  for (i=0; i<nG; i++){
    z    = T[i][i]; 
    d[i] = c_sqrt(z); 
  }

   
  parlett(nG, d, T, U, As); 

  for (i=0; i<nG; i++){
    for (j=0; j<nG; j++){
      CheckValue(FUNCTION_NAME, "As[i][j].re","", As[i][j].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "As[i][j].im","", As[i][j].im, -INFTY, INFTY);
    }
  }

  matProduct(nG,nG,nG,As,B,E); 

  for (i=0; i<nG; i++){
    for (j=0; j<nG; j++){
      CheckValue(FUNCTION_NAME, "E[i][j].re","", E[i][j].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "E[i][j].im","", E[i][j].im, -INFTY, INFTY);
    }
  }

  /* Free auxiliary variables */  

  Mem(MEM_FREE, d); 
  for (i=0; i<nG; i++){
    Mem(MEM_FREE,As[i]);
    Mem(MEM_FREE,B[i]);      
  }

  Mem(MEM_FREE,As); 
  Mem(MEM_FREE,B); 

  /* Return */

  return; 

}

/**************************************************************************************/
