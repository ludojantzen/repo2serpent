/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ludecomposition.c                              */
/*                                                                           */
/* Created:       2014/06/08 (MPu)                                           */
/* Last modified: 2014/06/08 (MPu)                                           */
/* Version:       2.1.22                                                      */
/*                                                                           */
/*                                                                           */
/* Description: Solve Ax = b and compute LU decomposition of A               */
/*                                                                           */
/* Comments: This is for dense matrices                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "LUdecomposition:"

/*****************************************************************************/
void LUdecomposition(long n, complex **A, complex *b, complex *x)
{

  long i,j,k,*p;
  complex akk, bk, z, sum; 
  complex *Ak; 


  /* Allocate auxiliary variables */

  p = (long *)Mem(MEM_ALLOC, n, sizeof(long)); 


  /* Gaussian elimination with partial pivoting */

  for (k=0; k<n-1; k++){ /* columns of A */


    /* --- Pivoting --- */

    akk.re = 0.0; 
    akk.im = 0.0; 
    j = k; 
    for (i=k; i<n; i++){
      if ( c_norm(A[i][k]) > c_norm(akk)){
        akk = A[i][k]; 

        CheckValue(FUNCTION_NAME, "A[i][k].re","", A[i][k].re, -INFTY, INFTY);
        CheckValue(FUNCTION_NAME, "A[i][k].im","", A[i][k].im, -INFTY, INFTY);

        j = i; 
      }
    }

    if (j != k){ /* swap rows j and k*/


      Ak   = A[k];  
      A[k] = A[j]; /* swap pointers*/
      A[j] = Ak; 

      bk   = b[k]; 
      b[k] = b[j]; 
      b[j] = bk; 

    }

    p[k] = j;  /* Keep list of permutations */

    /* Sanity check */

    if (akk.re != A[k][k].re && akk.im != A[k][k].im){
      Die(FUNCTION_NAME, "Something went wrong with pivoting?"); 
      return; 
    }

    if (akk.re == 0.0 && akk.im == 0.0){
      Die(FUNCTION_NAME, "Matrix singular?"); 
      return; 
    }

    /* Store L part */

    for (i=k+1;i<n; i++){
      A[i][k] = c_div(A[i][k],A[k][k]); 

      CheckValue(FUNCTION_NAME, "A[i][k].re","", A[i][k].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "A[i][k].im","", A[i][k].im, -INFTY, INFTY);

    }

    /* Update */

    for (i=k+1; i<n; i++){

      /* update b */

      z = c_mul(A[i][k], b[k]);  
      b[i] = c_sub(b[i], z); 

      CheckValue(FUNCTION_NAME, "b[i].re","", b[i].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "b[i].im","", b[i].im, -INFTY, INFTY);

      for (j=k+1;  j<n; j++){
        
        /* Update U part of A */
        
        z = c_mul( A[i][k], A[k][j]);
        A[i][j] = c_sub(A[i][j],z);  

        CheckValue(FUNCTION_NAME, "A[i][j].re","", A[i][j].re, -INFTY, INFTY);
        CheckValue(FUNCTION_NAME, "A[i][j].im","", A[i][j].im, -INFTY, INFTY);
      }
    }        
  }


  /* Solve x */

  x[n-1] = c_div(b[n-1],A[n-1][n-1]);

  for (i=n-2; i>=0; i--){
    sum.re = 0.0;
    sum.im = 0.0;
    for (j=i+1; j<n; j++){
      z   = c_mul(A[i][j], x[j]);
      sum = c_add(sum, z);
    }
    z    = c_sub(b[i], sum);
    x[i] = c_div(z, A[i][i]);
  }

  Mem(MEM_FREE, p); 

  /* Return */

  return;   
    
  
}
/*--------------------------------------------------------------------*/

