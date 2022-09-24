/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : parlett.c                                      */
/*                                                                           */
/* Created:       2014/06/08 (MPu)                                           */
/* Last modified: 2014/06/08 (MPu)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/*                                                                           */
/* Description: Compute matrix function based on Parlett method              */
/*              Input: Schur decomposition of the matrix                     */
/*                                                                           */
/* Comments:  -                                                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "parlett:"

/*****************************************************************************/
void parlett(long n, complex *d, complex **T, complex **U, complex **F)
{

  /* Compute matrix function F(A) based on Parlett algorithm */
  /* d = F(diag(A)) => this determines the function!!        */
  /* A = U * T * U' Schur factorization */
  /* T  = upper triangular, U = unitary   */          

  long i,j,ii,k;
  double eps;  
  complex diff,z,z1,z2; 
  complex **B; 

  eps = 1e-30; 

  /* ---------------------------------------------------------*/

  /* Allocate auxiliary variables */
  B = (complex **)Mem(MEM_ALLOC, n, sizeof(complex *)); 
  for (i=0; i<n; i++){
    B[i] = (complex *)Mem(MEM_ALLOC, n,sizeof(complex)); 
  }
  
  /* ---------------------------------------------------------*/

  
  /* Init. F:  */

  for (j=0; j<n; j++){ 

    F[j][j].re = d[j].re; 
    F[j][j].im = d[j].im; 
        
    CheckValue(FUNCTION_NAME, "d[j].re","", d[j].re, -INFTY, INFTY);
    CheckValue(FUNCTION_NAME, "d[j].im","", d[j].im, -INFTY, INFTY);

  }

  /* ---------------------------------------------------------*/

  /* Compute F(T) one super-diagonal at a time */


  for ( ii = 1; ii<n; ii++ ){ /* super-diagonals */


    /* diagonal: (j,j+i), j=1, n-i */

    for ( i = 0; i<n-ii; i++ ){   /* rows */

      j = i + ii; /*  */

      diff = c_sub(T[i][i], T[j][j]);

      if (c_norm( diff ) > eps){ /* OK */

        z1 = c_sub(F[i][i], F[j][j]);
        z2 = c_mul(z1,T[i][j]);
        z  = c_div(z2,diff);

        F[i][j].re = z.re;
        F[i][j].im = z.im;

        CheckValue(FUNCTION_NAME, "F[i][j].re","", F[i][j].re, -INFTY, INFTY);
        CheckValue(FUNCTION_NAME, "F[i][j].im","", F[i][j].im, -INFTY, INFTY);

        for (k=i+1; k<j; k++){

          z1 = c_mul(F[i][k], T[k][j]);
          z2 = c_mul(T[i][k], F[k][j]);

          z = c_sub(z1,z2);
          z = c_div(z,diff);
          
          F[i][j] = c_add(F[i][j],z);

          CheckValue(FUNCTION_NAME, "F[i][j].re","", F[i][j].re, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "F[i][j].im","", F[i][j].im, -INFTY, INFTY);
        }
      }
      else{
        Die(FUNCTION_NAME, "Multiple eigenvalue, division by zero");
        return;
      }
    }      
  } 
 
  /* ---------------------------------------------------------*/

  /* Compute F(A) = U * F(T) * U' */

  
  /* B = F * U' */

  for (i=0; i<n; i++){
    for (j=0; j<n; j++){
      for (k=i; k<n; k++){
        z1 = c_con(U[j][k]); /* transpose: j->k */
        z  = c_mul( F[i][k], z1);  
        B[i][j] = c_add(B[i][j],z); 

        CheckValue(FUNCTION_NAME, "B[i][j].re","", B[i][j].re, -INFTY, INFTY);
        CheckValue(FUNCTION_NAME, "B[i][j].im","", B[i][j].im, -INFTY, INFTY);

      }
    }
  }

  /* Clear F */
  for (i=0; i<n; i++){
    memset(F[i],0,n*sizeof(complex)); 
  }

  /* F = U * B */
  
  matProduct(n,n,n,U,B,F); 

  /* ---------------------------------------------------------*/

  /* Free auxiliary variables */

  for (i=0; i<n; i++){
    Mem(MEM_FREE, B[i]);     
  }
  Mem(MEM_FREE, B); 
  

  /* Return */

  return;   
    
  
}

/*--------------------------------------------------------------------*/
