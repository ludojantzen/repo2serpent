/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : schurfactorization.c                           */
/*                                                                           */
/* Created:       2014/06/08 (MPu)                                           */
/* Last modified: 2014/06/08 (MPu)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/*                                                                           */
/* Description: Complex Schur decomposition of complex matrix                */
/*                                                                           */
/* Comments:  Wilkinson shift                                                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "schurFactorization:"

/*****************************************************************************/
void schurFactorization(long n, complex **A, complex **T, complex **U)
{

  /* Schur factorization: A = U*T*U', T = upper triangular, U = unitary */
            
  long i,j,iter,maxIter;
  double tol, diff1,diff2; 
  complex T11, T12, T21, T22; 
  complex sigma1, sigma2, sigma; 
  complex z, z1, z2; 
  complex **P, **Q, **R;


  /* Allocate auxiliary matrices */

  P     = (complex **)Mem(MEM_ALLOC,n, sizeof(complex *)); 
  Q     = (complex **)Mem(MEM_ALLOC,n, sizeof(complex *)); 
  R     = (complex **)Mem(MEM_ALLOC,n, sizeof(complex *)); 

  for (i=0; i<n; i++){
    P[i]     = (complex *)Mem(MEM_ALLOC,n, sizeof(complex)); 
    Q[i]     = (complex *)Mem(MEM_ALLOC,n, sizeof(complex)); 
    R[i]     = (complex *)Mem(MEM_ALLOC,n, sizeof(complex)); 
  }

  /* ------------------------------------------------------------*/

  /* Parameters for iteration */

   maxIter = 500;
   tol     = 1E-30; 

  /* ------------------------------------------------------------*/

  /* Init U = eye(n) (identity matrix) */

  for (i=0; i<n; i++){
    U[i][i].re = 1.0; 
    U[i][i].im = 0.0; 
  }  
  
  /* ------------------------------------------------------------*/

  /* Reduce A to Hessenberg form */

  hessFactorization(n,A,P,T); 

  /* ------------------------------------------------------------*/

  /* Compute Schur factorization of Hessenberg matrix T */


   for (j=n-1; j>0; j--){ /* Main loop */

     for (iter=0; iter<maxIter; iter++){ /* Iteration loop */

       sigma.re = T[j][j].re;
       sigma.im = T[j][j].im; 


       /* -- Use Wilkinson shift -- */

       /* submatrix considered in the shift */

       T11 = T[j-1][j-1];
       T12 = T[j-1][j];
       T21 = T[j][j-1];
       T22 = T[j][j];

       /* Compute eigenvalues of submatrix */

       z.re  = 0.0;
       z.im  = 0.0;
       z2.re = 0.0;
       z2.im = 0.0;

       /* z = T11*T11 + T22*T22 - 2*T11*T22 + 4*T12*T21 */

       z1 = c_mul(T11,T11);

       z  = c_add(z ,z1);
       z2 = c_add(z2,z1);
       
       z1 = c_mul(T22,T22);

       z  = c_add(z ,z1);
       z2 = c_add(z2,z1);

       z1 = c_mul(T11,T22);

       z1.re = -2.0 * z1.re;
       z1.im = -2.0 * z1.im;

       z  = c_add(z,z1);

       z1 = c_mul(T12,T21);
       z1.re = 4.0 * z1.re;
       z1.im = 4.0 * z1.im;
       z = c_add(z,z1);

       /* Square root*/

       z = c_sqrt(z);

       /* Eigenvalues */
       
       sigma1 = c_add(z2,z);
       sigma2 = c_sub(z2,z);

/*        printf("sigma1 = %e %e\n", sigma1.re, sigma1.im); */
/*        printf("sigma2 = %e %e\n", sigma2.re, sigma2.im); */

       /* Select eigenvalue for shift*/

       diff1 = c_norm( c_sub(T[j][j], sigma1) );
       diff2 = c_norm( c_sub(T[j][j], sigma2) );

       if (diff1 < diff2){
         sigma.re = sigma1.re;
         sigma.im = sigma1.im;
       }else{
         sigma.re = sigma2.re;
         sigma.im = sigma2.im;
       }

       /* --- QR step with Wilkinson shift --- */

       /* Shift: T(1:j,1:j) = T(1:j,1:j) - sigma * eye(j) */

       for (i=0; i<j+1; i++){

         CheckValue(FUNCTION_NAME, "T[i][i].re","", T[i][i].re, -INFTY, INFTY);
         CheckValue(FUNCTION_NAME, "T[i][i].im","", T[i][i].im, -INFTY, INFTY);         

         T[i][i].re = T[i][i].re - sigma.re;   
         T[i][i].im = T[i][i].im - sigma.im;   
         
       }

       /* Compute QR factorization of shifted Hessenberg matrix */

       for (i=0; i<n; i++){
         memset(Q[i], 0, n*sizeof(complex));
         memset(R[i], 0, n*sizeof(complex));
       }

       QRfactorization(n,T,Q,R); 

       /* T = T_new = R * Q  */

       for (i=0; i<n; i++){
         memset(T[i], 0, n*sizeof(complex));
       }
       matProduct(n, n, n, R, Q, T);

       /* T(1:j,1:j) = T(1:j,1:j) + sigma * eye(j) */
       for (i=0; i<j+1; i++){
         T[i][i].re = T[i][i].re + sigma.re;   
         T[i][i].im = T[i][i].im + sigma.im;   
       }


       /* R =  U_new = U * Q */

       for (i=0; i<n; i++){
         memset(R[i], 0, n*sizeof(complex));
       }       
       matProduct(n,n,n,U,Q,R); 

       /* U = R */

       for (i=0; i<n; i++){
         memcpy(U[i],R[i], n*sizeof(complex));
       } 

       /* Check convergence */

       if (c_norm( T[j][j-1] ) <= tol * (c_norm(T[j-1][j-1]) + c_norm(T[j][j]))){
         T[j][j-1].re = 0.0;
         T[j][j-1].im = 0.0;
         break; 
       }
       
   
     }        /* end of iter loop */  
    
   } /* end of main loop */


  /* -------------------------------------------------------------*/

   /* U = P*U */

   for (i=0; i<n; i++){
     memset(U[i], 0, n*sizeof(complex));
   }
   matProduct(n,n,n,P,R,U);
   

  /* -------------------------------------------------------------*/
  /* Free auxiliary variables */

   for (i=0; i<n; i++){
     Mem(MEM_FREE,P[i]); 
     Mem(MEM_FREE,Q[i]); 
     Mem(MEM_FREE,R[i]); 
   }

   Mem(MEM_FREE,P); 
   Mem(MEM_FREE,Q); 
   Mem(MEM_FREE,R); 

  /* Return */

  return;   
    
  
}
/*--------------------------------------------------------------------*/

