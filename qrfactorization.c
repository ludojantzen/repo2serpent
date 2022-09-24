/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : qrfactorization.c                              */
/*                                                                           */
/* Created:       2014/06/07 (MPu)                                           */
/* Last modified: 2014/06/07 (MPu)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/*                                                                           */
/* Description: QR decomposition of complex matrix                           */
/*                                                                           */
/* Comments: based on Householder transformations                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "QRfactorization:"

/*****************************************************************************/

void QRfactorization(long n, complex **A, complex **Q, complex **R)
{

  /* QR factorization based on Householder transformations */
            
  long i,j,k,m;
  double nrm; 
  complex z,z1,z2;  
  complex *vj;

 /* Init. Q = eye(n) (identity matrix) */

  for (i=0; i<n; i++){
    Q[i][i].re = 1.0; 
    Q[i][i].im = 0.0; 
  }


  /* Init. R = A  */

  for(j=0; j<n; j++){
    for (i=0; i<n; i++){
      R[i][j].re = A[i][j].re;
      R[i][j].im = A[i][j].im;

      CheckValue(FUNCTION_NAME, "A[i][j].re","", A[i][j].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "A[i][j].im","", A[i][j].im, -INFTY, INFTY);

    }
  }

  /* Allocate auxiliary variables*/

  vj = (complex *)Mem(MEM_ALLOC, n, sizeof(complex));

 /*  printf("begin calc, n = %d\n", n);  */

    /* ------------------------------------------------------------*/

  for (j=0; j<n; j++){ /* Main loop */

    /* R(j:end, j)  */

    for (i=j; i<n; i++){
      vj[i-j].re = R[i][j].re;
      vj[i-j].im = R[i][j].im;
    }
    
    nrm = vectorNorm(n-j, vj); 

    /* v(1) = v(1) + sign(R(j,j)) * norm(v) */

    vj[0].re = vj[0].re + R[j][j].re / c_norm(R[j][j]) * nrm; 
    vj[0].im = vj[0].im + R[j][j].im / c_norm(R[j][j]) * nrm;


    /* Update norm */

    nrm = vectorNorm(n-j, vj);  

    /* v = v./norm(v) */

    for (i=0; i<n-j; i++){
      vj[i].re = vj[i].re / nrm; 
      vj[i].im = vj[i].im / nrm; 
    }

    /* Update */

    /* R(j:end, :) = R(j:end,:) - 2 * vj * vj' * R(j:end,:), : */

    /* Q(:,j:end)  = Q(:,j:end) - 2 * Q(:,j:end) * vj * vj^T */

    for (k=0; k<n; k++){

        /* (v * v' * A)_ik = v_i SUM_m Conj(v_m) A_mk */

        z.re = 0.0; 
        z.im = 0.0; 

        for (m=j; m<n; m++){

          z1 = c_con(vj[m-j]); 
          z1 = c_mul(z1, R[m][k]);           
          z  = c_add(z,z1); 
        }

        for (i=j; i<n; i++){ 

          z2    = c_mul(vj[i-j],z);

          /* Update R(i,k) */
          
          R[i][k].re = R[i][k].re - 2.0 * z2.re; 
          R[i][k].im = R[i][k].im - 2.0 * z2.im; 

          CheckValue(FUNCTION_NAME, "R[i][k].re","", R[i][k].re, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "R[i][k].im","", R[i][k].im, -INFTY, INFTY);

        } 

        /* (A * v * v^')_ki = v_i * SUM_m Conj(v_m) A_km */

        z.re = 0.0;
        z.im = 0.0;

        for (m=j; m<n; m++){

          z1 = vj[m-j]; 
          z1 = c_mul(z1, Q[k][m]);
          z  = c_add(z,z1);
        }

        for (i=j; i<n; i++){

          z1 = c_con(vj[i-j]); 
          z2 = c_mul(z1,z);

          /* Update Q(k,i)*/

          Q[k][i].re = Q[k][i].re - 2.0 * z2.re; 
          Q[k][i].im = Q[k][i].im - 2.0 * z2.im; 

          CheckValue(FUNCTION_NAME, "Q[k][i].re","", Q[k][i].re, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "Q[k][i].im","", Q[k][i].im, -INFTY, INFTY);
        }              
    }
  } /* End of main loop (j) */

  /* -------------------------------------------------------------*/

  /* Free auxiliary variables */

  Mem(MEM_FREE, vj); 
  return;   

    
  
}
/*--------------------------------------------------------------------*/


