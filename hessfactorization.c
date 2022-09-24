/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : hessFactorization.c                            */
/*                                                                           */
/* Created:       2014/06/07 (MPu)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Hessenberg decomposition of complex matrix                   */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "hessFactorization:"

/*****************************************************************************/
void hessFactorization(long n, complex **A, complex **P, complex **H)
{

  /* Hessenberg factorization based on Householder transformations */
  /* H = Hessenberg, P = unitary */
            
  long i,j,k,m;
  double nrm; 
  complex z,z1,z2;  
  complex *vj; 

  /* Allocate auxiliary variables */

  vj = (complex *)Mem(MEM_ALLOC, n,sizeof(complex));


 /* Init. P = eye(n) (identity matrix) */

  for (i=0; i<n; i++){
    P[i][i].re = 1.0; 
    P[i][i].im = 0.0; 
  }

  /* Init. H = A  */

  for(j=0; j<n; j++){
    for (i=0; i<n; i++){
      H[i][j].re = A[i][j].re;
      H[i][j].im = A[i][j].im;

      CheckValue(FUNCTION_NAME, "A[i][j].re","", A[i][j].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "A[i][j].im","", A[i][j].im, -INFTY, INFTY);

    }
  }

  /* ------------------------------------------------------------*/

  for (j=0; j<n-2; j++){ /* Main loop */

    /* H(j+1:end, j)  */

    for (i=j+1; i<n; i++){
      vj[i-j-1].re = H[i][j].re;
      vj[i-j-1].im = H[i][j].im;
    }

    nrm = vectorNorm(n-j-1, vj); 

    /* v(1) = v(1) + sign(H(j+1,j)) * norm(v) */

    /* Check division by zero (30.9.2016 / 2.1.28 / JLe) */

    if ((fabs(H[j+1][j].re) > ZERO) || (fabs(H[j+1][j].im) > ZERO))
      {
        vj[0].re = vj[0].re + H[j+1][j].re / c_norm(H[j+1][j]) * nrm; 
        vj[0].im = vj[0].im + H[j+1][j].im / c_norm(H[j+1][j]) * nrm;
      }

#ifdef DEBUG

    else
      Warn(FUNCTION_NAME, "Division by zero");

#endif

    /* Update norm */

    nrm = vectorNorm(n-j-1, vj);  

    /* v = v./norm(v) */

    for (i=0; i<n-j-1; i++){
      vj[i].re = vj[i].re / nrm; 
      vj[i].im = vj[i].im / nrm; 
    }

    /* Update */

    /* H(j+1:end,j:end) = H(j+1:end,j:end) - 2 * vj * vj' * H(j+1:end,j:end), : */

    for (k=j; k<n; k++){

        /* (v * v' * A)_ik = v_i SUM_m Conj(v_m) A_mk */

        z.re = 0.0; 
        z.im = 0.0; 

        for (m=j+1; m<n; m++){

          z1 = c_con(vj[m-j-1]); 
          z1 = c_mul(z1, H[m][k]);           
          z  = c_add(z,z1); 
        }

        for (i=j+1; i<n; i++){ 

          z2    = c_mul(vj[i-j-1],z);

          /* Update H(i,k) */
          
          H[i][k].re = H[i][k].re - 2.0 * z2.re; 
          H[i][k].im = H[i][k].im - 2.0 * z2.im; 

        } 
    }
    
    /* H(:, j+1:end) = H(:,j+1:end) - 2 * H(:,j+1:end) * v * v';  */

    for (k=0; k<n; k++){

        /* (A * v * v^')_ki = Conj(v_i) * SUM_m v_m A_km */

        z.re = 0.0;
        z.im = 0.0;

        for (m=j+1; m<n; m++){

          z1 = vj[m-j-1];
          z1 = c_mul(z1, H[k][m]);
          z  = c_add(z,z1);
        }

        for (i=j+1; i<n; i++){

          z1 = c_con(vj[i-j-1]);
          z2 = c_mul(z1,z);

          /* Update Q(k,i)*/

          H[k][i].re = H[k][i].re - 2.0 * z2.re;
          H[k][i].im = H[k][i].im - 2.0 * z2.im;
        }
    }

    /* P(:,j+1:end)  = P(:,j+1:end) - 2 * P(:,j+1:end) * vj * vj' */

    for (k=0; k<n; k++){        

        /* (A * v * v^')_ki = Conj(v_i) * SUM_m v_m A_km */

        z.re = 0.0;
        z.im = 0.0;

        for (m=j+1; m<n; m++){

          z1 = vj[m-j-1]; 
          z1 = c_mul(z1, P[k][m]);
          z  = c_add(z,z1);
        }

        for (i=j+1; i<n; i++){

          z1 = c_con(vj[i-j-1]); 
          z2 = c_mul(z1,z);

          /* Update P(k,i)*/

          P[k][i].re = P[k][i].re - 2.0 * z2.re; 
          P[k][i].im = P[k][i].im - 2.0 * z2.im; 
        }
    }
              
  } /* end of main loop */


  /* -------------------------------------------------------------*/

  /* Free auxiliary variables */

  Mem(MEM_FREE, vj); 

  /* Return */

  return;   
    
  
}
/*--------------------------------------------------------------------*/




