/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dfsol.c                                        */
/*                                                                           */
/* Created:       2014/06/07 (MPu)                                           */
/* Last modified: 2014/06/07 (MPu)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/*                                                                           */
/* Description: Solution to homogeneous diffusion eq. at r                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "dfSol:"

/*****************************************************************************/

double *dfSol(long nG, long nJ, double **surfs_n, double *r, 
              complex **T, complex **U, complex *c)
{
  /* Return Phi = solution to hom. DE at r */

  long i, j, ii, jj, dim; 
  double prod; 
  double *phi; 
  complex z; 
  complex *d, **E, **M; 

  /* Allocate Phi */

  phi = (double *)Mem(MEM_ALLOC, nG, sizeof(double)); 

  /* Allocate temporary variables */

  d    = (complex * )Mem(MEM_ALLOC, nG, sizeof(complex  )); 

  E    = (complex **)Mem(MEM_ALLOC, nG, sizeof(complex *));
  for(i=0; i<nG; i++)
    E [i] = (complex *)Mem(MEM_ALLOC, nG, sizeof(complex)); 


  /* M = nG x (nG * nJ)  */
  M    = (complex **) Mem(MEM_ALLOC, nG   ,sizeof(complex *));
  for(i=0; i<nG; i++)
    M[i] = (complex *)Mem(MEM_ALLOC, nG*nJ,sizeof(complex)); 

  /*-----------------------------------------------------------------------*/

  for(j=0; j<nJ; j++){ /* One basis function at a time */

    /* 3-D space */

    dim = 3;
    
    /*Compute inner product n' * r: */

    for (ii=0; ii<nG; ii++)
      memset(E[ii], 0, nG*sizeof(complex)); 

    prod = 0.0; 
    for (ii=0; ii<dim; ii++){
      prod = prod + surfs_n[j][ii] * r[ii]; 
      CheckValue(FUNCTION_NAME, "prod","", prod, -INFTY, INFTY); 
    }

    if ( prod > 1E-20 || prod <-1E-20 ){ 

      /* E = expm( prod * sqrtm(A) ), A = U*T*U' */ 


      for (ii=0; ii<nG; ii++){

        z    = T[ii][ii];
        z    = c_sqrt(z);
        z.re = z.re * prod;
        z.im = z.im * prod;

        CheckValue(FUNCTION_NAME, "z.re","", z.re, -INFTY, INFTY);
        CheckValue(FUNCTION_NAME, "z.im","", z.im, -INFTY, INFTY);
      
        d[ii] = c_exp(z);
      }

      /* Parlett algorithm */

      parlett(nG, d, T, U, E);

    }
    else{ /* expm of zero  matrix => identity matrix */
  
      for (ii=0; ii<nG; ii++){
        E[ii][ii].re = 1.0;
      }
    }

    /* Check matrix exponential : */

    for (ii=0; ii<nG; ii++){
      for (jj=0; jj<nG; jj++){
        CheckValue(FUNCTION_NAME, "E[ii][jj].re","", E[ii][jj].re, -INFTY, INFTY);
        CheckValue(FUNCTION_NAME, "E[ii][jj].im","", E[ii][jj].im, -INFTY, INFTY);
/*         printf("E(%ld,%ld) = %e %e\n", ii+1, jj+1, E[ii][jj].re, E[ii][jj].im);  */
      }
    }

    /* Store to matrix M */

    for (ii = 0; ii<nG; ii++){
      for (jj=0; jj<nG; jj++){
        M[ii][j*nG + jj] = E[ii][jj];
      }
    }
   
  } /* End of j loop */

 /*-----------------------------------------------------------------------*/

  /* Phi = M * c  */


    for (i=0; i<nG; i++){
      for (j=0; j<nG*nJ; j++){
        z = c_mul(M[i][j],c[j]);
        phi[i] += z.re;
    }
  }

  /* OK! */

 /*-----------------------------------------------------------------------*/

  /* Free auxiliary variables */  

  Mem(MEM_FREE, d); 

  for (i=0; i<nG; i++){
    Mem(MEM_FREE,E[i]);
  }
  for (i=0; i<nG; i++){
    Mem(MEM_FREE,M[i]);
  }

  Mem(MEM_FREE,E); 
  Mem(MEM_FREE,M); 

  /* Return */

  return phi; 

}

/*****************************************************************************/


