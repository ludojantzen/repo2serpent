/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dfsolver.c                                     */
/*                                                                           */
/* Created:       2014/06/07 (MPu)                                           */
/* Last modified: 2016/03/03 (JLe)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Solve homogeneous diffusion equation                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DFSolver:"

/*****************************************************************************/

void DFSolver(long ifCorn, long nG, long nJ, const double *ds0, double dc, 
              double **surfs_n, double **surfs_r, complex **T, complex **U, 
              complex *b, complex *c)
{
  long i, j, k, ii, jj, ind, ind1, ind2, Nr, nS;
  double r[3], tn[3], tn1[3], tn2[3]; 
  double *x,*y,*t, *x1, *x2, *y1, *y2, ds; 
  double del, nrm, nrm1, nrm2; 
  complex z; 
  complex **E, **Mr, **M; 

  /* Integration parameter */

  Nr = 100; 

  /* Number of boundary surfaces (basically checks if */
  /* corners are included) */

  if ((nJ == 8) || (nJ == 12))
    nS = nJ/2; 
  else
    nS = nJ; 
  
  /* Allocate temporary variables */

  x = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  x1 = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  x2 = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  y = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  y1 = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  y2 = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  t = (double *)Mem(MEM_ALLOC, Nr, sizeof(double));
  
  E  = (complex **)Mem(MEM_ALLOC, nG   , sizeof(complex));
  for (i = 0; i < nG; i++)
    E[i] = (complex *)Mem(MEM_ALLOC, nG, sizeof(complex)); 

  Mr = (complex **)Mem(MEM_ALLOC, nG*nG, sizeof(complex)); 
  for (i = 0; i < nG*nG; i++)
    Mr[i] = (complex *)Mem(MEM_ALLOC, Nr, sizeof(complex)); 

  M  = (complex **)Mem(MEM_ALLOC, nG*nJ, sizeof(complex));
  for (i = 0; i < nG*nJ; i++)
    M[i] = (complex *)Mem(MEM_ALLOC, nG*nJ, sizeof(complex)); 
  
  /***************************************************************************/

  /***** Surfaces ************************************************************/

  /* Loop over surfaces */

  for (i = 0; i < nS; i++)
    { 
      /* Surface area */

      ds = ds0[i];
      CheckValue(FUNCTION_NAME, "ds", "", ds, ZERO, INFTY);

      /* Integration variables */

      del = ds/((double)Nr - 1.0); 
      for(j = 1; j < Nr; j++)
        t[j] = t[j - 1] + del; 

      /* Tangent vector of the surface  */

      tn[0] = -surfs_n[i][1];
      tn[1] = surfs_n[i][0];
      tn[2] = 0.0; 

      /* Take Nr points from the boundary */
      
      for (k = 0; k < Nr; k++)
        {
          /* x(t) = x0 + tn_x * t */
          
          x[k] = surfs_r[i][0] + t[k]*tn[0]; 
          y[k] = surfs_r[i][1] + t[k]*tn[1]; 
        }

      /* Loop over basis functions */

      for (j = 0; j < nJ; j++)
        {
          /* surfs_n[j] = basis function vector */
          /* surfs_n[i] = surface normal vector */

          /* Compute inner product of vectors */

          nrm = -surfs_n[j][0]*surfs_n[i][0] - surfs_n[j][1]*surfs_n[i][1]; 

          /* Divide by surface area (this is for integration later) */
      
          nrm = nrm/ds;

          /* Integrate over surface i: (nrm = 0.0 => zero matrix) */

          if ((nrm > 1E-20) || (nrm < -1E-20))
            { 
              /* Loop over points */
              
              for (k = 0; k < Nr; k++)
                {
                    r[0] = x[k]; 
                  r[1] = y[k]; 
                  r[2] = 0.0; 

                  /* Normal derivative at r */

                  for (ii = 0; ii < nG; ii++)
                    memset(E[ii], 0, nG*sizeof(complex)); 
                  
                  BsfunN(nG, surfs_n[j], r, T, U, E); 

                  /* Store row by row */

                  ind = 0; 
                  for (ii = 0; ii < nG; ii++)
                    for (jj = 0; jj < nG; jj++)
                      {
                        Mr[ind][k] = E[ii][jj]; 
                        ind = ind + 1; 
                      }
                }  

              /* Integrate => E */

              ind = 0; 
              for (ii = 0; ii < nG; ii++)
                for (jj = 0; jj < nG; jj++)
                  {
                    z = trapz(Nr, t, Mr[ind]); 
                    E[ii][jj].re = z.re*nrm; 
                    E[ii][jj].im = z.im*nrm; 
                      ind = ind + 1; 
                  }
              
              /* Store to matrix M:  */

              /* -> functions */
              /* |           */
              /* V  surfaces */

              for (ii = 0; ii < nG; ii++)
                for (jj = 0; jj < nG; jj++)
                  M[i*nG + ii][j*nG + jj] = E[ii][jj]; 
            }
        }
    }

  /***************************************************************************/

  /***** Corners *************************************************************/

  /* Check if corners are included */

  if (ifCorn == 1)
    {
      /* Check corner area */

      CheckValue(FUNCTION_NAME, "dc", "", dc, ZERO, INFTY);

      /* Loop over corners (assume nC = nS) */

      for (i = 0; i < nS; i++)
        {
          /* Integration variables (half corner) */
          
          del = (dc/2.0)/((double)Nr - 1.0);
          for(j = 1; j < Nr; j++)
            t[j] = t[j - 1] + del; 
          
          /* First  part of corner on surface i, second part of */
          /* corner on surface i + 1 */

          ind1 = i; 

          /* Check for last corner */

          if (i == nS - 1)
            ind2 = 0; 
          else
            ind2 = i + 1; 

          /* Tangent vector for surface 1 */

          tn1[0] = -surfs_n[ind1][1];
          tn1[1] = surfs_n[ind1][0];
          tn1[2] = 0.0;

          /* Tangent vector for surface 2 */

          tn2[0] = -surfs_n[ind2][1];
          tn2[1] = surfs_n[ind2][0];
          tn2[2] = 0.0;

          /* Take Nr points  */

          for (k = 0; k < Nr; k++)
            {
              /*
              x1[k] = surfs_r[ind1][0] + (0.95*ds0[0] + t[k])*tn1[0]; 
              y1[k] = surfs_r[ind1][1] + (0.95*ds0[0] + t[k])*tn1[1]; 
              */

              x1[k] = surfs_r[ind2][0] - t[Nr - 1]*tn1[0] + t[k]*tn1[0]; 
              y1[k] = surfs_r[ind2][1] - t[Nr - 1]*tn1[1] + t[k]*tn1[1];

              x2[k] = surfs_r[ind2][0] + t[k]*tn2[0]; 
              y2[k] = surfs_r[ind2][1] + t[k]*tn2[1];
            }

          /* Loop over basis functions */

          for (j = 0; j < nJ; j++)
            {
              /* surfs_n[j] = basis function vector */
              /* surfs_n[i] = surface normal vector */

              /* Compute inner product of vectors */

              nrm1 = -surfs_n[j][0]*surfs_n[ind1][0] 
                - surfs_n[j][1]*surfs_n[ind1][1]; 

              nrm2 = -surfs_n[j][0]*surfs_n[ind2][0] 
                - surfs_n[j][1]*surfs_n[ind2][1] ; 

              /* Divide by surface area (this is for integration later) */
      
              nrm1 = nrm1/dc;
              nrm2 = nrm2/dc;

              /* Integrate first part (nrm1 = 0 => zero matrix) */
              
              if ((nrm1 > 1E-20) || (nrm1 < -1E-20))
                {
                  /* Loop over points */

                  for (k = 0; k < Nr; k++)
                    {
                      r[0] = x1[k]; 
                      r[1] = y1[k]; 
                      r[2] = 0.0;
                      
                      /* Normal derivative at r */

                      for (ii = 0; ii < nG; ii++)
                        memset(E[ii], 0, nG*sizeof(complex)); 

                      BsfunN(nG, surfs_n[j], r, T, U, E);

                      /* Store row by row */

                      ind = 0; 
                      for (ii = 0; ii < nG; ii++)
                        for (jj = 0; jj < nG; jj++)
                          {
                            Mr[ind][k] = E[ii][jj]; 
                            ind = ind + 1; 
                          }
                    }  

                  /* Integrate => E*/

                  ind = 0; 
                  for (ii = 0; ii < nG; ii++)
                    for (jj = 0; jj < nG; jj++)
                      {
                        z = trapz(Nr, t, Mr[ind]); 
                        E[ii][jj].re = z.re*nrm1; 
                        E[ii][jj].im = z.im*nrm1;   
                        ind = ind + 1; 
                      }
                  
                  /* Store to matrix M: */
     
                  for (ii = 0; ii < nG; ii++)
                    for (jj = 0; jj < nG; jj++)
                      M[nS*nG + i*nG + ii][j*nG + jj] = E[ii][jj]; 
                } 
              
              /* Integrate second part (nrm2 = 0 => zero matrix) */

              if ((nrm2 > 1E-20) || (nrm2 < -1E-20))
                {
                  /* Loop over points */
                  
                  for (k = 0; k < Nr; k++)
                    {
                      r[0] = x2[k]; 
                      r[1] = y2[k]; 
                      r[2] = 0.0;
                      
                      /* Normal derivative at r */

                      for (ii = 0; ii < nG; ii++)
                        memset(E[ii], 0, nG*sizeof(complex)); 

                      BsfunN(nG, surfs_n[j], r, T, U, E);

                      /* Store row by row */

                      ind = 0; 
                      for (ii = 0; ii < nG; ii++)
                        for (jj = 0; jj < nG; jj++)
                          {
                            Mr[ind][k] = E[ii][jj];        
                            ind = ind + 1; 
                          }                      
                    }

                  /* Integrate => E */

                  ind = 0; 
                  for (ii = 0; ii < nG; ii++)
                    for (jj = 0; jj < nG; jj++)
                      {
                        z = trapz(Nr, t, Mr[ind]); 
                        E[ii][jj].re = z.re*nrm2; 
                        E[ii][jj].im = z.im*nrm2; 
                        ind = ind + 1; 
                      }

                  /* Store to matrix M:  */

                  for (ii = 0; ii < nG; ii++)
                    for (jj = 0; jj < nG; jj++)
                      {
                        /* Add to existing data */

                        M[nS*nG + i*nG + ii][j*nG + jj] = 
                          c_add(M[nS*nG + i*nG + ii][j*nG + jj], E[ii][jj]); 

                        /* Check values */

                        CheckValue(FUNCTION_NAME, "M[ii][jj].re","", 
                                   M[nS*nG + i*nG + ii][j*nG + jj].re, 
                                   -INFTY, INFTY);

                        CheckValue(FUNCTION_NAME, "M[ii][jj].im","", 
                                   M[nS*nG + i*nG + ii][j*nG + jj].im, 
                                   -INFTY, INFTY);              
                      }
                }
            }
        }
    }

  /***************************************************************************/

  /* Solve c from Mc = b */
  
  LUdecomposition(nG*nJ, M, b, c);

  /* Free temporary variables */

  Mem(MEM_FREE, x); 
  Mem(MEM_FREE, x1); 
  Mem(MEM_FREE, x2); 
  Mem(MEM_FREE, y); 
  Mem(MEM_FREE, y1); 
  Mem(MEM_FREE, y2); 
  Mem(MEM_FREE, t); 

  for (i=0; i<nG; i++)
    Mem(MEM_FREE, E[i]);
  Mem(MEM_FREE, E);

  for (i=0; i<nG*nG; i++)
    Mem(MEM_FREE, Mr[i]);
  Mem(MEM_FREE, Mr);  

  for (i=0; i<nG*nJ; i++)
    Mem(MEM_FREE, M[i]);
  Mem(MEM_FREE, M);  

  /***************************************************************************/
}

/*****************************************************************************/
