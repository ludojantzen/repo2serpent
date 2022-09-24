/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : azirot.c                                       */
/*                                                                           */
/* Created:       2011/02/28 (JLe)                                           */
/* Last modified: 2016/02/23 (JLe)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Rotates direction cosines around a random azimuthal angle    */
/*                                                                           */
/* Comments: - From Serpent 1.1.16 (23.8.2011)                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "AziRot:"

/*****************************************************************************/

void AziRot(double mu, double *u, double *v, double *w, long id)
{
  double u0, v0, w0;
  double C1, C2, C3;
  double rnd1, rnd2;
  long i;

  /* Check normalization of incident direction cosines */
  
  C1 = (*u)*(*u) + (*v)*(*v) + (*w)*(*w);

  if (fabs(C1 - 1.0) > 1E-9)
    {
      /* Check error */

      if (fabs(C1 - 1.0) > 1E-4)
        Die(FUNCTION_NAME, "Error in incident direction cosines");

      /* Re-normalize */

      C1 = sqrt(C1);

      *u = *u/C1;
      *v = *v/C1;
      *w = *w/C1;
    }

  /* Sanity check for mu and direction vectors (for NAN's etc.) */

  CheckValue(FUNCTION_NAME, "mu", "", mu, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "u", "", *u, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "v", "", *v, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "w", "", *w, -1.01, 1.01);

  /* Truncate mu */

  if (mu > 1.0)
    mu = 1.0;
  else if (mu < -1.0)
    mu = -1.0;

  /* Remember incident direction cosines */

  u0 = *u;
  v0 = *v;
  w0 = *w;

  /* Re-sampling loop */

  for (i = 0; i < 5; i++)
    {
      /* Sample like in MCNP (MCNP4C manual p. 2-38). */

      if ((C1 = 1.0 - w0*w0) > 1E-9)
        {
          /* Select two random numbers using the rejection criterion. */
          
          do
            {
              rnd1 = 1.0 - 2.0*RandF(id);
              rnd2 = 1.0 - 2.0*RandF(id);
            }
          while ((C2 = rnd1*rnd1 + rnd2*rnd2) > 1.0);
          
          C3 = sqrt((1.0 - mu*mu)/(C1*C2));
          
          *u = u0*mu + C3*(rnd1*u0*w0 - rnd2*v0); 
          *v = v0*mu + C3*(rnd1*v0*w0 + rnd2*u0); 
          *w = w0*mu - rnd1*C1*C3;
        }
      else
        {
          /* Same as previous but swap v and w */
          
          C1 = 1.0 - v0*v0;
          
          /* Select two random numbers using the rejection criterion. */
          
          do
            {
              rnd1 = 1.0 - 2.0*RandF(id);
              rnd2 = 1.0 - 2.0*RandF(id);
            }
          while ((C2 = rnd1*rnd1 + rnd2*rnd2) > 1.0);
          
          C3 = sqrt((1.0 - mu*mu)/(C1*C2));
          
          *u = u0*mu + C3*(rnd1*u0*v0 - rnd2*w0); 
          *w = w0*mu + C3*(rnd1*w0*v0 + rnd2*u0); 
          *v = v0*mu - rnd1*C1*C3;
        }

      /* Check final direction cosines */
      
      if (fabs((*u)*(*u) + (*v)*(*v) + (*w)*(*w) - 1.0) < 1E-6)
        break;
    }
  
  /* Check some values */
  
  CheckValue(FUNCTION_NAME, "final direction cosines", "", 
             (*u)*(*u) + (*v)*(*v) + (*w)*(*w) - 1.0, -1E-4, 1E-4);
  
  CheckValue(FUNCTION_NAME, "rotated vector", "", 
             mu - (u0*(*u) + v0*(*v) + w0*(*w)), -1E-4, 1E-4); 
}

/*****************************************************************************/
