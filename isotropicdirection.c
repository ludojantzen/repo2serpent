/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : isotropicdirection.c                           */
/*                                                                           */
/* Created:       2010/11/12 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Samples isotropic direction cosines                          */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
  
#define FUNCTION_NAME "IsotropicDirection:"

/*****************************************************************************/

void IsotropicDirection(double *u, double *v, double *w, long id)
{
  double rand1, rand2, rand3, C1, C2;

  /* Use the rejection method described in Lux & Koblinger, pp. 21-22. */

  do
    {
      rand1 = 2.0*RandF(id) - 1.0;
      rand2 = 2.0*RandF(id) - 1.0;
      C1 = rand1*rand1 + rand2*rand2;
    }
  while (C1 > 1.0);
  
  rand3 = 2.0*RandF(id) - 1.0;
  
  C2 = sqrt(1 - rand3*rand3);
  
  *u = C2*(rand1*rand1 - rand2*rand2)/C1;
  *v = C2*2.0*rand1*rand2/C1;
  *w = rand3;

  /* Check value and exit. */

  CheckValue(FUNCTION_NAME, "r", "", *u**u+*v**v+*w**w - 1.0, -1E-10, 1E-10);
}

/*****************************************************************************/
