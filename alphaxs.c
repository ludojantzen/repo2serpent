/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : alphaxs.c                                      */
/*                                                                           */
/* Created:       2012/11/01 (JLe)                                           */
/* Last modified: 2012/11/01 (JLe)                                           */
/* Version:       2.1.10                                                     */
/*                                                                           */
/* Description: Returns cross section for alpha-eigenvalue mode              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AlphaXS:"

/*****************************************************************************/

double AlphaXS(double E)
{
  double alpha, xs, spd;

  /* Get speed */
  
  spd = Speed(PARTICLE_TYPE_NEUTRON, E);  

  /* Check eigenvalue */

  if ((alpha = RDB[DATA_ALPHA_EIG]) > 0.0)
    xs = alpha/spd;
  else if (alpha < 0.0)
    xs = -alpha/spd;
  else
    xs = 0.0;

  /* Return */

  return xs;
}

/*****************************************************************************/
