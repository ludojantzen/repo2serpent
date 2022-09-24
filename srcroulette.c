/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : srcroulette.c                                  */
/*                                                                           */
/* Created:       2017/09/08 (JLe)                                           */
/* Last modified: 2018/10/25 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Rejects source points when weight windows are in use         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SrcRoulette:"

/*****************************************************************************/

double SrcRoulette(long type, double x, double y, double z, double E, 
                   double wgt, double *p0, long id)
{
  double min, p;

  /* Get importance */
  
  p = WWImportance(type, x, y, z, 0.0, 0.0, 0.0, E, WWMESH_SRC);
  *p0 = p;

  /* Check */

  if (p == 0.0)
    return wgt;

  /* Lower weight window boundary */

  min = RDB[DATA_WWD_LOWER_BOUND]/p;

  /* Compare weight to boundary */
  
  if (wgt > min)
    return wgt;
  
  /* Calculate survival probability for russian roulette */
      
  p = wgt/min;

  /* Sample survival */
      
  if (RandF(id) < p)
    {
      /* Increase weight */
      
      wgt = wgt/p;
    }
  else
    {
      /* Kill particle */

      wgt = -1.0;
    }

  /* Return adjusted weight */

  return wgt;
}

/****************************************************************************/
