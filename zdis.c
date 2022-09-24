/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : zdis.c                                         */
/*                                                                           */
/* Created:       2015/10/07 (JLe)                                           */
/* Last modified: 2015/10/07 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Calculates distence to axial plane                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ZDis:"

/*****************************************************************************/

double ZDis(double z, double w, double d)
{  
  if (w == 0.0)
    return INFTY;
  else if ((d = -(z - d)/w) < 0.0)
    return INFTY;
  else
    return d;
}

/*****************************************************************************/
