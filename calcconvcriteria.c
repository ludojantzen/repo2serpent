/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calcconvcriteria.c                             */
/*                                                                           */
/* Created:       2016/05/04 (VVa)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Calculates convergence criteria for a distribution           */
/*                                                                           */
/*                                                                           */
/* Comments: -Used in relaxinterfacepower and updateinterface                */
/*           -                                                               */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalcConvCriteria:"

/*****************************************************************************/

void CalcConvCriteria(double new, double old, double *maxdiff, double *maxeps,
                      double *L2abs, double *L2rel)
{
  double eps, diff;

  /* Calculate pointwise relative convergence criterion */
  /* NB: If new == 0, we'll say the relative difference is 1.0 */
  /* which is not really true*/

  if (new != 0)
    eps = fabs(new - old)/fabs(new);
  else
    eps = 1.0;

  if (eps > *maxeps)
    *maxeps = eps;

  /* Calculate pointwise absolute convergence criterion */

  diff = fabs(new - old);

  if (diff > *maxdiff)
    *maxdiff = diff;

  /* Add to absolute L2 norm of difference */

  *L2abs = *L2abs + (new - old)*(new - old);

  /* Add to relative L2 norm of difference */
  /* NB: If new == 0, we'll say the relative difference is 1.0 */
  /* which is not really true*/

  if (new != 0)
    *L2rel = *L2rel + (new - old)*(new - old)/(new*new);
  else
    *L2rel = *L2rel + 1.0;

  return;

  /***************************************************************************/
}

/*****************************************************************************/
