/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stopatwwbound.c                                */
/*                                                                           */
/* Created:       2015/10/02 (JLe)                                           */
/* Last modified: 2019/02/21 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Stops track at weight window boundary                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StopAtWWBound:"

/*****************************************************************************/

long StopAtWWBound(long type, long trk, double *x, double *y, double *z,
                   double x0, double y0, double z0, double u, double v,
                   double w, double spd, double *t, double *lmax, long *cell,
                   long id)
{
  double d;

  /* Check for weight windows and response matrix calculation */

  if (((long)RDB[DATA_USE_WEIGHT_WINDOWS] == NO) &&
      ((long)RDB[DATA_RMTX_CALC] == NO))
    return trk;

  /* Distance to boundary */

  if ((d = WWDis(type, x0, y0, z0, u, v, w)) > *lmax)
    return trk;

  /* Check */

  CheckValue(FUNCTION_NAME, "d", "", d, ZERO, INFTY);

  /* Extrapolate */

  d = d + EXTRAP_L;

  /* Adjust time and distence */

  *t = *t - (*lmax - d)/spd;
  *lmax = d;

  /* Move over boundary */

  *x = x0 + d*u;
  *y = y0 + d*v;
  *z = z0 + d*w;

  /* Get new location */

  *cell = WhereAmI(*x, *y, *z, u, v, w, id);
  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

  /* Check type (added 7.7.2016 / 2.1.27) */

  if ((long)RDB[*cell + CELL_TYPE] == CELL_TYPE_OUTSIDE)
    return TRACK_END_SURF;
  else
    return TRACK_END_WWIN;
}

/*****************************************************************************/
