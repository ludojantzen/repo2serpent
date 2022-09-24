/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : randf.c                                        */
/*                                                                           */
/* Created:       2010/11/22 (JLe)                                           */
/* Last modified: 2017/03/24 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Stops the track at surfaces where detector flags are set     */
/*                                                                           */
/* Comments: - This is needed to avoid setting and checking flags for the    */
/*             same track                                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StopSurfDetFlg:"

/*****************************************************************************/

long StopSurfDetFlg(long trk, double *x, double *y, double *z, double x0,
                    double y0, double z0, double u, double v, double w, 
                    double spd, double *t, double *lmax, long *cell, long id)
{
  long det, surf, np, ptr, type;
  double d, min;

  /* Reset minimum distance */

  min = INFTY;

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
      /* Check flagging and pointer to surface bins */

      if (((long)RDB[det + DET_PTR_FLAGGING] > VALID_PTR) && 
          ((surf = (long)RDB[det + DET_PTR_SBINS]) > VALID_PTR))
        {
          /* Get surface pointer */

          surf = (long)RDB[surf + DET_SBIN_PTR_SURF];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);
          
          /* Get type */

          type = (long)RDB[surf + SURFACE_TYPE];
          
          /* Get number of parameters */

          np = (long)RDB[surf + SURFACE_N_PARAMS];
          
          /* Pointer to parameter list */
          
          ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          /* Get distance */
          
          d = SurfaceDistance(-1, &RDB[ptr], type, np, x0, y0, z0,
                              u, v, w, id);
          CheckValue(FUNCTION_NAME, "d", "1", d, 0.0, INFTY);

          /* Compare to minimum */

          if ((d > 0.0) && (d < min))
            min = d;
        }

      det = NextItem(det);
    }

  /* Check distance */

  if (min > *lmax)
    return trk;

  /* Extrapolate */

  min = min + EXTRAP_L;

  /* Adjust time and distence */

  *t = *t - (*lmax - min)/spd;
  *lmax = min;

  /* Move over boundary */

  *x = x0 + min*u;
  *y = y0 + min*v;
  *z = z0 + min*w;

  /* Get new location */
              
  *cell = WhereAmI(*x, *y, *z, u, v, w, id);
  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

  /* Check type */

  if ((long)RDB[*cell + CELL_TYPE] == CELL_TYPE_OUTSIDE)
    return TRACK_END_SURF;
  else
    return TRACK_END_FLAG;
}

/*****************************************************************************/
