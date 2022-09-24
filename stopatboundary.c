/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stopatboundary.c                               */
/*                                                                           */
/* Created:       2012/08/11 (JLe)                                           */
/* Last modified: 2015/07/23 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Stops particle at outer boundary                             */
/*                                                                           */
/* Comments: - Needed to account for handling leakage and repeated boundary  */
/*             conditions in delta-tracking mode                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StopAtBoundary:"

/*****************************************************************************/

void StopAtBoundary (long *cell0, double *x0, double *y0, double *z0, 
                     double *l0, double u, double v, double w, long id)
{
  double x, y, z, l, d;
  const double *params;
  long cell, loc0, surf, ptr, n, type, np;

  /* Check cell pointer */

  CheckPointer(FUNCTION_NAME, "(cell0)", DATA_ARRAY, *cell0);

  /* Move particle back to previous position */

  x = *x0 - *l0*u;
  y = *y0 - *l0*v;
  z = *z0 - *l0*w;

  /* Reset traveled distance */

  *l0 = 0.0;

  /* Move forward until particle is outside (tän pitäis olla turha) */

  for (n = 0; n < 1000; n++)
    {
      /* Reset minimum distance */

      l = INFTY;

      /* Loop over boundaries */

      loc0 = (long)RDB[DATA_PTR_OUTER_BOUNDS];
      while (loc0 > VALID_PTR)
        {
          /* Pointer to surface */
          
          surf = (long)RDB[loc0 + BOUNDS_PTR_SURF];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Get surface type */
              
          type = (long)RDB[surf + SURFACE_TYPE];

          /* Get number of parameters */
          
          np = (long)RDB[surf + SURFACE_N_PARAMS];

          /* Pointer to parameter list */
              
          ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          params = &RDB[ptr];

          /* Calculate distance */
                  
          d = SurfaceDistance(surf, params, type, np, x, y, z, u, v, w, id);

          /* Compare to minimum */

          if (d > 0.0)
            if (d < l)
              l = d;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Check distance */

      if (l == INFTY)
        Die(FUNCTION_NAME, "Crossed boundary not within line-of-sight");

      /* Update position */
          
      x = x + (l + EXTRAP_L)*u;
      y = y + (l + EXTRAP_L)*v;
      z = z + (l + EXTRAP_L)*w;

      /* Update traveled distance */

      *l0 = *l0 + l + EXTRAP_L;

      /* Find cell */
  
      cell = WhereAmI(x, y, z, u, v, w, id);
      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

      /* Break loop if outside */
      
      if ((long)RDB[cell + CELL_TYPE] == CELL_TYPE_OUTSIDE)
        break;
    }

  /* Check for infinite loop */

  if (n == 1000)
    {
      /* Tracking error */
      
      TrackingError(TRACK_ERR_OUTSIDE, -1, -1, -1, id);
    }

  /* Update coordinates */

  *x0 = x;
  *y0 = y;
  *z0 = z;

  /* Update cell pointer */

  *cell0 = cell;
}

/*****************************************************************************/
