/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : movedt.c                                       */
/*                                                                           */
/* Created:       2012/10/05 (JLe)                                           */
/* Last modified: 2019/03/30 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Moves particle forward using delta-tracking                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MoveDT:"

/*****************************************************************************/

long MoveDT(long part, double majorant, double minxs, long *cell, double *xs0,
            double *x, double *y, double *z, double *xt, double *yt,
            double *zt, double *l0, double *u, double *v, double *w, double E,
            long id)
{
  long ptr, type, mat, mat0, bc;
  double totxs, l, wgt;

  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Check coordinates and direction cosines */

  CheckValue(FUNCTION_NAME, "x", "", *x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", *y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", *z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "u", "", *u, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "v", "", *v, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "w", "", *w, -1.0, 1.0);

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Add to DT fraction counter */

  ptr = (long)RDB[RES_DT_TRACK_FRAC];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(1.0, 1.0, ptr, id, 2 - type);

  /* Double the majorant cross section in order to increase number of */
  /* collisions for sensitivity calculations. TODO: Add this to the   */
  /* extra majorants instead? */

  if ((long)RDB[DATA_SENS_MODE] != SENS_MODE_NONE)
    majorant *= 2;

  /* Check cross section and sample path length */

  CheckValue(FUNCTION_NAME, "majorant", "", majorant, ZERO, INFTY);
  l = -log(RandF(id))/majorant;

  /* Particles that came from another domain are identified from */
  /* a mismatch in MPI id. Path length is set to zero to enforce */
  /* virtual collision sampling on-site. */
  
  if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
      ((long)RDB[part + PARTICLE_MPI_ID] != mpiid))
    {
      /* Set the new MPI id */
      
      WDB[part + PARTICLE_MPI_ID] = mpiid;
      
      /* Set the path lenght to zero */
      
      l = 0.0;
    }
  
  /* Move particle to collision site */

  *x = *x + *u*l;
  *y = *y + *v*l;
  *z = *z + *w*l;

  /* Set path length */

  *l0 = l;

  /* Reset previous material pointer and total xs */

  mat0 = -1;
  totxs = 0.0;

  /* Find location */

  *cell = WhereAmI(*x, *y, *z, *u, *v, *w, id);
  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

  /* Check if point is outside the geometry (domain decomposition can */
  /* set l = 0) */

  if (((long)RDB[*cell + CELL_TYPE] == CELL_TYPE_OUTSIDE) && (l > 0.0))
    {
      /* Check if track is stopped at outer boundary */

      if ((long)RDB[DATA_STOP_AT_BOUNDARY] == NO)
        {
          /* Set weight to test that albedo boundary conditions were not */
          /* applied */

          wgt = 1.0;

          /* Apply repeated boundary conditions */

          bc = BoundaryConditions(cell, x, y, z, u, v, w, xt, yt, zt,
                                  &wgt, id);

          /* Check that condition was applied */

          if (bc != YES)
            Die(FUNCTION_NAME, "Repeated bc not apllied");

          /* Check change in weight */

          if (wgt != 1.0)
            Die(FUNCTION_NAME, "Change in weight (albedo)");

          /* Find location */

          *cell = WhereAmI(*x, *y, *z, *u, *v, *w, id);
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);
        }
      else
        {
          /* Stop track at boundary */

          StopAtBoundary(cell, x, y, z, l0, *u, *v, *w, id);

          /* Set cross section */

          *xs0 = -1.0;

          /* Check distance */

          CheckValue(FUNCTION_NAME, "l0", "", *l0, ZERO, INFTY);

          /* Return surface crossing */

          return TRACK_END_SURF;
        }
    }

  /* Get material pointer */

  mat = (long)RDB[*cell + CELL_PTR_MAT];
  mat = MatPtr(mat, id);

  /* Check domain */

  if (CheckDDDomain(mat) == NO)
    {
      /* Put material pointer and tracking mode */
      
      WDB[part + PARTICLE_PTR_MAT] = (double)mat;
      WDB[part + PARTICLE_DD_TRACK_MODE] = TRACK_MODE_DT;
      
      /* Exit */
      
      return TRACK_END_DD;
    }

  /* Check pointer */

  if (mat != mat0)
    {
      /* Get total cross section */

      totxs = TotXS(mat, type, E, id);

      /* Remember material */

      mat0 = mat;
    }

  /* Check total xs */

  CheckValue(FUNCTION_NAME, "totxs", "", totxs, 0.0, INFTY);

  /* Double the total cross section in order to increase number of collisions */
  /* for sensitivity calculations */
  /* TODO: Add this to the extra majorants instead (somehow) */

  if ((long)RDB[DATA_SENS_MODE] != SENS_MODE_NONE)
    totxs *= 2;

  /*****************************/

  /* Compare to minimum value */

  if (totxs > minxs)
    {
      /* Sample between real and virtual collision */

      if (RandF(id) < totxs/majorant)
        {
          /* Set cross section */

          *xs0 = totxs;

          /* Add to success counter */

          ptr = (long)RDB[RES_DT_TRACK_EFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY,ptr);
          AddBuf1D(1.0, 1.0, ptr, id, 2 - type);

          /* Check distance */

          CheckValue(FUNCTION_NAME, "l0", "", *l0, 0.0, INFTY);

          /* Return collision */

          return TRACK_END_COLL;
        }
      else
        {
          /* Set cross section */

          *xs0 = -1.0;

          /* Add to failure counter */

          ptr = (long)RDB[RES_DT_TRACK_EFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY,ptr);
          AddBuf1D(1.0, 1.0, ptr, id, 4 - type);

          /* Check distance */

          CheckValue(FUNCTION_NAME, "l0", "", *l0, 0.0, INFTY);

          /* Return virtual */

          return TRACK_END_VIRT;
        }
    }
  else
    {
      /* Sample scoring */

      if (RandF(id) < minxs/majorant)
        {
          /* Set cross section */

          *xs0 = minxs;

          /* Sample between real and virtual collision */

          if (RandF(id) < totxs/minxs)
            {
              /* Add to success counter */

              ptr = (long)RDB[RES_DT_TRACK_EFF];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY,ptr);
              AddBuf1D(1.0, 1.0, ptr, id, 2 - type);

              /* Check distance */

              CheckValue(FUNCTION_NAME, "l0", "", *l0, 0.0, INFTY);

              /* Return collision */

              return TRACK_END_COLL;
            }
          else
            {
              /* Add to failure counter */

              ptr = (long)RDB[RES_DT_TRACK_EFF];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY,ptr);
              AddBuf1D(1.0, 1.0, ptr, id, 4 - type);

              /* Check distance */

              CheckValue(FUNCTION_NAME, "l0", "", *l0, 0.0, INFTY);

              /* Return virtual */

              return TRACK_END_VIRT;
            }
        }
      else
        {
          /* Set cross section */

          *xs0 = -1.0;

          /* Add to failure counter */

          ptr = (long)RDB[RES_DT_TRACK_EFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY,ptr);
          AddBuf1D(1.0, 1.0, ptr, id, 4 - type);

          /* Check distance */

          CheckValue(FUNCTION_NAME, "l0", "", *l0, 0.0, INFTY);

          /* Return virtual */

          return TRACK_END_VIRT;
        }
    }
}

/*****************************************************************************/
