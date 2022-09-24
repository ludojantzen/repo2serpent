/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoredf.c                                      */
/*                                                                           */
/* Created:       2013/03/05 (JLe)                                           */
/* Last modified: 2018/12/20 (ARi)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores fluxes and currends used for calculating              */
/*              discontinuity factors                                        */
/*                                                                           */
/* Comments: - Noita n2, m2 ja l2 -muuttujia ei käytetä                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreDF:"

/*****************************************************************************/

void ScoreDF(double x0, double y0, double z0, double u, double v, double w,
             double lmax, double E, double wgt, long id)
{
  long adf, gcu, ntot, ng, surf, type, ptr, np, in0, in1, n1, n2, m1, m2;
  long l1, l2;
  double d, l, x, y, z, mu, u0, v0, w0, sgn;
  const double *params;

  /* Check that group constants are calculated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Check if active cycle */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Check if discontinuity factors are calculated and get pointer */

  if ((adf = (long)RDB[DATA_PTR_ADF0]) < VALID_PTR)
    return;

  /* Check for corrector calculation */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    return;

  /* Get pointer to microgroup energy grid */

  ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Number of groups */

  ntot = (long)RDB[ptr + ENERGY_GRID_NE] - 1;

  /* Get group index */

  if ((ng = GridSearch(ptr, E)) < 0)
    return;
  else
    ng = ntot - ng - 1;

  /* Loop over discontinuity factors */

  while (adf > VALID_PTR)
    {
      /* Pointer to gcu universe */

      gcu = (long)RDB[adf + ADF_PTR_GCU];
      CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);

      /* Get surface pointer */

      surf = (long)RDB[adf + ADF_PTR_SURF];
      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

      /* Get surface type */

      type = (long)RDB[surf + SURFACE_TYPE];

      /* Get number of parameters */

      np = (long)RDB[surf + SURFACE_N_PARAMS];

      /* Pointer to parameter list */

      ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      params = &RDB[ptr];

      /* Move to starting position */

      x = x0;
      y = y0;
      z = z0;

      /* Reset total distance */

      l = 0.0;

      /* Get initial position */

      in0 = TestSurface(surf, x, y, z, NO, id);

      /* Loop over all surfaces in track */

      while (l < lmax)
        {
          /* Get distance */

          d = SurfaceDistance(surf, params, type, np, x, y, z, u, v, w, id);

          /* Extrapolate */

          d = d + EXTRAP_L;

          /* Update coordinates */

          x = x + d*u;
          y = y + d*v;
          z = z + d*w;

          /* Test position */

          in1 = TestSurface(surf, x, y, z, NO, id);

          /* Check with maximum */

          if (l + d > lmax)
            {
              /* Score TLE of cell flux */

              if (in0 == YES)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_ADF_CELL_FLUX];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  AddPrivateRes(ptr + ng, (lmax - l)*wgt, id);
                }

              /* Break loop */

              break;
            }
          else
            {
              /* Score TLE of cell flux */

              if (in0 == YES)
                {
                  ptr = (long)RDB[gcu + GCU_MICRO_ADF_CELL_FLUX];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  AddPrivateRes(ptr + ng, d*wgt, id);
                }

              /* Update distance */

              l = l + d;
            }

          /* Check if surface was crossed */

          if (in0 != in1)
            {
              /* Get surface, corner, normal and sign */

              DFPos(surf, x, y, z, &n1, &n2, &m1, &m2, &l1, &l2,
                    &u0, &v0, &w0, &sgn);

              /* Check values */

              CheckValue(FUNCTION_NAME, "n1", "", n1, -1, 5);
              CheckValue(FUNCTION_NAME, "m1", "", m1, -1, 5);
              CheckValue(FUNCTION_NAME, "l1", "", l1, -1, 5);

              /* Calculate cosine */

              if ((mu = fabs(u*u0 + v*v0 + w*w0)) < 0.1)
                mu = 0.05;

              /* Check surface index */

              if (n1 > -1)
                {
                  /* Score surface flux */

                  ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_FLUX];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  AddPrivateRes(ptr + ng + n1*ntot, wgt/mu, id);

                  /* Sign moments of discontinuity factors */

                  ptr = (long)RDB[gcu + GCU_MICRO_ADF_SGN_SURF_FLUX];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  AddPrivateRes(ptr + ng + n1*ntot, sgn*wgt/mu, id);

                  /* Score inward or outward surface average current */

                  if ((in0 == NO) && (in1 == YES))
                    {
                      ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_IN_CURR];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddPrivateRes(ptr + ng + n1*ntot, wgt, id);

                      /* Sign moments of discontinuity factors */

                      ptr = (long)RDB[gcu + GCU_MICRO_ADF_SGN_SURF_IN_CURR];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddPrivateRes(ptr + ng + n1*ntot, sgn*wgt, id);
                    }
                  else if ((in0 == YES) && (in1 == NO))
                    {
                      ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_OUT_CURR];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddPrivateRes(ptr + ng + n1*ntot, wgt, id);

                      /* Sign moments of discontinuity factors */

                      ptr = (long)RDB[gcu + GCU_MICRO_ADF_SGN_SURF_OUT_CURR];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddPrivateRes(ptr + ng + n1*ntot, sgn*wgt, id);
                    }
                }

              /* Score inward or outward surface mid-point current */

              if (l1 > -1)
                {
                  if ((in0 == NO) && (in1 == YES))
                    {
                      ptr = (long)RDB[gcu + GCU_MICRO_ADF_MID_IN_CURR];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddPrivateRes(ptr + ng + l1*ntot, wgt, id);
                    }
                  else if ((in0 == YES) && (in1 == NO))
                    {
                      ptr = (long)RDB[gcu + GCU_MICRO_ADF_MID_OUT_CURR];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddPrivateRes(ptr + ng + l1*ntot, wgt, id);
                    }
                }

              /* Check Corner index */

              if (m1 > -1)
                {
                  /* Score corner flux */

                  ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_FLUX];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  AddPrivateRes(ptr + ng + m1*ntot, wgt/mu, id);

                  /* Score inward or outward corner current */

                  if ((in0 == NO) && (in1 == YES))
                    {
                      ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_IN_CURR];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddPrivateRes(ptr + ng + m1*ntot, wgt, id);
                    }
                  else if ((in0 == YES) && (in1 == NO))
                    {

                      ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_OUT_CURR];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      AddPrivateRes(ptr + ng + m1*ntot, wgt, id);
                    }
                }
            }

          /* Put previous position */

          in0 = in1;
        }

      /* Next adf */

      adf = NextItem(adf);
    }
}

/*****************************************************************************/
