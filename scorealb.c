/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorealb.c                                     */
/*                                                                           */
/* Created:       2014/06/25 (JLe)                                           */
/* Last modified: 2018/12/19 (ARi)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores currents used for calculating albedos and partial     */
/*              albedos                                                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreAlb:"

/*****************************************************************************/

void ScoreAlb(long part, double x0, double y0, double z0, double u, double v,
              double w, double lmax, double E, double wgt, long id)
{
  long alb, gcu, ntot, ng, surf, type, ptr, np, in0, in1, n1, n2, m1, m2;
  long l1, l2, gcu0, ng0, idx0, gcu1, ng1, idx1, ncross, ns, flip;
  double d, l, x, y, z, u0, v0, w0, dum;
  const double *params;

  /* Check that group constants are calculated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Check if albedos are calculated and get pointer */

  if ((alb = (long)RDB[DATA_PTR_ALB0]) < VALID_PTR)
    return;

  /* Number of energy groups */

  ntot = (long)RDB[DATA_ERG_FG_NG];

  /* Get pointer to few-group structure */

  ptr = (long)RDB[DATA_ERG_FG_PTR_GRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get few-group index */

  if ((ng = GridSearch(ptr, E)) > -1)
    {
      /* Convert index */

      ng = ntot - ng - 1;
      CheckValue(FUNCTION_NAME, "ng", "", ng, 0, ntot - 1);
    }
  else
    return;

  /* Reset new data */

  gcu1 = -1;
  idx1 = -1;
  ng1 = -1;

  /* Loop over albedos */

  while (alb > VALID_PTR)
    {
      /* Get original particle Albedo data */

      gcu0 = (long)RDB[part + PARTICLE_ALB_PTR_GCU];
      idx0 = (long)RDB[part + PARTICLE_ALB_SURF_IDX];
      ng0 = (long)RDB[part + PARTICLE_ALB_G];

      /* Get number of surfaces */

      ns = (long)RDB[alb + ALB_NSURF];
      CheckValue(FUNCTION_NAME, "ns", "", ns, 1, 8);

      /* Pointer to gcu universe */

      gcu = (long)RDB[alb + ALB_PTR_GCU];
      CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);

      /* Get surface pointer */

      surf = (long)RDB[alb + ALB_PTR_SURF];
      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

      /* Get direction flip flag (TODO: muuta toi nimi) */

      flip = (long)RDB[alb + ALB_DIR];
      CheckValue(FUNCTION_NAME, "flip", "", flip, 0, 1);

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

      /* Reset total distance and number of crossings */

      l = 0.0;
      ncross = 0;

      /* Get initial position */

      in0 = TestSurface(surf, x, y, z, NO, id);
      in0 = in0^flip;
      in1 = -1;

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
          in1 = in1^flip;

          /* Check with maximum */

          if (l + d > lmax)
            {
              /* Go back to collision site (for debugging) */

              x = x0 + lmax*u;
              y = y0 + lmax*v;
              z = z0 + lmax*w;

              /* Cancel crossing */

              in1 = in0;

              /* Break loop */

              break;
            }
          else
            {
              /* Update distance */

              l = l + d;
            }

          /* Check if surface was crossed */

          if (in0 != in1)
            {
              /* Add counter */

              ncross++;

              /* Get surface, corner and normal (most values not used */
              /* here, because the subroutine was written for ADF's). */

              DFPos(surf, x, y, z, &n1, &n2, &m1, &m2, &l1, &l2,
                    &u0, &v0, &w0, &dum);

              /* Check values */

              CheckValue(FUNCTION_NAME, "n1", "", n1, -1, 5);
              CheckValue(FUNCTION_NAME, "n2", "", n2, -1, 5);
              CheckValue(FUNCTION_NAME, "m1", "", m1, -1, 5);
              CheckValue(FUNCTION_NAME, "m2", "", m2, -1, 5);
              CheckValue(FUNCTION_NAME, "l1", "", l1, -1, 5);
              CheckValue(FUNCTION_NAME, "l2", "", l2, -1, 5);

              /* Check surface index */

              if (n1 > -1)
                {
                  /* Check direction */

                  if ((in0 == NO) && (in1 == YES))
                    {
                      /* Score inward current */

                      if (RDB[DATA_CYCLE_IDX] > RDB[DATA_CRIT_SKIP] - 1.0)
                        {
                          ptr = (long)RDB[gcu + GCU_MICRO_ALB_IN_CURR];
                          CheckPointer(FUNCTION_NAME, "(ptr)",
                                       RES2_ARRAY, ptr);
                          AddPrivateRes(ptr + n1*ntot + ng, wgt, id);
                        }

                      /* Put new data */

                      gcu0 = gcu;
                      idx0 = n1;
                      ng0 = ng;
                    }
                  else if ((in0 == YES) && (in1 == NO))
                    {
                      /* Score outward current */

                      if (RDB[DATA_CYCLE_IDX] > RDB[DATA_CRIT_SKIP] - 1.0)
                        {
                          /* Check index and group */

                          if ((idx0 > -1) && (ng0 > -1))
                            {
                              ptr = (long)RDB[gcu + GCU_MICRO_ALB_OUT_CURR];
                              CheckPointer(FUNCTION_NAME, "(ptr)",
                                           RES2_ARRAY, ptr);
                              AddPrivateRes(ptr + idx0*ntot*ntot*ns +
                                            ng0*ntot*ns + n1*ntot + ng,
                                            wgt, id);
                            }
                          else if ((ng0 > -1) && (RDB[DATA_CYCLE_IDX] >
                                                  RDB[DATA_CRIT_SKIP] + 20.0))
                            Warn(FUNCTION_NAME, "Problem in geometry?");
                        }

                      /* Reset values */

                      gcu0 = -1;
                      idx0 = -1;
                      ng0 = -1;
                    }
                }
            }

          /* Put previous position */

          in0 = in1;
        }

      /* Check if particle is in */

      if (in1 == YES)
        {
          /* Update new values */

          if (gcu1 > VALID_PTR)
            Die(FUNCTION_NAME, "Point (%E, %E, %E) is in multiple regions",
                x, y, z);
          else
            {
              gcu1 = gcu0;
              idx1 = idx0;
              ng1 = ng0;
            }

          /* Check number of crossings */

          if (ncross == 0)
            {
              /* Surface index cannot have changed */

              if (idx0 != (long)RDB[part + PARTICLE_ALB_SURF_IDX])
                Die(FUNCTION_NAME, "This is impossible");

              /* Break loop */

              break;
            }
        }

      /* Next albedo */

      alb = NextItem(alb);
    }

  /* Store values */

  WDB[part + PARTICLE_ALB_PTR_GCU] = (double)gcu1;
  WDB[part + PARTICLE_ALB_SURF_IDX] = (double)idx1;
  WDB[part + PARTICLE_ALB_G] = (double)ng1;
}

/*****************************************************************************/
