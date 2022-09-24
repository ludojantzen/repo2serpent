/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : superdet.c                                     */
/*                                                                           */
/* Created:       2011/09/19 (JLe)                                           */
/* Last modified: 2019/12/02 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Scores super-imposed detectors                               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SuperDet:"

/*****************************************************************************/

void SuperDet(long part, double x0, double y0, double z0, double u, double v,
              double w, double lmax, double E, double t, double wgt, long id)
{
  long det, loc0, idx, surf, type, ptr, ptr1;
  long np, in0, in1, norm, loc1, rbin, i;
  long msh, fail;
  double d, b, l, val, cross, x, y, z, f, spd, un, vn, wn;
  const double *params;

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Calculate speed */

  spd = Speed(type, E);
  CheckValue(FUNCTION_NAME, "spd", "", spd, ZERO, INFTY);

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
      /* Compare particle types and check super-imposed */

      if (((long)RDB[part + PARTICLE_TYPE] != (long)RDB[det + DET_PARTICLE]) ||
          ((loc0 = (long)RDB[det + DET_PTR_SBINS]) < VALID_PTR))
        {
          /* Next detector */

          det = NextItem(det);

          /* Cycle loop */

          continue;
        }

      /* Loop over surface bins */

      while (loc0 > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[loc0 + DET_SBIN_TYPE] == SUPERDET_TYPE_CURRENT)
            {
              /***************************************************************/

              /***** Current detector ****************************************/

              /* Reset index */

              idx = -1;

              /* Get normal */

              norm = (long)RDB[loc0 + DET_SBIN_SURF_NORM];

              /* Start at previous position */

              x = x0;
              y = y0;
              z = z0;

              /* Reset distance */

              l = 0.0;

              /* Get surface pointer */

              surf = (long)RDB[loc0 + DET_SBIN_PTR_SURF];
              CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

              /* Get surface type */

              type = (long)RDB[surf + SURFACE_TYPE];

              /* Get number of parameters */

              np = (long)RDB[surf + SURFACE_N_PARAMS];

              /* Pointer to parameter list */

              ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              params = &RDB[ptr];

              /* Get initial position */

              in0 = TestSurface(surf, x, y, z, NO, id);

              /* Reset result */

              val = 0.0;

              /* Loop over all surfaces in track */

              while (l < lmax)
                {
                  /* Get distance */

                  d = SurfaceDistance(surf, params, type, np, x, y, z,
                                      u, v, w, id);

                  /* Check infinity */

                  if (d == INFTY)
                    break;

                  /* Extrapolate */

                  d = d + EXTRAP_L;

                  /* Update distance */

                  l = l + d;

                  if (l > lmax)
                    break;

                  /* Update coordinates */

                  x = x + d*u;
                  y = y + d*v;
                  z = z + d*w;

                  /* Test position */

                  in1 = TestSurface(surf, x, y, z, NO, id);

                  /* Reset crossing */

                  cross = 0.0;

                  /* Add to number of crossings */

                  if (norm == 0)
                    {
                      /* Net current */

                      if ((in0 == NO) && (in1 == YES))
                        cross = 1.0;
                      else if ((in0 == YES) && (in1 == NO))
                        cross = -1.0;
                    }
                  else if (norm == 1)
                    {
                      /* Outward current */

                      if ((in0 == YES) && (in1 == NO))
                        cross = -1.0;
                    }
                  else if (norm == -1)
                    {
                      /* Inward current */

                      if ((in0 == NO) && (in1 == YES))
                        cross = 1.0;
                    }
                  else if (norm == -2)
                    {
                      /* Check that surface was crossed */

                      if (in0 != in1)
                        {
                          /* Flux, get surface normal */

                          SurfaceNormal(surf, x, y, z, &un, &vn, &wn, id);

                          /* Calculate scalar product */

                          cross = fabs(u*un + v*vn + w*wn);
                          CheckValue(FUNCTION_NAME, "cross", "", cross,
                                     0.0, 1.0);

                          /* Check for zero */

                          if (cross > 0.0)
                            cross = 1.0/cross;
                          else
                            cross = 0.0;
                        }
                    }
                  else
                    Die(FUNCTION_NAME, "Invalid normal");

                  /* Check bin */

                  if ((i = DetBin(det, -1, part, x, y, z, E, t, id)) > -1)
                    {
                      /* Add to total */

                      val = val + cross;

                      /* Write to point to source file */

                      if (cross != 0.0)
                        WriteSourceFile(det, x, y, z, u, v, w, E, wgt, t,
                                        -1.0, id);

                      /* Set index and check */

                      if (idx < 0)
                        idx = i;
                      else if (idx != i)
                        Die(FUNCTION_NAME, "Change in bin index");
                    }

                  /* Put previous position */

                  in0 = in1;
                }

              /* Check value */

              if (val == 0.0)
                {
                  /* Pointer to next surface bin */

                  loc0 = NextItem(loc0);

                  /* Cycle loop */

                  continue;
                }

              /* Check index (should be set if value is non-zero) */

              if (idx < 0)
                Die(FUNCTION_NAME, "Index not set");

              /* Reset response index */

              rbin = 0;

              /* Get pointer to response functions */

              loc1 = (long)RDB[det + DET_PTR_RBINS];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

              /* Loop over responses */

              while (loc1 > VALID_PTR)
                {
                  /* Get response function */

                  f = DetResponse(det, loc1, part, -1, E, 1.0, id);

                  /* Get pointer to statistics */

                  ptr = (long)RDB[det + DET_PTR_STAT];
                  CheckPointer(FUNCTION_NAME, "(stat)", DATA_ARRAY, ptr);

                  /* Tää on uusien kirjastojen validointiin. Ensimmäiselle  */
                  /* ja viimeiselle detektorille binit menee kosinin mukaan */
                  /* energiajakaumaan. */

#ifdef mmmmmmmmmmmmmmmmmmmmmmmm

                  if ((NextItem(det) < VALID_PTR) ||
                      (PrevItem(det) < VALID_PTR))
                    {
                      idx = (long)(200.0*(w + 1.0)*0.5);
                      AddBuf(val, wgt, ptr, id, -1, idx, 0);
                    }
                  else

#endif
                    {
                      /* Add to statistics */

                      if ((long)RDB[det + DET_TYPE] == DETECTOR_TYPE_IMP_WGT)
                        AddBuf(val*f, 1.0, ptr, id, -1, idx, rbin);
                      else
                        AddBuf(val*f, wgt, ptr, id, -1, idx, rbin);
                    }

                  /* Contribution to response matrix */

                  ScoreRMXResp(part, det, val*f*wgt, id);

                  /* Score bilinear detector ratio */

                  if ((ptr1 = (long)RDB[det + DET_PTR_SENS_STAT_ARRAY])
                      > VALID_PTR)
                    {
                      /* Get fraction of events to skip */

                      f = RDB[det + DET_SKIP_FRAC];

                      /* Get correct stat from stat array */

                      ptr1 = (long)RDB[ptr1 + idx];
                      CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr1);

                      /* Score indirect and direct part */

                      if (RandF(id) > f)
                        EventsToSensitivity(part, wgt, ptr1, val, id);
                    }

                  /* Set flags */

                  SetDetFlags(det, part);

                  /* Update bin index */

                  rbin++;

                  /* Next response */

                  loc1 = NextItem(loc1);
                }

              /***************************************************************/
            }
          else if ((long)RDB[loc0 + DET_SBIN_TYPE] == SUPERDET_TYPE_TLEFLUX)
            {
              /***************************************************************/

              /***** Track-length detector within surface ********************/

              /* Start at previous position */

              x = x0;
              y = y0;
              z = z0;

              /* Reset distance */

              l = 0.0;

              /* Get surface pointer */

              surf = (long)RDB[loc0 + DET_SBIN_PTR_SURF];
              CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

              /* Get surface type */

              type = (long)RDB[surf + SURFACE_TYPE];

              /* Get number of parameters */

              np = (long)RDB[surf + SURFACE_N_PARAMS];

              /* Pointer to parameter list */

              ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              params = &RDB[ptr];

              /* Check initial position */

              in0 = TestSurface(surf, x, y, z, NO, id);

              /* Get pointer to mesh */

              msh = (long)RDB[det + DET_PTR_MESH];

              /* Loop over track */

              while (l < lmax)
                {
                  /* Get distance */


                  d = SurfaceDistance(surf, params, type, np, x, y, z,
                                      u, v, w, id);

                  /* Check if detector is associated with a mesh */

                  if (msh > VALID_PTR)
                    {
                      /* Distance to mesh boundaries */

                      b = NearestMeshBoundary(msh, x, y, z, u, v, w, &fail);

                      /* Compare */

                      if (b < d)
                        d = b;
                    }

                  /* Check infinity */

                  if (d == INFTY)
                    break;

                  /* When in0 == YES, the neutron is inside the surface */
                  /* and track-length estimator is scored */

                  if (in0 == YES)
                    {
                      /* Check that lmax is not exceeded */

                      if (l + d < lmax)
                        val = d + EXTRAP_L;
                      else
                        val = lmax - l;

                      /* Reset response index */

                      rbin = 0;

                      /* Get pointer to response functions */

                      loc1 = (long)RDB[det + DET_PTR_RBINS];
                      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

                      /* Loop over responses */

                      while (loc1 > VALID_PTR)
                        {
                          /* Get bin index */

                          if ((idx = DetBin(det, -1, part, x, y, z, E, t, id))
                              < 0)
                            break;

                          /* Skip recoil energy and source point types */

                          if (((long)RDB[loc1 + DET_RBIN_MT] !=
                               MT_SOURCE_RATE) &&
                              ((long)RDB[loc1 + DET_RBIN_MT] !=
                               MT_MACRO_RECOILE) &&
                              ((long)RDB[loc1 + DET_RBIN_MT] !=
                               MT_MACRO_HEATPHOTANA) &&
                              ((long)RDB[loc1 + DET_RBIN_MT] != MT_LEAK_RATE))
                            {
                              /* Get response function */

                              f = DetResponse(det, loc1, part, -1, E, 1.0, id);

                              /* Get pointer to statistics */

                              ptr = (long)RDB[det + DET_PTR_STAT];
                              CheckPointer(FUNCTION_NAME, "(stat)",
                                           DATA_ARRAY, ptr);

                              /* Add to statistics */

                              if ((long)RDB[det + DET_TYPE] ==
                                  DETECTOR_TYPE_IMP_WGT)
                                AddBuf(val*f, 1.0, ptr, id, -1, idx, rbin);
                              else
                                AddBuf(val*f, wgt, ptr, id, -1, idx, rbin);

                              /* Contribution to response matrix */

                              ScoreRMXResp(part, det, val*f*wgt, id);

                              /* Score bilinear detector ratio */

                              if ((ptr1 = (long)RDB[det + DET_PTR_SENS_STAT_ARRAY]) > VALID_PTR)
                                {
                                  /* Get fraction of events to skip */

                                  f = RDB[det + DET_SKIP_FRAC];

                                  /* Get correct stat from stat array */

                                  ptr1 = (long)RDB[ptr1 + idx];
                                  CheckPointer(FUNCTION_NAME, "(ptr1)",
                                               DATA_ARRAY, ptr1);

                                  /* Score indirect and direct part */

                                  if (RandF(id) > f)
                                    EventsToSensitivity(part, wgt, ptr1, val,
                                                        id);
                                }

                              /* Set flags */

                              SetDetFlags(det, part);
                            }

                          /* Update bin index */

                          rbin++;

                          /* Next response */

                          loc1 = NextItem(loc1);
                        }
                    }

                  /* Extrapolate */

                  d = d + EXTRAP_L;

                  /* Update distance */

                  l = l + d;

                  if (l > lmax)
                    break;

                  /* Update coordinates */

                  x = x + d*u;
                  y = y + d*v;
                  z = z + d*w;

                  /* Test position */

                  in0 = TestSurface(surf, x, y, z, NO, id);
                }

              /***************************************************************/
            }
          else
            Die(FUNCTION_NAME, "Invalid type");

          /* Next surface bin */

          loc0 = NextItem(loc0);
        }

      /* Next detector */

      det = NextItem(det);
    }
}

/*****************************************************************************/
