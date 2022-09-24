/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoreicmtrk.c                                  */
/*                                                                           */
/* Created:       2013/09/11 (JLe)                                           */
/* Last modified: 2019/09/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Scores ICM currents for coupling coefficients                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreICMTrk:"

/*****************************************************************************/

void ScoreICMTrk(long part, double x0, double y0, double z0, double u,
                 double v, double w, double lmax, double E, double wgt,
                 long id)
{
  long icm, ntot, nmu0, nmu1, nmu2, nseg, surf, type, ptr, np, in0, in1;
  long idx, ng, ncross, icm0, icm1, idx0, idx1, ng0, ng1;
  long mua, mus, mua0, mus0, mua1, mus1;
  long s0, s1, s2;
  double d, l, x, y, z, un0, vn0, wn0, un1, vn1, wn1, un2, vn2, wn2;
  double wgt0, wgt1, mu0, mu1, mu2;
  const double *params;

  /* Check mode */

  if ((long)RDB[DATA_ICM_CALC] == NO)
    return;

  /* Get sizes */

  ntot = (long)RDB[DATA_ICM_NG0];
  nmu0 = (long)RDB[DATA_ICM_NMU0];
  nmu1 = (long)RDB[DATA_ICM_NMU1];
  nmu2 = (long)RDB[DATA_ICM_NMU2];
  nseg = (long)RDB[DATA_ICM_NSEG];

  /* Get pointer to few-group structure */

  ptr = (long)RDB[DATA_ICM_PTR_ENE0];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get few-group index */

  if ((ng = GridSearch(ptr, E)) < 0)
    return;

  /* Convert index */

  ng = ntot - ng - 1;
  CheckValue(FUNCTION_NAME, "ng1", "", ng, 0, ntot - 1);

  /* Reset new data */

  icm1 = -1;
  idx1 = -1;
  mua1 = -1;
  mus1 = -1;
  ng1 = -1;
  wgt1 = -1;

  /* Loop over data */

  icm = (long)RDB[DATA_PTR_ICM0];
  while (icm > VALID_PTR)
    {
      /* Get original particle Albedo data */

      icm0 = (long)RDB[part + PARTICLE_ICM_PTR_ICM];
      idx0 = (long)RDB[part + PARTICLE_ICM_IDX];
      mua0 = (long)RDB[part + PARTICLE_ICM_MUA];
      mus0 = (long)RDB[part + PARTICLE_ICM_MUS];
      ng0 = (long)RDB[part + PARTICLE_ICM_G];
      wgt0 = RDB[part + PARTICLE_ICM_WGT];

      /* Get surface pointer */

      surf = (long)RDB[icm + ICM_PTR_SURF];
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

      /* Toi ekstrapolointi lienee turhaa sen jälkeen kun */
      /* trackin.c:tä muutettiin 14.8.2014 / 2.1.22 */

      /*Die(FUNCTION_NAME, "HUOM!!!");*/
      /*/
      x = x0 - 2.0*EXTRAP_L*u;
      y = y0 - 2.0*EXTRAP_L*v;
      z = z0 - 2.0*EXTRAP_L*w;
      */

      x = x0;
      y = y0;
      z = z0;

      /* Reset total distance and number of crossings */

      l = 0.0;
      ncross = 0;

      /* Get initial position */

      in0 = TestSurface(surf, x, y, z, NO, id);
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

              /* Get surface index */

              if ((idx = ICMIdx(surf, x, y, z, &un0, &vn0, &wn0, &un1, &vn1,
                                &wn1, &un2, &vn2, &wn2)) > -1)
                {
                  /* Calculate cosines */

                  mu0 = un0*u + vn0*v + wn0*w;
                  mu1 = un1*u + vn1*v + wn1*w;
                  mu2 = un2*u + vn2*v + wn2*w;

                  /* Get bins */

                  ptr = (long)RDB[DATA_ICM_PTR_MU0];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  s0 = SearchArray(&RDB[ptr], fabs(mu0), nmu0 + 1);
                  CheckValue(FUNCTION_NAME, "s0", "", s0, 0, nmu0 - 1);

                  ptr = (long)RDB[DATA_ICM_PTR_MU1];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  s1 = SearchArray(&RDB[ptr], mu1, nmu1 + 1);
                  CheckValue(FUNCTION_NAME, "s1", "", s1, 0, nmu1 - 1);

                  ptr = (long)RDB[DATA_ICM_PTR_MU2];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  s2 = SearchArray(&RDB[ptr], mu2, nmu2 + 1);
                  CheckValue(FUNCTION_NAME, "s2", "", s2, 0, nmu2 - 1);

                  /* Calculate asymmetric and symmetric part */

                  mua = s1;
                  mus = s0*nmu2 + s2;

                  /* Check direction */

                  if ((in0 == NO) && (in1 == YES))
                    {
                      /* Check bins */

                      CheckValue(FUNCTION_NAME, "idx", "(1)", idx, 0,
                                 nseg - 1);
                      CheckValue(FUNCTION_NAME, "mua", "(1)", mua, 0,
                                 nmu1 - 1);
                      CheckValue(FUNCTION_NAME, "mus", "(1)", mus, 0,
                                 nmu0*nmu2 - 1);
                      CheckValue(FUNCTION_NAME, "ng", "(1)", ng, 0,
                                 ntot - 1);

                      /* Score inward current */

                      if (RDB[DATA_CYCLE_IDX] > RDB[DATA_CRIT_SKIP] - 1.0)
                        {
                          ptr = (long)RDB[icm + ICM_RES_CURR0];
                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                          AddBuf(1.0, wgt, ptr, id, -1, idx, mua, mus, ng);
                        }

                      /* Put new data and reset weight */

                      icm0 = icm;
                      idx0 = idx;
                      mua0 = mua;
                      mus0 = mus;
                      ng0 = ng;
                      wgt0 = 1.0;
                    }
                  else if ((in0 == YES) && (in1 == NO))
                    {
                      /* Check active cycles */

                      if (RDB[DATA_CYCLE_IDX] > RDB[DATA_CRIT_SKIP] - 1.0)
                        {
                          /* Check index, group and angular bins */

                          if ((idx0 > -1) && (ng0 > -1) && (mua0 > -1) &&
                              (mus0 > -1))
                            {
                              /* Check bins */

                              CheckValue(FUNCTION_NAME, "idx0", "(2)", idx0, 0,
                                         nseg - 1);
                              CheckValue(FUNCTION_NAME, "mua0", "(2)", mua0, 0,
                                         nmu1 - 1);
                              CheckValue(FUNCTION_NAME, "mus0", "(2)", mus0, 0,
                                         nmu0*nmu2 - 1);
                              CheckValue(FUNCTION_NAME, "ng0", "(2)", ng0, 0,
                                         ntot - 1);
                              CheckValue(FUNCTION_NAME, "idx", "(2)", idx, 0,
                                         nseg - 1);
                              CheckValue(FUNCTION_NAME, "mua", "(2)", mua, 0,
                                         nmu1 - 1);
                              CheckValue(FUNCTION_NAME, "mus", "(2)", mus, 0,
                                         nmu0*nmu2 - 1);

                              /* Score outward current */

                              ptr = (long)RDB[icm + ICM_RES_CC1];
                              CheckPointer(FUNCTION_NAME, "(ptr)",
                                           DATA_ARRAY, ptr);
                              AddBuf(1.0, wgt, ptr, id, -1, idx0, mua0, mus0,
                                     ng0, idx, mua, mus, ng);

                              ptr = (long)RDB[icm + ICM_RES_CC2];
                              CheckPointer(FUNCTION_NAME, "(ptr)",
                                           DATA_ARRAY, ptr);
                              AddBuf(wgt0, wgt, ptr, id, -1, idx0, mua0, mus0,
                                     ng0, idx, mua, mus, ng);
                            }
                          else if ((ng0 > -1) && (RDB[DATA_CYCLE_IDX] >
                                                  RDB[DATA_CRIT_SKIP] + 20.0))
                            Warn(FUNCTION_NAME, "Problem in geometry?");
                        }

                      /* Reset values */

                      icm0 = -1;
                      idx0 = -1;
                      mua0 = -1;
                      mus0 = -1;
                      ng0 = -1;
                      wgt0 = 1.0;
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

          if (icm1 > VALID_PTR)
            Die(FUNCTION_NAME, "Point (%E, %E, %E) is in multiple regions",
                x, y, z);
          else
            {
              icm1 = icm0;
              idx1 = idx0;
              mua1 = mua0;
              mus1 = mus0;
              ng1 = ng0;
              wgt1 = wgt0;
            }

          /* Check number of crossings */

          if (ncross == 0)
            {
              /* Surface index cannot have changed */

              if (idx0 != (long)RDB[part + PARTICLE_ICM_IDX])
                Die(FUNCTION_NAME, "This is impossible");

              /* No other crossings are possible, add to counter */

              ptr = (long)RDB[icm + ICM_BREAK_PTR_COUNT];
              CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
              AddPrivateData(ptr, 1, id);

              /* Break loop */

              break;
            }
        }

      /* Next */

      icm = NextItem(icm);
    }

  /* Store values */

  WDB[part + PARTICLE_ICM_PTR_ICM] = (double)icm1;
  WDB[part + PARTICLE_ICM_IDX] = (double)idx1;
  WDB[part + PARTICLE_ICM_MUA] = (double)mua1;
  WDB[part + PARTICLE_ICM_MUS] = (double)mus1;
  WDB[part + PARTICLE_ICM_G] = (double)ng1;
  WDB[part + PARTICLE_ICM_WGT] = wgt1;
}

/*****************************************************************************/
