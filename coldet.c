/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : coldet.c                                       */
/*                                                                           */
/* Created:       2011/03/03 (JLe)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores collision flux detectors                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ColDet:"

/*****************************************************************************/

void ColDet(long part, long mat0, double flx0, double x0, double y0,
            double z0, double u0, double v0, double w0, double E0, double t0,
            double wgt0, double g0, long id)
{
  long det, loc0, idx, rbin, ptr, ptr1, type, mt, loc1;
  double f0, val0, u, v, w, f;

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
      /*********************************************/

      /* Tää skippaa ensimmäisen detekorin kytketyn n/p-moodin */
      /* testausvaiheessa */

#ifdef WANHAAA

      if (((long)RDB[DATA_PHOTON_PRODUCTION] != NO) &&
          (PrevItem(det) < VALID_PTR))
        {
          /* Next detector */

          det = NextItem(det);

          /* Cycle loop */

          continue;
        }

#endif

      /*********************************************/

      /* Compare particle types and check super-imposed */

      if ((type != (long)RDB[det + DET_PARTICLE]) ||
          ((long)RDB[det + DET_PTR_SBINS] > VALID_PTR))
        {
          /* Next detector */

          det = NextItem(det);

          /* Cycle loop */

          continue;
        }

      /* Get bin index */

      if (1 == 2)
        {
          /* Tää on collision binnauksen testausta varten (JLe 28.7.2015) */

          if ((idx = (long)RDB[part + PARTICLE_COL_IDX]) >
              (long)RDB[det + DET_N_TOT_BINS] - 1)
            return;

          flx0 = 1.0;
        }
      else
        {
          /* Get bin */

          if ((idx = DetBin(det, mat0, part, x0, y0, z0, E0, t0, id)) < 0)
            {
              /* Next detector */

              det = NextItem(det);

              /* Cycle loop */

              continue;
            }
        }

      /* Reset response index */

      rbin = 0;

      /* Get pointer to response functions */

      loc0 = (long)RDB[det + DET_PTR_RBINS];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Loop over responses */

      while (loc0 > VALID_PTR)
        {
          /* Get mt */

          mt = (long)RDB[loc0 + DET_RBIN_MT];

          /* Check special cases */

          if ((mt == MT_MACRO_RECOILE) || (mt == MT_SOURCE_RATE) ||
              (mt == MT_MACRO_HEATPHOTANA) || (mt == MT_PHOTON_PULSE_HEIGHT) ||
              (mt == MT_LEAK_RATE) ||
              ((mt >= MT_SECONDARY_DN_SOURCE_G8) && (mt <= MT_PRIMARY_LIVE_SOURCE)) ||
              ((mt >= MT_ELECTRON_AUGER) && (mt <= MT_ELECTRON_PE)))
            {
              /* Update response index */

              rbin++;

              /* Next response */

              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }

          /* Get response */

          if (mt == MT_MACRO_HEATTOT)
            {
              /* Get fission response */

              loc1 = NextItem(loc0);
              if (loc1 < VALID_PTR)
                Die(FUNCTION_NAME, "Should not go here.");
              f0 = DetResponse(det, loc1, part, mat0, E0, g0, id);

              /* Get kerma response */

              if ((long)RDB[DATA_EDEP_MODE] > EDEP_MODE_MT458)
                {
                  loc1 = NextItem(loc1);
                  if (loc1 < VALID_PTR)
                    Die(FUNCTION_NAME, "Should not go here.");
                  f0 = f0 + DetResponse(det, loc1, part, mat0, E0, g0, id);
                }
            }
          else
            f0 = DetResponse(det, loc0, part, mat0, E0, g0, id);

          /* Calculate value */

          if ((val0 = f0*flx0) == 0.0)
            {
              /* Update response index */

              rbin++;

              /* Next response */

              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }

          /* Get direction vector */

          u = RDB[det + DET_DIRVEC_U];
          v = RDB[det + DET_DIRVEC_V];
          w = RDB[det + DET_DIRVEC_W];

          /* Check if non-zero */

          if ((u != 0.0) || (v != 0.0) || (w != 0.0))
            {
              /* Calculate scalar product */

              val0 = val0*(u*u0 + v*v0 + w*w0);

              /* Check negative */

              if (val0 < 0.0)
                val0 = 0.0;
            }

          /*******************************************************************/

          /***** Store result ************************************************/

          /* Get pointer to statistics */

          ptr = (long)RDB[det + DET_PTR_STAT];
          CheckPointer(FUNCTION_NAME, "(stat)", DATA_ARRAY, ptr);

          /* Check if FET type */

          if ((ptr1 = (long)RDB[det + DET_FET_PTR_PARAMS]) > VALID_PTR)
            ScoreFET(&RDB[ptr1], ptr, part, val0, wgt0, x0, y0, z0, rbin, id);

          /* Check index */

          CheckValue(FUNCTION_NAME, "idx", "", idx, 0, 1000000000);

          /* Score */

          if ((long)RDB[det + DET_TYPE] == DETECTOR_TYPE_IMP_WGT)
            AddBuf(val0, 1.0, ptr, id, -1, idx, rbin);
          else
            AddBuf(val0, wgt0, ptr, id, -1, idx, rbin);

          /* Write to point to source file */

          WriteSourceFile(det, x0, y0, z0, u0, v0, w0, E0, wgt0, t0, flx0, id);

          /* Score bilinear detector ratio */

          if ((ptr1 = (long)RDB[det + DET_PTR_SENS_STAT_ARRAY]) > VALID_PTR)
            {
              /* Get fraction of events to skip */

              f = RDB[det + DET_SKIP_FRAC];

              /* Get correct stat from stat array */

              ptr1 = (long)RDB[ptr1 + idx];
              CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr1);

              /* Score indirect and direct part */

              if (RandF(id) > f)
                {
                  EventsToSensitivity(part, wgt0, ptr1, val0, id);
                  ScoreDirectSensitivity(part, wgt0, idx, loc0, mat0,
                                         E0, flx0, g0, id);
                }
            }

          /* Contribution to response matrix */

          ScoreRMXResp(part, det, val0*wgt0, id);

          /* Set flags */

          SetDetFlags(det, part);

          /*******************************************************************/

          /* Update response index */

          rbin++;

          /* Next response */

          loc0 = NextItem(loc0);
        }

      /* Next detector */

      det = NextItem(det);
    }
}

/*****************************************************************************/
