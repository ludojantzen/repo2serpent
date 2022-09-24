/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collectdet.c                                   */
/*                                                                           */
/* Created:       2011/03/03 (JLe)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Collects detector results                                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectDet:"

/*****************************************************************************/

void CollectDet()
{
  long det0, ptr, erg, type, stp, ebins0, ubins0, cbins0, mbins0, lbins0;
  long rbins0, zbins0, ybins0, xbins0, tbins0, eb0, ub0, cb0, mb0, lb0;
  long rb0, zb0, yb0, xb0, tb0, idx0, idx2, rb2, tb, tme;
  double val, norm, div, mul, cum;
  double tmin, tmax, tming, tmaxg;

  /* Reduce scoring buffer */

  ReduceBuffer();

  /* Avoid compiler warning */

  idx2 = -1;
  rb2 = -1;

  /* Get transport time interval */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
    {
      /* Get time bin index */

      tb = (long)RDB[DATA_DYN_TB];

      /* Get pointer to time bin structure */

      tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);

      /* Get pointer to bins */

      tme = (long)RDB[tme + TME_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(tme ptr)", DATA_ARRAY, tme);

      /* Get transport time interval */

      tming = RDB[tme + tb];
      tmaxg = RDB[tme + tb + 1];
    }
  else
    {
      /* Criticality source mode  */
      /* Set interval to infinity */

      tming = -INFTY;
      tmaxg = INFTY;
    }

  /* Loop over detectors */

  det0 = (long)RDB[DATA_PTR_DET0];
  while(det0 > VALID_PTR)
    {
      /* Get detector type */

      type = (long)RDB[det0 + DET_TYPE];

      /* Reset cumulative */

      cum = 0.0;

      /* Get normalization factor */

      norm = NormCoef((long)RDB[det0 + DET_PARTICLE]);
      CheckValue(FUNCTION_NAME, "norm", "", norm, 0.0, INFTY);

      /* Get number of bins */

      ebins0 = (long)RDB[det0 + DET_N_EBINS];
      ubins0 = (long)RDB[det0 + DET_N_UBINS];
      cbins0 = (long)RDB[det0 + DET_N_CBINS];
      mbins0 = (long)RDB[det0 + DET_N_MBINS];
      lbins0 = (long)RDB[det0 + DET_N_LBINS];
      rbins0 = (long)RDB[det0 + DET_N_RBINS];
      tbins0 = (long)RDB[det0 + DET_N_TBINS];

      /* Mesh bins */

      if ((ptr = (long)RDB[det0 + DET_PTR_MESH]) > VALID_PTR)
        {
          xbins0 = (long)RDB[ptr + MESH_N0];
          ybins0 = (long)RDB[ptr + MESH_N1];
          zbins0 = (long)RDB[ptr + MESH_N2];
        }
      else
        {
          xbins0 = 1;
          ybins0 = 1;
          zbins0 = 1;
        }

      /* Pointer to statistics */

      stp = (long)RDB[det0 + DET_PTR_STAT];
      CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);

      /* Pointer to energy distribution */

      if ((erg = (long)RDB[det0 + DET_PTR_EGRID]) > VALID_PTR)
        {
          /* Get pointer to data */

          erg = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
          CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
        }

      /* Get pointer to time binning (might not exist) */

      tme = (long)RDB[det0 + DET_PTR_TME];

      /* Check whether detector is using FETs */

      if ((ptr = (long)RDB[det0 + DET_FET_PTR_PARAMS]) > VALID_PTR)
        CollectFET(&RDB[ptr], stp, 0);
      else
        {
          /* Loop over bins */

          for (eb0 = 0; eb0 < ebins0; eb0++)
          for (ub0 = 0; ub0 < ubins0; ub0++)
          for (cb0 = 0; cb0 < cbins0; cb0++)
          for (mb0 = 0; mb0 < mbins0; mb0++)
          for (lb0 = 0; lb0 < lbins0; lb0++)
          for (rb0 = 0; rb0 < rbins0; rb0++)
          for (zb0 = 0; zb0 < zbins0; zb0++)
          for (yb0 = 0; yb0 < ybins0; yb0++)
          for (xb0 = 0; xb0 < xbins0; xb0++)
          for (tb0 = 0; tb0 < tbins0; tb0++)
            {
              /* Get index */

              idx0 = DetIdx(det0, eb0, ub0, cb0, mb0, lb0, zb0, yb0, xb0, tb0);

              /* Get detector time bin limits */

              if (tme > VALID_PTR)
                {
                  tmin = RDB[tme + tb0];
                  tmax = RDB[tme + tb0 + 1];
                }
              else
                {
                  /* Criticality source mode or no time bin limits */

                  tmin = -INFTY;
                  tmax = INFTY;
                }

              /* In dynamic (coupled) simulation mode, do not update */
              /* detector bins that are not in the current time bin. */
              /* Bins that are only partly in this time interval are */
              /* problematic. */

              if (((tmin < tming - 1E-15) || (tmax > tmaxg + 1E-15)) &&
                  ((tmax - tmin) < INFTY))
                {

                  /* Check if detector time binning does not match */
                  /* transport time binning */

                  if (((tmin < tming-1E-15) && (tmax > tming+1E-15)) ||
                      ((tmin < tmaxg-1E-15) && (tmax > tmaxg+1E-15)))
                    Warn(FUNCTION_NAME,
                         "Detector time bin (%E %E) split by "
                         "transport time bin (%E %E) (%E %E)",
                         tmin, tmax, tming, tmaxg);

                  /* Cycle time bin loop */

                  continue;
                }

              /* Divide by volume */

              div = RDB[det0 + DET_VOL];

              /* Energy and lethargy divider */

              if (type == DETECTOR_TYPE_UNI_E)
                {
                  /* Check pointer */

                  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

                  /* Divide by energy interval */

                  div = div*(RDB[erg + eb0 + 1] - RDB[erg + eb0]);
                }
              else if (type == DETECTOR_TYPE_UNI_L)
                {
                  /* Check pointer */

                  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

                  /* Divide by lethargy interval */

                  div = div*(log(RDB[erg + eb0 + 1]) - log(RDB[erg + eb0]));
                }

              /* Time divider */
              /*
                ptr = (long)RDB[det0 + DET_PTR_TME];
                if (ptr > VALID_PTR)
                div = div*(RDB[ptr + tb0 + 1] - RDB[ptr + tb0]);
              */
              /* Reset multiplier */

              mul = 1.0;

              /* Detector divider and multiplier */

              if ((type == DETECTOR_TYPE_MULTI) || (type == DETECTOR_TYPE_DIVI))
                {
                  /* Pointer to second detector */

                  ptr = (long)RDB[det0 + DET_PTR_MUL];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Check number of values */

                  if ((long)RDB[ptr + DET_N_TOT_BINS] == 1)
                    {
                      /* Single-valued divisor */

                      idx2 = 0;
                      rb2 = 0;
                    }
                  else if ((long)RDB[ptr + DET_N_TOT_BINS] ==
                           (long)RDB[det0 + DET_N_TOT_BINS])
                    {
                      /* Equal number of bins */

                      idx2 = idx0;
                      rb2 = rb0;
                    }
                  else
                    Die(FUNCTION_NAME, "Mismatch in number of bins");

                  /* Pointer to statistics */

                  ptr = (long)RDB[ptr + DET_PTR_STAT];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Get multiplier or divider */

                  if (type == DETECTOR_TYPE_MULTI)
                    mul = norm*BufVal(ptr, idx2, rb2);
                  else if (type == DETECTOR_TYPE_DIVI)
                    div = norm*div*BufVal(ptr, idx2, rb2);
                  else
                    Die(FUNCTION_NAME, "Error in type");
                }

              /* Get buffer value */

              if (type == DETECTOR_TYPE_SUM_SCORES)
                val = BufN(stp, idx0, rb0);
              else
                val = BufVal(stp, idx0, rb0);

              /* Normalize */

              val = val*norm;

              /* Check divisor */

              if (div != 0.0)
                {
                  /* Check cumulative */

                  if (type == DETECTOR_TYPE_CUMU)
                    {
                      cum = cum + mul*val/div;
                      AddStat(cum, stp, idx0, rb0);
                    }
                  else
                    AddStat(mul*val/div, stp, idx0, rb0);
                }
            }
        }

      /* Next detector */

      det0 = NextItem(det0);
    }
}

/*****************************************************************************/
