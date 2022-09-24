/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collectsensresults.c                           */
/*                                                                           */
/* Created:       2017/04/05 (VVa)                                           */
/* Last modified: 2018/10/19 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Collects sensitivity calculation results                     */
/*                                                                           */
/* Comments:   -Called from CollectResults()                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectSensResults:"

/*****************************************************************************/

void CollectSensResults()
{
  long g, i, maxi, nmat, nzai, nrea, nene, nmu;
  long det1, det2, loc0, detstat1, detstat2, type;
  long dir1ptr, dir2ptr, ptr, binidx;
  long izai, irea, iene, imat, mingen, stp, stp1, stp2, resp;
  double det1_val, val1;
  double det2_val, val2;
  double wgt1, wgt2, div1, div2;
  double P1, P2, plus;
  double hitmiss;

  /* No need to collect sensitivity results during inactive cycles */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Get pointer to the sensitivity block or exit the subroutine */

  if ((loc0 = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  /* Reduce scoring buffer */

  ReduceBuffer();

  /* Calculate maximum label index */

  nmat = (long)RDB[loc0 + SENS_N_MAT];
  nzai = (long)RDB[loc0 + SENS_N_ZAI];
  nrea = (long)RDB[loc0 + SENS_N_PERT] + 1;
  nene = (long)RDB[loc0 + SENS_N_ENE] + 1;
  nmu = (long)RDB[loc0 + SENS_N_MU];

  maxi = (1 + (nmat*nzai*nrea*nene*nmu));

  /* Get minimum latent generation */

  if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_SCORE_FLAG_HIS)
    mingen = 0;
  else
    mingen = (long)RDB[DATA_SENS_LAST_GEN] - 1;

  /* Loop over generations to collect results with various */
  /* numbers of latent generations */

  for (g = mingen; g < (long)RDB[DATA_SENS_LAST_GEN]; g++)
    {
      /* Loop over basic responses */

      resp = (long)RDB[loc0 + SENS_PTR_RESP0];

      while (resp > VALID_PTR)
        {
          /* Do not collect responses that have partials */

          if ((long)RDB[resp + SENS_RESP_HAS_PARTIALS] == YES)
            {
              /* Process next response */

              resp = NextItem(resp);
              continue;
            }

          /* Get response type */

          type = (long)RDB[resp + SENS_RESP_TYPE];

          /* Check response type */

          if (type != SENS_RESP_TYPE_RATIO)
            {
              /* Basic responses */

              /* Get basic wgt and divider wgt here */

              stp1 = (long)RDB[resp + SENS_RESP_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(stp1)", DATA_ARRAY, stp1);

              /* NB: void sensitivity uses "val" not "wgt" (avoid negative wgt) */

              if (type != SENS_RESP_TYPE_VOID)
                wgt1 = BufWgt(stp1, g, 0);
              else
                wgt1 = BufVal(stp1, g, 0);

              /* Divider */

              stp2 = (long)RDB[resp + SENS_RESP_PTR_STAT_DIVIDER];

              /* K-eff doesn't have a divider, set divider weight to zero here */

              if (stp2 > VALID_PTR)
                wgt2 = BufWgt(stp2, g, 0);
              else
                wgt2 = 0.0;

              div1 = wgt1;
              div2 = wgt2;

              det1_val = 0.0;
              det2_val = 0.0;

              dir1ptr = 0;
              dir2ptr = 0;
            }
          else
            {
              /* Detector responses */

              binidx = (long)RDB[resp + SENS_RESP_DET_BIN_IDX];

              /* Get pointer to det 1 */

              det1 = (long)RDB[resp + SENS_RESP_PTR_DET1];
              CheckPointer(FUNCTION_NAME, "(det1)", DATA_ARRAY, det1);

              /* Get pointer to statistics of detector 1 */

              detstat1 = (long)RDB[det1 + DET_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(detstat1)", DATA_ARRAY, detstat1);

              /* Get pointer to stat 1 array */

              stp1 = (long)RDB[det1 + DET_PTR_SENS_STAT_ARRAY];
              CheckPointer(FUNCTION_NAME, "(stp1)", DATA_ARRAY, stp1);

              /* Get pointer to stat1 */

              stp1 = (long)RDB[stp1 + binidx];
              CheckPointer(FUNCTION_NAME, "(stp1)", DATA_ARRAY, stp1);

              /* Get pointer to det 2 */

              det2 = (long)RDB[resp + SENS_RESP_PTR_DET2];
              CheckPointer(FUNCTION_NAME, "(det2)", DATA_ARRAY, det2);

              /* Get pointer to statistics of detector 2 */

              detstat2 = (long)RDB[det2 + DET_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(detstat2)", DATA_ARRAY, detstat2);

              /* Get pointer to stat 2 array */

              stp2 = (long)RDB[det2 + DET_PTR_SENS_STAT_ARRAY];
              CheckPointer(FUNCTION_NAME, "(stp2)", DATA_ARRAY, stp2);

              /* Get pointer to stat2 */

              stp2 = (long)RDB[stp2 + binidx];
              CheckPointer(FUNCTION_NAME, "(stp2)", DATA_ARRAY, stp2);

              /* Get detector values */

              det1_val = BufVal(detstat1, binidx, 0);
              det2_val = BufVal(detstat2, binidx, 0);

              /* Get probability to score each detector (needed for normalizing buffers) */

              if ((P1 = RDB[resp + SENS_RESP_DET1_P]) == 0.0)
                P1 = 1.0;

              if ((P2 = RDB[resp + SENS_RESP_DET2_P]) == 0.0)
                P2 = 1.0;

              /* Prepare divisors */

              div1 = det1_val*P1;
              div2 = det2_val*P2;

              /******************************************************/
              /*** The direct term(s) could be added here somehow ***/
              /******************************************************/

              /* We'll need (detector value due to imat, izai, irea, iene)/(detector value) */

              dir1ptr = 0;
              dir2ptr = 0;

              /* Get pointers to direct stats (might not exist) */

              if ((ptr = (long)RDB[det1 + DET_PTR_RBINS]) > VALID_PTR)
                {
                  /* Get pointer to array */

                  dir1ptr = (long)RDB[ptr + DET_RBIN_PTR_SENS_DIRECT_STAT_ARRAY];

                  if (dir1ptr > VALID_PTR)
                    {
                      /* If array exists, get pointer to correct bin */

                      dir1ptr = (long)RDB[dir1ptr + binidx];
                      CheckPointer(FUNCTION_NAME, "(dir1ptr)", DATA_ARRAY, dir1ptr);
                    }
                }

              /* Get pointers to direct stats (might not exist) */

              if ((ptr = (long)RDB[det2 + DET_PTR_RBINS]) > VALID_PTR)
                {
                  /* Get pointer to array */

                  dir2ptr = (long)RDB[ptr + DET_RBIN_PTR_SENS_DIRECT_STAT_ARRAY];

                  if (dir2ptr > VALID_PTR)
                    {
                      /* If array exists, get pointer to correct bin */

                      dir2ptr = (long)RDB[dir2ptr + binidx];
                      CheckPointer(FUNCTION_NAME, "(dir2ptr)", DATA_ARRAY, dir2ptr);
                    }
                }
            }

          /* Get pointer to result statistic here */

          stp = (long)RDB[resp + SENS_RESP_PTR_STAT];
          CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);

          /* Loop over perturbations (could probably be an OMP for if AddStat for */
          /* different bins doesn't interfere) */

          /* Start parallel timer */

          StartTimer(TIMER_OMP_PARA);

#ifdef OPEN_MP
#pragma omp parallel private(i, imat, izai, irea, iene, hitmiss, val1, val2, plus)
#endif
          {

#ifdef OPEN_MP
#pragma omp for schedule(dynamic)
#endif
            for (i = 1; i < maxi; i++)
              {
                /*************************************************************/

                /* Extract data from label */

                ExtractSensLabel(i, &imat, &izai, &irea, &iene, &hitmiss);

                /* Get first value */

                if (div1 != 0)
                  val1 = SensBufVal(stp1, g, imat, izai, irea, iene)/div1;
                else
                  val1 = 0.0;

                /* Get second value */

                if (div2 != 0)
                  val2 = SensBufVal(stp2, g, imat, izai, irea, iene)/div2;
                else
                  val2 = 0.0;

                /* Modify val1 and val2 with the direct term */

                if (dir1ptr > VALID_PTR)
                  {
                    if (det1_val != 0)
                      plus = SensBufVal(dir1ptr, -1, imat, izai, irea, iene)/det1_val;
                    else
                      plus = 0.0;

                    val1 = val1 + plus;
                  }

                if (dir2ptr > VALID_PTR)
                  {
                    if (det2_val != 0)
                      plus = SensBufVal(dir2ptr, -1, imat, izai, irea, iene)/det2_val;
                    else
                      plus = 0.0;

                    val2 = val2 + plus;
                  }

                /* Store value to statistics */

                AddStat(val1 - val2, stp, g, i);
              }
          }

          /* Stop parallel timer */

          StopTimer(TIMER_OMP_PARA);

          /* Next response */

          resp = NextItem(resp);
        }
    }
}

/*****************************************************************************/
