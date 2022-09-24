/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collectuncresults.c                            */
/*                                                                           */
/* Created:       2018/06/19 (VVa)                                           */
/* Last modified: 2019/01/18 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Collects uncertainty propagation results                     */
/*                                                                           */
/* Comments:   -Actually also does the propagation if needed                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectUncResults:"

/*****************************************************************************/

void CollectUncResults()
{
  long ptr;
  long gen, i, j;
  long mingen, resp, type, binidx;
  long loc1, sens, block, stp, stv, stp1, stp2, det1, det2, detstat1, detstat2;
  long izai, irea, iene, imat, idx, N, ng, dir1ptr, dir2ptr, blockidx;
  double val1, val2, wgt1, wgt2, div1, div2, P1, P2, det1_val, det2_val;
  double val, totval;
  double hitmiss, *Sarr;

  /* No need to collect uncertainty results during inactive cycles */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Get pointer to the sensitivity block or exit the subroutine */

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  /* Get pointer to the first covariance block or exit the subroutine */

  if ((block = (long)RDB[DATA_PTR_COVBLOCK0]) < VALID_PTR)
    return;

  /* Reduce scoring buffer */

  ReduceBuffer();

  /* We'll need to allocate memory for working arrays in this subroutine */

  Mem(MEM_ALLOW);

  /* Get minimum latent generation */

  if ((long)RDB[sens + SENS_RESP_FLAGS] & SENS_SCORE_FLAG_HIS)
    mingen = 0;
  else
    mingen = (long)RDB[DATA_SENS_LAST_GEN] - 1;

  /* Loop over generations to collect results with various */
  /* numbers of latent generations */

  for (gen = mingen; gen < (long)RDB[DATA_SENS_LAST_GEN]; gen++)
    {
      /* Loop over basic responses */

      resp = (long)RDB[sens + SENS_PTR_RESP0];

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
                wgt1 = BufWgt(stp1, gen, 0);
              else
                wgt1 = BufVal(stp1, gen, 0);

              /* Divider */

              stp2 = (long)RDB[resp + SENS_RESP_PTR_STAT_DIVIDER];

              /* K-eff doesn't have a divider, set divider weight to zero here */

              if (stp2 > VALID_PTR)
                wgt2 = BufWgt(stp2, gen, 0);
              else
                wgt2 = 0.0;

              div1 = wgt1;
              div2 = wgt2;

              det1_val = 1.0;
              det2_val = 1.0;

              dir1ptr = 0.0;
              dir2ptr = 0.0;
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
              /*
              fprintf(outp, "Response: %s\n", GetText(resp + SENS_RESP_PTR_NAME));

              fprintf(outp, "DET1: %s, DET2: %s\n", GetText((long)RDB[resp + SENS_RESP_PTR_DET1] + DET_PTR_NAME), GetText((long)RDB[resp + SENS_RESP_PTR_DET2] + DET_PTR_NAME));

              fprintf(outp, "divvals: %E %E\n", div1, div2);
              fprintf(outp, "stpptrs: %ld %ld\n", stp1, stp2);
              fprintf(outp, "dirptrs: %ld %ld\n", dir1ptr, dir2ptr);
              */
            }

          /* Get pointer to response uncertainty statistic */

          stp = (long)RDB[resp + SENS_RESP_PTR_UNC_STAT];
          stv = (long)RDB[resp + SENS_RESP_PTR_VAR_STAT];

          /* Loop over the connected covariance blocks */

          block = (long)RDB[DATA_PTR_COVBLOCK0];
          totval = 0.0;
          blockidx = 0;
          while (block > VALID_PTR)
            {
              /* Get block order */

              N = (long)RDB[block + COVBLOCK_ORDER];

              /* Get number of energy groups */

              ng = (long)RDB[block + COVBLOCK_NG];

              /* Allocate temporary array for the sensitivity array */

              Sarr = Mem(MEM_ALLOC, N*ng, sizeof(double));
              /*
              fprintf(outp, "Block contains the following covariances: ");
              ptr = (long)RDB[block + COVBLOCK_PTR_ZAIMT_ARRAY];
              for (i = 0; i< N; i++)
                fprintf(outp, "%ld ", (long)RDB[ptr+i]);
              fprintf(outp, "\n");
              */
              /**********************************/
              /* Build up the sensitivity array */
              /**********************************/

              loc1 = (long)RDB[block + COVBLOCK_PTR_SENS_INDICES];

              /* Start parallel timer */

              StartTimer(TIMER_OMP_PARA);

#ifdef OPEN_MP
#pragma omp parallel private(i, idx, imat, izai, irea, iene, hitmiss, val1, val2)
#endif
              {

#ifdef OPEN_MP
#pragma omp for schedule(dynamic)
#endif
                for (i = 0; i < N*ng; i++)
                  {
                    /* Get direct sensitivity label for i */

                    idx = (long)RDB[loc1 + i];

                    /* Extract data from label */

                    ExtractSensLabel(idx, &imat, &izai, &irea, &iene, &hitmiss);

                    /* Get first value */

                    if (div1 != 0)
                      val1 = SensBufVal(stp1, gen, imat, izai, irea, iene)/div1;
                    else
                      val1 = 0.0;

                    /* Get second value */

                    if (div2 != 0)
                      val2 = SensBufVal(stp2, gen, imat, izai, irea, iene)/div2;
                    else
                      val2 = 0.0;

                    /* Modify val1 and val2 with the direct term */

                    if ((dir1ptr > VALID_PTR) && (det1_val != 0.0))
                      val1 = val1 + SensBufVal(dir1ptr, -1, imat, izai, irea, iene)/det1_val;

                    if ((dir2ptr > VALID_PTR) && (det2_val != 0.0))
                      val2 = val2 + SensBufVal(dir2ptr, -1, imat, izai, irea, iene)/det2_val;

                    /* Store value to sensitivity array */

                    Sarr[i] = val1 - val2;
                  }
              }

              /* Stop parallel timer */

              StopTimer(TIMER_OMP_PARA);

              /* Get pointer to the covariance matrix of the block */

              ptr = RDB[block + COVBLOCK_PTR_MATRIX_DATA];

              /* Apply the sandwich rule S*COV*S^T */

              val = 0;

              for (i = 0; i < N*ng; i++)
                for (j = 0; j < N*ng; j++)
                  val = val + Sarr[i]*RDB[ptr + N*ng*i + j]*Sarr[j];

              /* The uncertainty may be negative due to non PSD covariance matrix */
              /* Now truncate it from below to zero here, in future preprocess */
              /* matrices to be PSD by removing negative eigenpairs */

              if (val < 0)
                val = 0.0;

              totval += val;

              AddStat(val, stv, gen, blockidx+1);

              /*
              fprintf(outp, "Block uncertainty value is %f percent\n", sqrt(val)*100.0);
              */
              /* Free the temporary array for the sensitivity array */

              Mem(MEM_FREE, Sarr);

              /* Store block uncertainty to statistic */
              val = sqrt(val);
              AddStat(val, stp, gen, blockidx+1);

              /* Next sub-block of the full matrix */

              block = NextItem(block);
              blockidx++;
            }

          AddStat(totval, stv, gen, 0);

          /* Need to take a square root of the variance */

          totval = sqrt(totval);
          /*
          fprintf(outp, "%.5f %% Total uncertainty (percent)\n", totval*100.0);
          */
          /* Store total uncertainty to statistic */

          AddStat(totval, stp, gen, 0);

          /* Get next response */

          resp = NextItem(resp);
        }
    }

  /* We'll need to deny memory allocation for future routines */

  Mem(MEM_DENY);

  return;

}

/*****************************************************************************/
