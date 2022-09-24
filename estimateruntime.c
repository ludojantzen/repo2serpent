/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : estimateruntime.c                              */
/*                                                                           */
/* Created:       2011/03/12 (JLe)                                           */
/* Last modified: 2019/03/13 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Estimates running time                                       */
/*                                                                           */
/* Comments: - Predictor-corrector -ym. askeleet pitää ottaa tohon mukaan    */
/*           - Coupled calculation iteraatioissa pitää arvioida jotenkin     */
/*             eri tavalla                                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "EstimateRuntime:"

/*****************************************************************************/

void EstimateRuntime()
{
  long i, tot, skip, np, nc, niter, ncoe, ncoet;
  long tb, nt, tme;
  double t0, t1, t2, t3, ct, tt;
  double tming, tmaxg, tmint, tmaxt;

  /* Get cycle index and number of skip cycles */

  i = (long)RDB[DATA_CYCLE_IDX] + 1;

  /* Set number of inactive batches */

  if(RDB[DATA_USE_FSP] == (double)NO)
    {
      /* No fission source passing*/

      skip = (long)RDB[DATA_CRIT_SKIP];
    }
  else if(RDB[DATA_BURN_STEP] + RDB[DATA_SOL_REL_ITER] == 0.0)
    {
      /* First transportcycle with fission source passing */

      skip = (long)RDB[DATA_CRIT_SKIP];
    }
  else
    {
      /* Subsequent transportcycle with fission source passing */

      skip = (long)RDB[DATA_FSP_CRIT_SKIP];
    }

  /* Check mode */

  if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC) ||
      ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DELDYN))
    {
      /* Previous cycle time */
      
      t1 = TimerVal(TIMER_TRANSPORT) - RDB[DATA_PREV_CYCLE_T0];
      WDB[DATA_PREV_CYCLE_T0] = TimerVal(TIMER_TRANSPORT);
      
      /* Use previous estimate if sort cycle */

      if (((long)RDB[DATA_CYCLE_IDX] - 1)*3 == (long)RDB[DATA_SORT_COUNT])
        ct = RDB[DATA_ESTIM_CYCLE_TIME];
      else
        {
          /* Time per cycle */
          
          if (i > 1)
            t1 = 0.5*RDB[DATA_PREV_CYCLE_T1] + 0.5*t1;
          
          WDB[DATA_PREV_CYCLE_T1] = t1;

          /* Get number of batches */
      
          tot = (long)RDB[DATA_SRC_BATCHES];
          
          /* Running time so far */
      
          t0 = TimerVal(TIMER_RUNTIME);

          /* Estimate total */
          
          ct = t0 + (double)(tot - i)*t1;
        }
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
    {

      /* Get current timestep */

      tb = (long)RDB[DATA_DYN_TB];

      /* Get total number of timesteps */

      nt = (long)RDB[DATA_DYN_NB];

      /* Get pointer to time bin structure */

      tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);

      /* Get pointer to bins */

      tme = (long)RDB[tme + TME_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(tme ptr)", DATA_ARRAY, tme);

      /* Get transport time interval */

      tming = RDB[tme + tb];

      tmaxg = RDB[tme + tb + 1];

      /* Get total time limits */

      tmint = RDB[tme + 0];

      tmaxt = RDB[tme + nt];

      /* Calculate simulated time before this interval */

      t0 = tming - tmint;

      /* Get number of batches */

      tot = (long)RDB[DATA_SRC_BATCHES];

      /* Get portion of simulated time on this interval */

      t0 = t0 + (tmaxg - tming)*((double)i/(double)tot);

      /* Get total time to be simulated */

      t1 = tmaxt - tmint;

      /* Init time */

      t2 = TimerVal(TIMER_RUNTIME) - TimerVal(TIMER_TRANSPORT);

      /* Total transport time (corresponds to t0) */

      t3 = TimerVal(TIMER_TRANSPORT);

      /* Proportion t0/t1 has been simulated */
      /* Estimate total time */

      ct = t2 + t3*t1/t0;

    }
  else if (i > skip)
    {
      /* Get number of criticality cycles */

#ifdef MPI_MODE2

      /* Check number of tasks */

      if (mpitasks > 1)
        {
          /* Calculate number of particles per task */

          tot = (long)(RDB[DATA_CRIT_CYCLES]/((double)mpitasks));

        }
      else
        tot = (long)RDB[DATA_CRIT_CYCLES];

#else

      tot = (long)RDB[DATA_CRIT_CYCLES];

#endif

      /* Init time */

      t0 = TimerVal(TIMER_RUNTIME) - TimerVal(TIMER_TRANSPORT_ACTIVE);

      /* Check gpop mode */

      if ((long)RDB[DATA_GROW_POP_SIM] == NO)
        {
          /* Time per cycle */
          
          t1 = TimerVal(TIMER_TRANSPORT_ACTIVE)/((double)(i - skip));
          
          /* Estimate total time */
          
          ct = t0 + (double)tot*t1;
        }
      else
        {
          /* Time per history */
          
          t1 = TimerVal(TIMER_TRANSPORT_ACTIVE)/RDB[DATA_NHIST_CYCLE];
          
          /* Estimate total time */
          
          ct = t0 + RDB[DATA_GROW_POP_MAX_HIS]*t1;
        }
    }
  else
    ct = 0.0;

  /* Put time */

  WDB[DATA_ESTIM_CYCLE_TIME] = ct;

  /* Number of predictor and corrector cycle left */

  np = (long)RDB[DATA_BURN_TOT_STEPS] - (long)RDB[DATA_BURN_PRED_STEP];

  if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_NONE)
    nc = 0;
  else
    {
      /* MODIFIED (AIs) */
      /* this only works if MAX iterations = MIN iterations */
      /* when convergence check is added this must change*/

      nc = (long)RDB[DATA_BURN_TOT_STEPS]*(long)RDB[DATA_BURN_CI_MAXI]
        - (long)RDB[DATA_BURN_CORR_STEP];
    }

  /* Multiply number of predictors and correctors with the number of */
  /* coupled calculation iterations */

  if ((long)RDB[DATA_RUN_CC] == YES)
    {
      /* Get number of iterations */

      if ((niter = (long)RDB[DATA_SOL_REL_MAX_ITER]) == 0)
        {
          /* Calculate number of iterations based on population */
          /* TODO: */

          niter = 1;
        }

      /* Modify number of predictors */

      np = np*niter;

      /* Modify number of correctors */

      nc = nc*niter;
    }

  /* Number of coefficient calculations left */

  ncoe = (long)RDB[DATA_TOT_COEF_CALC] - (long)RDB[DATA_COEF_CALC_IDX];
  ncoet = (long)RDB[DATA_COEF_CALC_TOT_RUNS] -
    (long)RDB[DATA_COEF_CALC_RUN_IDX];

  /* Estimate total running time with burnup and coefficient calculations */

  if (((long)RDB[DATA_BURN_STEP] > 0.0) ||
      ((long)RDB[DATA_COEF_CALC_IDX] > -1))
    {
      /* Add estimated time to complete this step */

      tt = ct;

      /* Add estimated processing time */

      tt = tt + (double)(np + nc)*TimerVal(TIMER_PROCESS);

      /* Add estimated burnup time */

      tt = tt + (double)(np + nc)*TimerVal(TIMER_BURNUP);

      /* Add estimated transport time */

      tt = tt + (double)np*RDB[DATA_PRED_TRANSPORT_TIME];
      tt = tt + (double)nc*RDB[DATA_CORR_TRANSPORT_TIME];

      /* Check coefficient calculation */

      if (((long)RDB[DATA_COEF_CALC_SPECIAL_MODE] !=
           SPECIAL_COEF_MODE_HIS_ONLY) &&
          ((long)RDB[DATA_PTR_COEF0] > VALID_PTR))
        {
          /* Add estimated transport time */

          tt = tt + (double)ncoet*RDB[DATA_COEF_TRANSPORT_TIME];

          /* Add estimated processing time */

          tt = tt + (double)ncoet*TimerVal(TIMER_PROCESS);

          /* Add estimated init time */

          tt = tt + (double)ncoe*TimerVal(TIMER_INIT);
        }
    }
 else
   tt = 0.0;

  /* Put time */

  WDB[DATA_ESTIM_TOT_TIME] = tt;
}

/*****************************************************************************/
