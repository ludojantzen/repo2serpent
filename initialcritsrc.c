/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : initialcritsrc.c                               */
/*                                                                           */
/* Created:       2017/10/05 (JLe)                                           */
/* Last modified: 2020/05/25 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Iterates an initial guess for criticality source             */
/*              distribution by applying the response-matrix method.         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitialCritSrc:"

/*****************************************************************************/

void InitialCritSrc()
{
  long no, ni, nn, idx, id, rr;
  long no_particles_left, maybe_done, done;

  /* Check if response matrix calculation is on */

  if ((long)RDB[DATA_RMTX_CALC] == NO)
    Die(FUNCTION_NAME, "No response matrix calculation");

  fprintf(outp, "Iterating initial fission source:\n\n");

  /* Set parameters needed to run the simulation */

  WDB[DATA_SIMULATION_MODE] = (double)SIMULATION_MODE_SRC;
  WDB[DATA_SRC_POP] = RDB[DATA_CRIT_POP];
  WDB[DATA_NPHYS_SAMPLE_FISS] = (double)NO;
  WDB[DATA_CYCLE_IDX] = RDB[DATA_CRIT_SKIP] + 1;
  WDB[DATA_RMTX_CALC] = YES;

  /* Set implicit reaction rates on */

  rr = (long)RDB[DATA_OPTI_IMPLICIT_RR];
  WDB[DATA_OPTI_IMPLICIT_RR] = (double)YES;

  /* Outer iteration loop */

  for (no = 0; no < (long)RDB[DATA_RMX_CONVG_OUT_N]; no++)
    {
      fprintf(outp,
              "Running source batches for iteration loop %2ld / %2ld...\n",
              no + 1, (long)RDB[DATA_RMX_CONVG_OUT_N]);

      /* Start timer */

      ResetTimer(TIMER_MISC);
      StartTimer(TIMER_MISC);

      /* Put index */

      WDB[DATA_RMX_CONVG_ITER_IDX] = (double)no;

      /* Loop over source batches */

      for (ni = 0; ni < (long)RDB[DATA_RMX_CONVG_SB]; ni++)
        {
          /* MGa: Reset the DD communications  */

          ResetDDComms();

          /* Parallel loop over histories */

#ifdef OPEN_MP
#pragma omp parallel private(id, idx, nn)
#endif
          {
            /* Get Open MP thread id */

            id = OMP_THREAD_NUM;

#ifdef OPEN_MP
#pragma omp for schedule(dynamic)
#endif
            /* Loop over source neutrons */

            for (nn = 0; nn < (long)RDB[DATA_CRIT_POP]; nn++)
              {
                /* Calculate particle index */

                idx = (long)RDB[DATA_NHIST_TOT];
                idx = idx + (long)RDB[DATA_CRIT_POP]*ni + nn;

                /* Sample source points */

                if (SampleSrcPoint(id, nn, idx) > VALID_PTR)
                  Tracking(id);
              }
          }

          /* MGa: Reset the flags used to manage termination in DD mode */

          no_particles_left = 0;
          maybe_done = 0;

          if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
            done = 0;
          else
            done = 1;

          /* MGa: This loop continues until all particles in the whole */
          /* system have been tracked */

          while (!done)
            {
              /* MGa: Send all the particles going to other domains */

              DistributeDDLimbo();

              /* MGa: Recieve all the particles coming from other domains */
              /* Note: if no particles are received termination is checked */

              ReceiveDDParticles(&no_particles_left);

              /* MGa: Check for termination if no particles have arrived */
              /* Note: first the termination is estimated using asynchronous */
              /* communications, so not an exact global check, then the */
              /* synchronous check is used if the tracking might be finished */

              if (no_particles_left)
                {
                  /* Check using asynchronous MPI functions */

                  CheckFinishDDAsynch(&maybe_done);

                  /* Check using a synchronous global reduction */

                  if (maybe_done)
                    CheckFinishDDSynch(&done);

                  /* Break or check for incoming particles again */

                  continue;
                }

              /* JLe: all particles should now be in the que  */

              /* Loop until source is empty */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
              {

                /* Get Open MP thread id */

                id = OMP_THREAD_NUM;

                /* Track particles in the que */

                Tracking(id);
              }
            }
        }

      /* Stop timer */

      StopTimer(TIMER_MISC);

      /* Print */

      fprintf(outp, "Completed in %s.\n\n",
              TimeIntervalStr(TimerVal(TIMER_MISC)));

      /* Update number of neutron histories run */

      WDB[DATA_NHIST_TOT] = RDB[DATA_NHIST_TOT] +
        RDB[DATA_RMX_CONVG_SB]*RDB[DATA_CRIT_POP];

      /* Collect data from MPI tasks */

      CollectParallelData();

      fprintf(outp, "Running response matrix solver...\n");

      /* Start timer */

      ResetTimer(TIMER_MISC);
      StartTimer(TIMER_MISC);

      /* Call response-matrix solver */

      ImportanceSolver();

      /* Stop timer */

      StopTimer(TIMER_MISC);

      /* Print */

      fprintf(outp, "Completed in %s.\n\n",
              TimeIntervalStr(TimerVal(TIMER_MISC)));

      /* Clear statistics */

      ClearStat(-1);
    }

  /* Put index */

  WDB[DATA_RMX_CONVG_ITER_IDX] = (double)no;

  /* Reset parameters to original values */

  WDB[DATA_SIMULATION_MODE] = (double)SIMULATION_MODE_CRIT;
  WDB[DATA_SRC_POP] = 0.0;
  WDB[DATA_NPHYS_SAMPLE_FISS] = (double)YES;
  WDB[DATA_CYCLE_IDX] = 0.0;
  WDB[DATA_RMTX_CALC] = NO;
  WDB[DATA_OPTI_IMPLICIT_RR] = (double)rr;
}

/*****************************************************************************/
