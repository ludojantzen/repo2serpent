/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : burnupcycle.c                                  */
/*                                                                           */
/* Created:       2011/05/23 (JLe)                                           */
/* Last modified: 2020/05/28 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Loop over burnup history                                     */
/*                                                                           */
/* Comments: the "step" variable indexes steps within current burnup interval*/
/*           while DATA_BURN_STEP  indexes total steps.                      */
/*                                                                           */
/*           See:                                                            */
/*                                                                           */
/*           A.E Isotalo & P.A. Aarnio "Higher Order methods for burnup      */
/*           calculations with Bateman solutions" Ann.Nucl.Energy 38 (2011)  */
/*           pp. 1987-1995.                                                  */
/*                                                                           */
/*           and                                                             */
/*                                                                           */
/*           A.E. Isotalo & P.A. Aarnio "Substep methods for burnup          */
/*           calculations with Bateman solutions" Ann.Nucl.Energy            */
/*           (Submitted)                                                     */
/*                                                                           */
/* for description of the burnup calculation methods used.                   */
/*                                                                           */
/* TODO? tätä vois ehkä  muuttaa niin että asteluvut, aliaskeleet ja tuleeko */
/*       correctori määritellään askeleen aluksi erillisessä funktiossa      */
/*       niin ettei tarkastuksia olisi siellä täällä                         */
/*                                                                           */
/*      -Voisi vaihtaa niin, että viimeinen correctori lasketaan             */
/*       paremmilla statistiikoilla, ehkä                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "BurnupCycle:"

/*****************************************************************************/

void BurnupCycle()
{
  long dep, type, step, steps, pc;
  char tmpstr[MAX_STR];

  /* Remove binary restart file */

  if (((long)RDB[DATA_WRITE_RESTART_FILE] == YES) &&
      ((long)RDB[DATA_RESTART_READ_CONTINUE] == NO))
    {
      /* Get file name */

      if ((long)RDB[DATA_RESTART_WRITE_PTR_FNAME] > VALID_PTR)
        sprintf(tmpstr, "%s", GetText(DATA_RESTART_WRITE_PTR_FNAME));
      else
        sprintf(tmpstr, "%s.wrk", GetText(DATA_PTR_INPUT_FNAME));

      /* Remove */

      remove(tmpstr);
    }

  /* First loop is over intervals */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Store pointer */

      WDB[DATA_BURN_PTR_CURRENT_DEP] = (double)dep;

      /* Set normalization */

      SetNormalization(dep);

      /* Get step type and number of steps */

      type = (long)RDB[dep + DEP_HIS_STEP_TYPE];
      steps = (long)RDB[dep + DEP_HIS_N_STEPS];

      /* Put type */

      WDB[DATA_BURN_STEP_TYPE] = (double)type;

      /* Add final step for last interval */

      if (NextItem(dep) < VALID_PTR)
        steps++;

      /* Reprocess */

      Reprocess(dep);

      /* Second loop is over steps */

      for (step = 0; step < steps; step++)
        {
          /* Check for history variation break point */

          if (TestHisvBreak() == YES)
            {
              /* Break loop */

              break;
            }

          /* Store step */

          WDB[DATA_BURN_PTR_CURRENT_STEP] = (double)step;

          /* Initialize corrector iteration index */

          WDB[DATA_BURN_CI_I] = (double)0;

          /* Reset corrector iteration stopping flag */

          if (RDB[DATA_BURN_SIE] == (double)YES)
            WDB[DATA_BURN_CI_LAST] = (double)NO;
          else
            WDB[DATA_BURN_CI_LAST] = (double)YES;

          /* Predictor-corrector -loop */

          for (pc = 0; pc < 2; pc++)
            {
             /* Set the predictor/corrector status. pc is not used directly */

              if (pc == 0)
                WDB[DATA_BURN_STEP_PC] = PREDICTOR_STEP;
              else
                {
                  WDB[DATA_BURN_STEP_PC] = CORRECTOR_STEP;

                  /* Signal the external program about moving */
                  /* to next burnup point */

                  if(RDB[DATA_BURN_CI_I] == 0)
                    SignalExternal(SIGUSR2);

                  /* Read data interfaces for end-of-step */

                  ReadDataInterfaces();
                }

              /* Alternate neutron population params for */
              /* Corrector iteration */

              SetCIPopulation();

              /* If we are not running SIE we'll run the transport cycle in */
              /* a normal fashion. If we are running SIE, we'll only run it */
              /* if we're on the first predictor or a corrector             */

              if ((RDB[DATA_BURN_SIE] == (double)NO) ||
                 (!((RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP) &&
                    (RDB[DATA_BURN_STEP] != 0.0))))
                {
                  /* Prepare transport cycle */

                  PrepareTransportCycle();

                  /* Transport calculation cycle if not decay step */

                  if ((type == DEP_STEP_DEC_STEP) ||
                      (type == DEP_STEP_DEC_TOT)  ||
                      (type == DEP_STEP_ACT_STEP) ||
                      (type == DEP_STEP_ACT_TOT))
                    {
                      /* Decay or activation step, set completed flag */

                      WDB[DATA_SIMULATION_COMPLETED] = (double)YES;

                      /* Calculate activities */

                      CalculateActivities();

                      /* Make exponential fit for decay heat */

                      ExpoDecFit();

                      /* Print output */

                      MatlabOutput();
                    }
                  else
                    {
                      /* Run transportcycle(s) */

                      do
                        {

                          /* Prepare coupled calculation iteration if needed */

                          PrepareCCIter();

                          /* Make exponential fit for decay heat */

                          ExpoDecFit();

                          /* Run transport cycle */

                          TransportCycle();

                          /* Iterate coupled calculation routines */

                          IterateCC();

                          /* Repeat if needed */
                        }
                      while (RDB[DATA_ITERATE] == (double)YES);
                    }
                }

              /* Print material compositions (tää siirrettiin tuolta */
              /* transportcyclen edestä tähän 5.12.2012 / 2.1.10)    */

              PrintCompositions((long)RDB[DATA_BURN_STEP]);

              /* Start burnup timers */

              ResetTimer(TIMER_BURNUP);
              StartTimer(TIMER_BURNUP);
              StartTimer(TIMER_BURNUP_TOTAL);

              /* Print output-interval flag */

              if ((long)RDB[DATA_BURN_DECAY_CALC] == NO)
                {
                  if ((long)RDB[DATA_BURN_PRINT_STEP] == BURN_OUT_PRINT_FINAL)
                    {
                      if (step == (long)RDB[dep + DEP_HIS_N_STEPS])
                        WDB[DATA_BURN_PRINT_OUTPUT] = (double)YES;
                      else
                        WDB[DATA_BURN_PRINT_OUTPUT] = (double)NO;
                    }
                  else if ((long)RDB[DATA_BURN_PRINT_STEP] ==
                           BURN_OUT_PRINT_ALL)
                    WDB[DATA_BURN_PRINT_OUTPUT] = (double)YES;
                  else if ((long)RDB[DATA_BURN_PRINT_STEP] ==
                           BURN_OUT_PRINT_NONE)
                    WDB[DATA_BURN_PRINT_OUTPUT] = (double)NO;
                  else
                    Die(FUNCTION_NAME, "Invalid output print-interval");
                }
              else
                {
                  if (step == RDB[dep + DEP_HIS_N_STEPS])
                    WDB[DATA_BURN_PRINT_OUTPUT] = (double)YES;
                  else
                    WDB[DATA_BURN_PRINT_OUTPUT] = (double)NO;
                }

              /* output only at predictor (corresponding to BOS and earlier) */

              if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
                {
                  /* Write binary depletion file */

                  WriteDepFile();

                  /* Refresh inventory (add nuclides in case top inventory */
                  /* option is used) */

                  if (((long)RDB[DATA_BURN_STEP] > 0) ||
                      ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
                    RefreshInventory();

                  /* Print depletion output */

                  if ((long)RDB[DATA_BURN_PRINT_OUTPUT] == YES)
                    {
                      fprintf(outp, "Writing depletion output...\n");

                      PrintDepOutput();

                      fprintf(outp, "OK.\n\n");
                    }
                }

              /* Add step to SIE counter */

              WDB[DATA_BURN_STEP_TOT] = RDB[DATA_BURN_STEP_TOT] + 1.0;

              /* Break here if final step */

              if (step == RDB[dep + DEP_HIS_N_STEPS])
                {
                  /* Stop burnup timers */

                  StopTimer(TIMER_BURNUP);
                  StopTimer(TIMER_BURNUP_TOTAL);

                  /* Break loop */

                  break;
                }

              /* Calculate coefficients for the fit to xs/flux/power */

              DepletionPolyFit(dep, step);

              /* Set depletion step size */

              SetDepStepSize(dep, step);

              /* Burnup calculation */

              BurnMaterials(dep, step);

              /* Collect material compositions from MPI parallel tasks */

              CollectBurnData();

              /* Copy compositions to parent materials (tää siirrettiin */
              /* burnmaterial.c:n lopusta tähän 15.7.2013 / 2.1.15 että */
              /* kaikilla MPI taskeilla olisi käytössään sama data). */

              SumDivCompositions();

              /* Stop burnup timers */

              StopTimer(TIMER_BURNUP);
              StopTimer(TIMER_BURNUP_TOTAL);

              /* Add to number of predictor and corrector cycles */

              if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
                WDB[DATA_BURN_PRED_STEP] = RDB[DATA_BURN_PRED_STEP] + 1.0;
              else
                WDB[DATA_BURN_CORR_STEP] = RDB[DATA_BURN_CORR_STEP] + 1.0;

              /* Check for iterating the corrector */

              /* Update CI stopping criterion*/
              if (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) && ((long)RDB[DATA_BURN_SIE] == YES))
                StopCI();

              /* If further iterations are needed */
              if(((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
                 ((long)RDB[DATA_BURN_CI_LAST] == NO) &&
                 ((long)RDB[DATA_BURN_CI_MAXI] > 1.0) )
                {
                  /* increment iteration count (indexed 0,1,2...) */

                  WDB[DATA_BURN_CI_I] = RDB[DATA_BURN_CI_I] + 1.0;

                  /* Repeat corrector */
                  pc--;

                }

              /* Break from predictor corrector loop if we should not run */
              /* corrector */

              if (NoCorrector())
                {
                  /* Signal the external program about moving */
                  /* to next burnup point */

                  if(RDB[DATA_BURN_CI_I] == 0)
                    SignalExternal(SIGUSR2);

                  /* Read data interfaces for end-of-step */

                  ReadDataInterfaces();

                  break;
                }
            }

          /* Update cumulative burnup and time */

          WDB[DATA_BURN_CUM_BURNTIME0] = RDB[DATA_BURN_CUM_BURNTIME];
          WDB[DATA_BURN_CUM_BURNTIME] = RDB[DATA_BURN_CUM_BURNTIME]
            + RDB[DATA_BURN_TIME_INTERVAL];
          WDB[DATA_BURN_CUM_BURNUP] = RDB[DATA_BURN_CUM_BURNUP]
            +  RDB[DATA_BURN_BURNUP_INTERVAL];

          /* Update burnup step (total, not in this interval) */

          WDB[DATA_BURN_STEP] = RDB[DATA_BURN_STEP] + 1.0;
        }

      /* Check for history variation break point */

      if (TestHisvBreak() == YES)
        {
          /* Calculate activities */

          CalculateActivities();

          /* Write binary depletion file */

          WriteDepFile();

          fprintf(outp, "Breaking burnup cycle for variations...\n\n");

          /* Break loop */

          break;
        }

      /* Next interval */

      dep = NextItem(dep);
    }

  /* Signal externally coupled program to end calculation */

  SignalExternal(SIGTERM);

  /* Check total time */

  if (RDB[DATA_BURN_CUM_BURNTIME] == 0.0)
    Die(FUNCTION_NAME, "No burnup calculation performed");
}

/*****************************************************************************/
