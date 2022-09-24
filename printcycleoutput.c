/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printcycleoutput.c                             */
/*                                                                           */
/* Created:       2011/04/03 (JLe)                                           */
/* Last modified: 2019/03/13 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Prints cycle-wise data & info to standard output             */
/*                                                                           */
/* Comments:  - Kytketyssä laskennassa tätä pitää päivittää                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintCycleOutput:"

/*****************************************************************************/

void PrintCycleOutput()
{
  long i, cycles, skip, skip1, pop, ptr, tb, tbmax, tme;
  double estimt, tott, tming, tmaxg;
  char tmpstr[MAX_STR];

  /* Check mpi task */

  if (mpiid > 0)
    return;

  if ((long)RDB[DATA_RMTX_MFP_CALC] == YES)
    return;

  /* Get cycle index, number of cycles and population */

  i = (long)RDB[DATA_CYCLE_IDX] + 1;

  /* Set number of inactive batches */

  if(RDB[DATA_USE_FSP] == (double)NO)
    {
      /* No fission source passing*/

      skip = (long)RDB[DATA_CRIT_SKIP];
      skip1 = skip;
    }
  else if ((RDB[DATA_BURN_STEP] + RDB[DATA_SOL_REL_ITER] == 0.0)
         && (RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP))
    {
      /* First transportcycle with fission source passing */

      skip = (long)RDB[DATA_CRIT_SKIP];
      skip1 = skip;
    }
  else
    {
      /* Subsequent transportcycle with fission source passing */

      skip = (long)RDB[DATA_CRIT_SKIP];
      skip1 = (long)RDB[DATA_FSP_CRIT_SKIP];
    }

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {

#ifdef MPI_MODE2

      /* Check number of tasks */

      if (mpitasks > 1)
        {
          /* Calculate number of particles per task */

          cycles = (long)(RDB[DATA_CRIT_CYCLES]/((double)mpitasks));

        }
      else
        cycles = (long)RDB[DATA_CRIT_CYCLES];

#else

      cycles = (long)RDB[DATA_CRIT_CYCLES];

#endif

      pop = (long)RDB[DATA_CRIT_POP];
    }
  else
    {
      cycles = (long)RDB[DATA_SRC_BATCHES];
      pop = (long)RDB[DATA_SRC_POP];
    }

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
    {
      /* Get time bin index */

      tb = (long)RDB[DATA_DYN_TB];
      tbmax = (long)RDB[DATA_DYN_NB];

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

      tb = 0;
      tbmax = 1;

      tming = -INFTY;
      tmaxg = INFTY;
    }

  /* Get estimater running time */

  EstimateRuntime();

  estimt = RDB[DATA_ESTIM_CYCLE_TIME];
  tott = RDB[DATA_ESTIM_TOT_TIME];

  /***************************************************************************/

  /***** Print inactive cycle output *****************************************/

  if (i < skip + 1)
    {
      if (RDB[DATA_USE_FSP] == (double)NO)
        fprintf(outp, "Inactive cycle %3ld / %3ld: ", i, skip);
      else
        fprintf(outp, "Inactive cycle %3ld / %3ld: ", i + skip1 - skip, skip1);

      if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
        fprintf(outp, "k-eff = %1.5f\n", RDB[DATA_WIELANDT_KP]);
      else
        fprintf(outp, "k-eff = %1.5f\n", RDB[DATA_CYCLE_KEFF]);
    }

  /***************************************************************************/

  /***** Print active cycle output *******************************************/

  if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT) &&
      (i == skip + 1))
    fprintf(outp, "\n----- Begin active cycles -----\n\n");

  if (i > skip)
    {
      fprintf(outp, "------------------------------------------------------------\n");

      fprintf(outp, "\nSerpent %s", CODE_VERSION);

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        fprintf(outp, " -- Static criticality source simulation\n");
      else if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
        {
          if ((long)RDB[DATA_DYN_NB] == 1)
            fprintf(outp, " -- Neutron external source simulation\n");
          else
            fprintf(outp, " -- Dynamic neutron external source simulation\n");
        }
      else if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
        fprintf(outp, " -- Photon transport simulation\n");

      if ((long)RDB[DATA_PTR_TITLE] > VALID_PTR)
        fprintf(outp, "\nTitle: \"%s\" ", GetText(DATA_PTR_TITLE));
      else
        fprintf(outp, "\nInput file: \"%s\" ", GetText(DATA_PTR_INPUT_FNAME));

      if ((long)RDB[DATA_CONFIDENTIAL] == YES)
        fprintf(outp, "(CONFIDENTIAL)\n");
      else
        fprintf(outp, "\n");

      if ((long)RDB[DATA_RUN_VR_ITER] == YES)
        fprintf(outp, "\nVR iteration: cycle = %ld / %ld\n",
                (long)RDB[DATA_VR_ITER_IDX] + 1,
                (long)RDB[DATA_TOT_VR_ITER]);
      if ((long)RDB[DATA_COEF_CALC_IDX] > 0)
        fprintf(outp, "\nCoefficient calculation: restart = %ld / %ld\n",
                (long)RDB[DATA_COEF_CALC_RUN_IDX],
                (long)RDB[DATA_COEF_CALC_TOT_RUNS]);
      else if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES) &&
               ((long)RDB[DATA_BURN_PTR_DEP] > VALID_PTR))
        {
          fprintf(outp, "\nTransport calculation: step = %ld / %ld ",
                 (long)RDB[DATA_BURN_STEP] + 1,
                 (long)RDB[DATA_BURN_TOT_STEPS] + 1);

          if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_NONE)
            fprintf(outp, "\n");
          else if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
            fprintf(outp, "(predictor)\n");

          else if ((long)RDB[DATA_BURN_CI_MAXI] > 1)

            fprintf(outp, "(corrector %ld/%ld)\n",
                    (long)RDB[DATA_BURN_CI_I]+1, (long)RDB[DATA_BURN_CI_MAXI]);

          else
            fprintf(outp, "(corrector)\n");

          if (RDB[DATA_INI_FMASS] > 0.0)
            fprintf(outp, "                       BU   = %1.2f MWd/kgU\n",
                   RDB[DATA_BURN_CUM_BURNUP]);

          if (RDB[DATA_BURN_CUM_BURNTIME] == 0.0)
            fprintf(outp, "                       time = 0.00 days\n");
          else
            fprintf(outp, "                       time = %s\n",
                   TimeIntervalStr(RDB[DATA_BURN_CUM_BURNTIME]));
        }

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        {
          if ((long)RDB[DATA_GROW_POP_SIM] == NO)
            fprintf(outp, 
                    "\nActive cycle %4ld / %ld  Source neutrons : %5ld\n\n",
                    i - skip, cycles, (long)RDB[DATA_CYCLE_BATCH_SIZE]);
          
          else
            {
              fprintf(outp, "\nSimulated histories per cycle   : %ld\n",
                      (long)RDB[DATA_CYCLE_BATCH_SIZE]);
              fprintf(outp, "Fraction of histories completed : %1.1f%% \n\n", 
                      100.0*RDB[DATA_NHIST_CYCLE]/RDB[DATA_GROW_POP_MAX_HIS]);
            }

          if(RDB[DATA_RUN_CC] == (double)YES)
            {
              fprintf(outp, "Coupled calculation iteration:");
              fprintf(outp, "%3ld / %3ld \n\n",
                      (long)RDB[DATA_SOL_REL_ITER] + 1,
                      (long)RDB[DATA_SOL_REL_MAX_ITER]);
            }
        }
      else
        {
          fprintf(outp, "\nSource batch %2ld / %ld (%ld histories per batch)\n\n",
                  i, cycles, pop);


          if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
            {
              /* Print current time interval */

              fprintf(outp, "Time interval %ld / %ld from %6E s to %6E s\n\n", 
                      tb + 1, tbmax, tming, tmaxg);
            }

        }

      /* Print running time */

      fprintf(outp, "Running time :                %9s\n",
              TimeStr((long)TimerVal(TIMER_RUNTIME)));

      if (tott > 0.0)
        {
          if (TimerVal(TIMER_TRANSPORT_ACTIVE) > 2.0)
            {
              strcpy(tmpstr, TimeStr((long)tott));
              fprintf(outp, "Estimated running time :      %9s %9s\n",
                     TimeStr((long)estimt), tmpstr);

              strcpy(tmpstr, TimeStr((long)(tott- TimerVal(TIMER_RUNTIME))));
              fprintf(outp, "Estimated running time left : %9s %9s\n",
                     TimeStr((long)(estimt - TimerVal(TIMER_RUNTIME))),
                     tmpstr);
            }
          else
            {
              fprintf(outp, "Estimated running time :        -:--:--   -:--:--\n");
              fprintf(outp, "Estimated running time left :   -:--:--   -:--:--\n");
            }
        }
      else
        {
          if (TimerVal(TIMER_TRANSPORT_ACTIVE) > 2.0)
            {
              fprintf(outp, "Estimated running time :      %9s\n",
                     TimeStr((long)estimt));
              fprintf(outp, "Estimated running time left : %9s\n",
                     TimeStr((long)(estimt - TimerVal(TIMER_RUNTIME))));
            }
          else
            {
              fprintf(outp, "Estimated running time :        -:--:--\n");
              fprintf(outp, "Estimated running time left :   -:--:--\n");
            }
        }

      fprintf(outp, "\nEstimated relative CPU usage : %7.1f%%\n",
             100.0*TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
             TimerVal(TIMER_TRANSPORT_CYCLE));

      /* K-eff estimates */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        {
          ptr = (long)RDB[RES_ANA_KEFF];
          fprintf(outp, "\nk-eff (analog)    = %1.5f +/- %1.5f  [%1.5f  %1.5f]\n",
                 Mean(ptr, 0), StdDev(ptr, 0),
                 Mean(ptr, 0) - 1.96*StdDev(ptr, 0),
                 Mean(ptr, 0) + 1.96*StdDev(ptr, 0));

          if ((long)RDB[DATA_OPTI_IMPLICIT_RR] == YES)
            {
              ptr = (long)RDB[RES_IMP_KEFF];
              fprintf(outp, "k-eff (implicit)  = %1.5f +/- %1.5f  [%1.5f  %1.5f]\n",
                     Mean(ptr, 0), StdDev(ptr, 0),
                     Mean(ptr, 0) - 1.96*StdDev(ptr, 0),
                     Mean(ptr, 0) + 1.96*StdDev(ptr, 0));
            }
          else
            {
              ptr = (long)RDB[RES_COL_KEFF];
              fprintf(outp, "k-eff (collision) = %1.5f +/- %1.5f  [%1.5f  %1.5f]\n",
                     Mean(ptr, 0), StdDev(ptr, 0),
                     Mean(ptr, 0) - 1.96*StdDev(ptr, 0),
                     Mean(ptr, 0) + 1.96*StdDev(ptr, 0));
            }
        }
      else if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
        {
          ptr = (long)RDB[RES_SRC_MULT];
          fprintf(outp, "\nMultiplication   = %1.4E (%1.5f)\n",
                  Mean(ptr, 0), RelErr(ptr, 0));

          ptr = (long)RDB[RES_ANA_KEFF];
          fprintf(outp, "k-eff (analog)   = %1.5f +/- %1.5f  [%1.5f  %1.5f]\n",
                  Mean(ptr, 0), StdDev(ptr, 0),
                  Mean(ptr, 0) - 1.96*StdDev(ptr, 0),
                  Mean(ptr, 0) + 1.96*StdDev(ptr, 0));

          ptr = (long)RDB[RES_EXT_K];
          fprintf(outp, "k0    (source)   = %1.5f +/- %1.5f  [%1.5f  %1.5f]\n",
                  Mean(ptr, 0),
                  StdDev(ptr, 0),
                  Mean(ptr, 0) - 1.96*StdDev(ptr, 0),
                  Mean(ptr, 0) + 1.96*StdDev(ptr, 0));

          if (RDB[DATA_DYN_TMAX] != INFTY)
            {
              if (RDB[DATA_DYN_DT] > ZERO)
                fprintf(outp, "\nTime cut-off at %1.2E seconds, dt = %s\n",
                        RDB[DATA_DYN_TMAX], TimeIntervalStr(RDB[DATA_DYN_DT]));
              else if ((long)RDB[DATA_DYN_NB] > 1)
                fprintf(outp, "\nTime cut-off at %1.2E seconds, %ld intervals\n",
                        RDB[DATA_DYN_TMAX], (long)RDB[DATA_DYN_NB]);
              else
                fprintf(outp, "\nTime cut-off at %1.2E seconds\n",
                        RDB[DATA_DYN_TMAX]);
            }

        }

      fprintf(outp, "\n");

      /* Options */

#ifdef DEBUG

      fprintf(outp, "(DBG) ");

#endif

      if ((long)RDB[DATA_OPTI_REPLAY] == YES)
        fprintf(outp, "(R) ");

      fprintf(outp, "(O%ld) ", (long)RDB[DATA_OPTI_MODE]);

      if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
        fprintf(outp, "(W) ");

      if ((long)RDB[DATA_ITER_MODE] == ITER_MODE_ALBEDO)
        fprintf(outp, "(IA) ");

      if ((long)RDB[DATA_B1_CALC] == YES)
        fprintf(outp, "(B1) ");

      if ((long)RDB[DATA_OPT_IMPL_CAPT] == YES)
        fprintf(outp, "(IC) ");
      /*
      if ((long)RDB[DATA_OPT_IMPL_NXN] == YES)
        fprintf(outp, "(IX) ");

      if ((long)RDB[DATA_OPT_IMPL_FISS] == YES)
        fprintf(outp, "(IF) ");
      */
      if ((long)RDB[DATA_USE_URES] == YES)
        fprintf(outp, "(UNR) ");

      if ((long)RDB[DATA_USE_DBRC] == YES)
        fprintf(outp, "(DBRC) ");

      if ((long)RDB[DATA_TMS_MODE] != TMS_MODE_NONE)
        fprintf(outp, "(TMS) ");

      if ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] > -1)
        fprintf(outp, "(XE) ");

      if ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] > -1)
        fprintf(outp, "(SM) ");

      if ((long)RDB[DATA_PHOTON_PRODUCTION] != NO)
        fprintf(outp, "(N/P) ");

      if ((long)RDB[DATA_SENS_MODE] != SENS_MODE_NONE)
        fprintf(outp, "(SENS) ");

#ifdef MPI

      fprintf(outp, "(MPI=%d) ", mpitasks);

#endif

#ifdef OPEN_MP

      fprintf(outp, "(OMP=%ld) ", (long)RDB[DATA_OMP_MAX_THREADS]);

#endif

      if((long)RDB[DATA_GROW_POP_SIM] == YES)
        fprintf(outp, "(GPOP) ");

      if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
        fprintf(outp, "(DD) ");

      if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)
        {
          if ((long)RDB[DATA_BURN_PRED_TYPE] == PRED_TYPE_CONSTANT)
            fprintf(outp, "(CE");
          else
            fprintf(outp, "(LE");

          if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_CONSTANT)
            fprintf(outp, "/CE");
          else if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_LINEAR)
            fprintf(outp, "/LI");
          else if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_QUADRATIC)
            fprintf(outp, "/QI");

          if ((long)RDB[DATA_OTF_BURN_MODE] == YES)
            fprintf(outp, " + OTF)");
          else if (((long)RDB[DATA_BURN_PRED_NSS] > 1) &&
                   ((long)RDB[DATA_BURN_CORR_NSS] < 1))
            fprintf(outp, " + %ld SS) ", (long)RDB[DATA_BURN_PRED_NSS]);
          else if (((long)RDB[DATA_BURN_PRED_NSS] == 1) &&
                   ((long)RDB[DATA_BURN_CORR_NSS] > 1))
            fprintf(outp, " + %ld SS) ", (long)RDB[DATA_BURN_CORR_NSS]);
          else if (((long)RDB[DATA_BURN_PRED_NSS] > 1) &&
                   ((long)RDB[DATA_BURN_CORR_NSS] > 1))
            fprintf(outp, " + %ld/%ld SS) ", (long)RDB[DATA_BURN_PRED_NSS],
                   (long)RDB[DATA_BURN_CORR_NSS]);
          else
            fprintf(outp, ") ");
        }

      if ((long)RDB[DATA_RUN_CC] == YES)
        fprintf(outp, "(CC) ");

      if ((long)RDB[DATA_RMTX_CALC] == YES)
        fprintf(outp, "(WWG) ");

      if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == YES)   
        fprintf(outp, "(WWIN) ");

      fprintf(outp, "\n");

      fprintf(outp, "------------------------------------------------------------\n");

    }

  /* Check if completed */

  if ((i == cycles + skip) &&
      ((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_DYN))
    {
      fprintf(outp, "\nTransport cycle completed in %s.\n\n",
              TimeIntervalStr(TimerVal(TIMER_TRANSPORT)));

      PrintTMSDiagnostics();
    }

  /* For EDo (18.1.2014) */

  DiffCoefED(5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

  /***************************************************************************/
}

/*****************************************************************************/
