/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : transportcycle.c                               */
/*                                                                           */
/* Created:       2011/05/23 (JLe)                                           */
/* Last modified: 2020/06/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Prepares and runs the main transport cycle                   */
/*                                                                           */
/* Comments: - Tää on ihan kauhea sekasotku!!!                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TransportCycle:"

/*****************************************************************************/

void TransportCycle()
{
  long nb, nb0, maxb, skip, nn, nt, maxt, id, idx, ptr, tme, loc0, sz;
  long nsrc, tosimulate, sign;
  long no_particles_left, maybe_done, done;
  double t0, c0, tmax;

  /***************************************************************************/

  /***** Main transport cycle ************************************************/

  /* Check norun option */

  if ((long)RDB[DATA_STOP_AFTER_PROCESSING] == YES)
    {
      fprintf(outp, "Processing completed, ");
      fprintf(outp, "calculation terminated before transport cycle.\n\n");
      exit(0);
    }

  /* Check if starting batch or history is given */

  if ((nb0 = (long)RDB[DATA_NHIST_BATCH0]) != 0)
    {
      /* Check mode */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        Note(0, "Starting batch has no effect in criticality source mode");
      else if (nb0 > 0)
        {
          /* Print warning */

          Note(0, "Starting simulation from batch %ld", nb0);

          /* Set value */

          WDB[DATA_NHIST_TOT] = ((double)nb0)*RDB[DATA_SRC_POP];
        }
      else
        {
          /* Print warning */

          Note(0, "Starting simulation from history %ld", -nb0);

          /* Set value */

          WDB[DATA_NHIST_TOT] = -((double)nb0);
        }
    }

  /* Avoid compiler warning */

  loc0 = -1;
  sz = -1;

  /* Reset total and active timer */

  ResetTimer(TIMER_TRANSPORT);
  ResetTimer(TIMER_TRANSPORT_ACTIVE);

  /* Start transport timer */

  StartTimer(TIMER_TRANSPORT);
  StartTimer(TIMER_TRANSPORT_TOTAL);

  /* Reset completed flag and sort counter */

  WDB[DATA_SIMULATION_COMPLETED] = (double)NO;
  WDB[DATA_SORT_COUNT] = 1.0;

  /* Put on-the-fly burnup concentrations */

  PutOTFBurnConc();

  /* Set maximum running time */

  tmax = RDB[DATA_MAX_TRANSPORT_RUNTIME];

  /* Check mode */

  if ((long)RDB[DATA_SRC_IMP_CALC] == YES)
    {
      /***********************************************************************/

      /***** Get source point importances ************************************/

      fprintf(outp,
              "Sampling %ld source points to find maximum importance...\n",
              (long)RDB[DATA_SRC_IMP_NPS]);

#ifdef OPEN_MP
#pragma omp parallel private(id, idx, nn, ptr)
#endif
      {
        /* Get Open MP thread id */

        id = OMP_THREAD_NUM;

#ifdef OPEN_MP
#pragma omp for schedule(dynamic)
#endif
        /* Loop over source neutrons */

        for (nn = 0; nn < (long)RDB[DATA_SRC_IMP_NPS]; nn++)
          {
            /* Calculate particle index */

            idx = (long)RDB[DATA_NHIST_TOT];
            idx = idx + nn;

            /* Sample source point */

            if (SampleSrcPoint(id, nn, idx) > VALID_PTR)
              {
                /* Get particle from que and put it back */

                ptr = FromQue(id);
                ToStack(ptr, id);
              }
          }
      }

      /* Print results (terminates run) */

      PrintSrcImp();

      /***********************************************************************/
    }

  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC)
    {
      /***********************************************************************/

      /***** External source simulation **************************************/

      /* Reset skip cycles */

      WDB[DATA_CRIT_SKIP] = 0.0;

      /* Set cycle-wise batch size */

      WDB[DATA_CYCLE_BATCH_SIZE] = RDB[DATA_SRC_POP];

      /* Set number of source points in batch */

      WDB[DATA_NHIST_CYCLE] = RDB[DATA_SRC_POP];

      /* Get number of bins */

      maxt = (long)RDB[DATA_DYN_NB];
      CheckValue(FUNCTION_NAME, "maxt", "", maxt, 1, 1000000000);

      /* Get pointer to time bin structure */

      tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);

      /* Get pointer to bins */

      tme = (long)RDB[tme + TME_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);

      /* Start active transport timer */

      StartTimer(TIMER_TRANSPORT_ACTIVE);

      /* Start cycle-wise transport timer */

      ResetTimer(TIMER_TRANSPORT_CYCLE);
      StartTimer(TIMER_TRANSPORT_CYCLE);

      /* Reset batch counters */

      WDB[DATA_BATCH_COUNT] = 0.0;
      WDB[DATA_MICRO_CALC_BATCH_COUNT] = 0.0;
      WDB[DATA_MICRO_CALC_BATCH_NUM] = 0.0;

      /* Check batching interval */

      if ((long)RDB[DATA_SRC_BATCHES] % (long)RDB[DATA_BATCH_INTERVAL])
        Error(0,
            "Total number of batches %ld is not a multiple of interval %ld",
            (long)RDB[DATA_SRC_BATCHES], (long)RDB[DATA_BATCH_INTERVAL]);

      /* Start simulation */

      if ((long)RDB[DATA_RMTX_MFP_CALC] == NO)
        fprintf(outp, "Starting external source simulation...\n\n");

      /* Loop over batches */

      for (nb = 0; nb < (long)RDB[DATA_SRC_BATCHES]; nb++)
        {
          /*******************************************************************/

          /***** Loop over source batches ************************************/

          /* Reset time bin index */

          WDB[DATA_DYN_TB] = 0.0;

          /* Set previous population size in dynamic calculation */

          WDB[DATA_DYN_SRC_PREV_POP] = RDB[DATA_SRC_POP];

          /* Set time cut-off for first cycle (intervallimoodissa   */
          /* noi on 0 ja se mitä set tcut -parametrilla on annettu) */

          WDB[DATA_TIME_CUT_TMIN] = RDB[tme];
          WDB[DATA_TIME_CUT_TMAX] = RDB[tme + 1];

          /* Stop tracks at outer boundary if time cut-off is used */

          if (RDB[DATA_TIME_CUT_TMAX] < INFTY)
            WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

          /* Reset cycle k-eff and put starting weight */

          WDB[DATA_CYCLE_KEFF] = 1.0;
          WDB[DATA_DYN_WGT0] = RDB[DATA_SRC_POP];

          /* Put cycle index */

          WDB[DATA_CYCLE_IDX] = (double)nb;

          /* Get beginning time */

          t0 = TimerVal(TIMER_TRANSPORT);
          c0 = TimerCPUVal(TIMER_TRANSPORT);

          /* Reset completed flags */

          OMPResetComp(-1);

          /* Start parallel timer */

          StartTimer(TIMER_OMP_PARA);

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

            for (nn = 0; nn < (long)RDB[DATA_SRC_POP]; nn++)
              {
                /* Calculate particle index */

                idx = (long)RDB[DATA_NHIST_TOT];
                idx = idx + (long)RDB[DATA_SRC_POP]*nb + nn;

                /* Sample source particle and track */

                if (SampleSrcPoint(id, nn, idx) > VALID_PTR)
                  Tracking(id);

                /* Check if common que is used */

                if ((long)RDB[DATA_COMMON_QUE_LIM] != 0)
                  while (OMPTestComp() == NO)
                    {
                      /* Call tracking routine */

                      Tracking(id);

                      /* Set completed flag */

                      OMPSetComp(id);
                    }
              }
          }

          /* Move collected pulse data to statistics */

          if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
            PulseDet(-1, -1, -1.0, 0.0, 0.0, 0.0, -1.0, -1);

          /* Stop parallel timer */

          StopTimer(TIMER_OMP_PARA);

          /* Collect data from interval */

          CollectDynData();

          /* Reset normalization coefficients (pitää kutsua tässä) */

          WDB[DATA_NORM_COEF_N] = -1.0;
          WDB[DATA_NORM_COEF_G] = -1.0;

          /* Update batch count (JLe: siirretty tuolta alempaa tänne */
          /* 1.9.2015 / 2.1.24) */

          WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_COUNT] + 1.0;

          /* Check for remaining time intervals */

          if (maxt > 1)
            {
              /* Fix normalization to first step */

              if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
                NormCoef(PARTICLE_TYPE_NEUTRON);
              if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
                NormCoef(PARTICLE_TYPE_GAMMA);

              /* Re-open buffer for writing */

              WDB[DATA_BUF_REDUCED] = (double)NO;

              /* Loop over remaining time intervals */

              for (nt = 1; nt < maxt; nt++)
                {
                  /* Re-open buffer for writing */

                  WDB[DATA_BUF_REDUCED] = (double)NO;

                  /* Set time bin index */

                  WDB[DATA_DYN_TB] = (double)nt;

                  /* Set new time cut-off (not set if interval mode) */

                  if (RDB[DATA_DYN_DT] < ZERO)
                    {
                      WDB[DATA_TIME_CUT_TMIN] = RDB[tme + nt];
                      WDB[DATA_TIME_CUT_TMAX] = RDB[tme + nt + 1];
                    }

                  /* Normalize source */

                  if (NormalizeDynSrc() < 0)
                    break;

                  /* Start parallel timer */

                  StartTimer(TIMER_OMP_PARA);

                  /* Loop until source is empty */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
                  {
                    /* Get Open MP thread id */

                    id = OMP_THREAD_NUM;

                    /* Loop over source */

                    while (FromSrc(id) > VALID_PTR)
                      Tracking(id);
                  }

                  /* Collect data from interval */

                  CollectDynData();

                  /* Stop parallel timer */

                  StopTimer(TIMER_OMP_PARA);
                }
            }

          /* Get end time */

          t0 = TimerVal(TIMER_TRANSPORT) - t0;
          c0 = TimerCPUVal(TIMER_TRANSPORT) - c0;

          /* Score time */

          ptr = (long)RDB[RES_CYCLE_RUNTIME];
          AddStat(t0, ptr, 0);
          AddStat(c0, ptr, 1);

          /* Add to micro batch counter */

          WDB[DATA_MICRO_CALC_BATCH_COUNT] =
            RDB[DATA_MICRO_CALC_BATCH_COUNT] + 1.0;

          /* Check batch interval */

          if ((long)RDB[DATA_BATCH_COUNT] == (long)RDB[DATA_BATCH_INTERVAL])
            {
              /* Collect and clear buffered results */

              CollectResults();
              CollectDet();
              PoisonEq();
              SolveOTFBurn();
              CalcMicroGroupXS();
              ClearBuf();

              /* Reset batch counter */

              WDB[DATA_BATCH_COUNT] = 0.0;
            }

          /* Flush bank */

          FlushBank();

          /* Stop cycle-wise transport timer */

          StopTimer(TIMER_TRANSPORT_CYCLE);

          if (TimerVal(TIMER_TRANSPORT_CYCLE) > 0.0)
            t0 = TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
              TimerVal(TIMER_TRANSPORT_CYCLE);
          else
            t0 = 0.0;

          /* CPU usage */

          ptr = (long)RDB[RES_CPU_USAGE];
          AddStat(t0, ptr, 0);

          /* Print cycle-wise output */

          PrintCycleOutput();

          /* Reset and restart cycle-wise transport timer */

          ResetTimer(TIMER_TRANSPORT_CYCLE);
          StartTimer(TIMER_TRANSPORT_CYCLE);

          /* Sort lists */

          SortAll();

          /* Print results */

          if (!((nb + 1) % (long)RDB[DATA_PRINT_INTERVAL]))
            {
              MatlabOutput();
              MicroDepOutput();
              DetectorOutput();
              MeshPlotter();
              PrintCoreDistr();
              PrintHistoryOutput();
              PrintPBData();
              PrintInterfaceOutput();
              PrintFinix();
              RROutput();
              /*
              GeometryPlotter(NO);
              */
              FissMtxOutput();
              MORAOutput();
              WriteICMData();
            }

          /* Add to time vector */

          ptr = (long)RDB[DATA_BTCH_TIME];
          WDB[ptr + nb] = TimerVal(TIMER_RUNTIME);

          /* Check maximum time */

          if ((tmax > 0.0) && (TimerVal(TIMER_RUNTIME) > tmax))
            {
              /* Print message */

              fprintf(outp, "\nMaximum running time %s reached.\n\n",
                      TimeStr((long)tmax));

              /* Terminate loop */

              break;
            }

          /******************************************************************/
        }

      /* Update number of neutron histories run */

      WDB[DATA_NHIST_TOT] = RDB[DATA_NHIST_TOT] +
        RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

      /* Stop cycle-wise transport timer */

      StopTimer(TIMER_TRANSPORT_CYCLE);

      /***********************************************************************/
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DELDYN)
    {
      /***********************************************************************/

      /***** Time dependent simulation with delayed neutron ******************/
      /* NB: This will be combined with SIMULATION_MODE_DYN                  */
      /* This is only set on, if a precursor detector is explicitly set up in*/
      /* input and the mode would otherwise be SIMULATION_MODE_SRC           */

      /* Reset skip cycles */

      WDB[DATA_CRIT_SKIP] = 0.0;

      /* Set number of batches */

#ifdef MPI_MODE1

      maxb = (long)RDB[DATA_SRC_BATCHES];
      WDB[DATA_SRC_POP] = RDB[DATA_SRC_POP]/(double)mpitasks;

#else

      maxb = (long)(RDB[DATA_SRC_BATCHES]/(double)mpitasks);

#endif

      /* Set cycle-wise batch size */

      WDB[DATA_CYCLE_BATCH_SIZE] = RDB[DATA_SRC_POP];

      /* Set number of source points in batch */

      WDB[DATA_NHIST_CYCLE] = RDB[DATA_SRC_POP];

      /* Get number of bins */

      maxt = (long)RDB[DATA_DYN_NB];
      CheckValue(FUNCTION_NAME, "maxt", "", maxt, 1, 1000000000);

      /* Get pointer to time bin structure */

      tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);

      /* Get pointer to bins */

      tme = (long)RDB[tme + TME_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);

      /* Start active transport timer */

      StartTimer(TIMER_TRANSPORT_ACTIVE);

      /* Start cycle-wise transport timer */

      ResetTimer(TIMER_TRANSPORT_CYCLE);
      StartTimer(TIMER_TRANSPORT_CYCLE);

      /* Start simulation */

      fprintf(outp, "Starting time dependent simulation with delayed neutrons...\n\n");

      /* Loop over batches */

      for (nb = 0; nb < maxb; nb++)
        {

          /* Reset time bin index */

          WDB[DATA_DYN_TB] = 0.0;

          /* Set time cut-off for first cycle */

          WDB[DATA_TIME_CUT_TMIN] = RDB[tme];
          WDB[DATA_TIME_CUT_TMAX] = RDB[tme + 1];

          /* Plot geometry for first time interval */

          if (nb == 0)
            GeometryPlotter(NO);

          /* Initialize values in precursor detectors */
          /* This initializes the values in buffers for the first */
          /* time bin, currently this has to be done separately for */
          /* each batch */

          InitPrecDet();

          /* Set simulation normalization and calculate weights */
          /* of live neutrons and precursors */

          NormalizePrecDet();

          /* Do population control for initial source */

          PrecursorPopControl();

          /* Handle decay of precursor during interval */
          /* and calculate weight to emit */

          DecayMeshPrecDet();

          /* Stop tracks at outer boundary if time cut-off is used */

          if (RDB[DATA_TIME_CUT_TMAX] < INFTY)
            WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

          /* Reset cycle k-eff and put starting weight */

          WDB[DATA_CYCLE_KEFF] = 1.0;
          WDB[DATA_DYN_WGT0] = RDB[DATA_SRC_POP];

          /* Set batch counters */

          WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_INTERVAL];
          WDB[DATA_MICRO_CALC_BATCH_COUNT] =
            RDB[DATA_MICRO_CALC_BATCH_SIZE];

          /* Put cycle index */

          WDB[DATA_CYCLE_IDX] = (double)nb;

          /* Sample delayed neutrons for the interval */

          SampleDelnu();

          /* Handle decay of pointwise precursors over interval */

          DecayPointPrecDet();

          /* Re-open buffer for writing */
          /* Was closed in sampledelnu  */

          WDB[DATA_BUF_REDUCED] = (double)NO;

          /* Get number of source neutrons */

          nsrc = (long)RDB[DATA_SRC_POP];
          sign = 1;

          /* Try to override with number of live neutrons */

          if ((ptr = (long)RDB[DATA_PTR_PREC_DET]) > VALID_PTR)
            {
              nsrc = (long)RDB[ptr + PRECDET_N_LIVE];
              sign = -1;
            }

          /* Get beginning time */

          t0 = TimerVal(TIMER_TRANSPORT);
          c0 = TimerCPUVal(TIMER_TRANSPORT);

          /* Start parallel timer */

          StartTimer(TIMER_OMP_PARA);

#ifdef DNPRINT
          fprintf(outp, "Sampling %ld live neutrons (printing from transportcycle.c)\n", nsrc);
#endif


          /* Parallel loop for sampling live source */
          /* TODO: set random number indexes in sampledelnu */

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

            for (nn = 0; nn < nsrc; nn++)
              {
                /* Calculate particle index */

                idx = (long)RDB[DATA_NHIST_TOT];
                idx = idx + (long)RDB[DATA_SRC_POP]*nb + nn;

                /* Sample source point */

                SampleSrcPoint(id, sign*nn, idx);
              }
          }

#ifdef DNPRINT
          fprintf(outp, "Moving to transport neutrons (printing from transportcycle.c)\n");
#endif

          /* Parallel loop over histories */
          /* Track neutrons */

          tosimulate = ReDistributeQues();

          while (tosimulate > 0)
            {

              /* Track current generation */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
              {
                /* Get Open MP thread id */

                id = OMP_THREAD_NUM;

                /* Track neutrons */

                Tracking(id);
              }

              /* Even out ques for next generation */

              tosimulate = ReDistributeQues();

            }

          /* Stop parallel timer */

          StopTimer(TIMER_OMP_PARA);

          /* Collect data from interval */

          CollectDynData();

          /* Get banked precursors */

          GetBankedPrecursors();

          /* Reset normalization coefficients if normalization was not */
          /* set by dynamic source */

          if ((long)RDB[DATA_PTR_PREC_DET] < VALID_PTR)
          {
            WDB[DATA_NORM_COEF_N] = -1.0;
            WDB[DATA_NORM_COEF_G] = -1.0;
          }

          /* Check for remaining time intervals */

          if (maxt > 1)
            {
              /* Fix normalization to first step */
              /* With dynamic source, normalization is fixed by source */

              if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
                NormCoef(PARTICLE_TYPE_NEUTRON);
              if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
                NormCoef(PARTICLE_TYPE_GAMMA);

              /* Re-open buffer for writing */

              WDB[DATA_BUF_REDUCED] = (double)NO;

              /* Loop over remaining time intervals */

              for (nt = 1; nt < maxt; nt++)
                {

                  /* Re-open buffer for writing */

                  WDB[DATA_BUF_REDUCED] = (double)NO;

                  /* Set time bin index */

                  WDB[DATA_DYN_TB] = (double)nt;

                  /* Set new time cut-off */

                  WDB[DATA_TIME_CUT_TMIN] = RDB[tme + nt];
                  WDB[DATA_TIME_CUT_TMAX] = RDB[tme + nt + 1];

                  /* Plot geometry for subsequent time interval */

                  if (nb == 0)
                    GeometryPlotter(NO);

                  /* Do population control on precursors */

                  PrecursorPopControl();

                  /* Count live neutrons and their weight */

                  CountDynSrc();

                  /* Handle decay of precursor during interval */
                  /* and calculate weight to emit from mesh    */

                  DecayMeshPrecDet();

                  /* Resize the number of live neutrons */
                  /* This is the new population control */

                  ResizeDynSrc();

                  /* Sample delayed neutron source for time interval */
                  /* After this we should have RDB[DATA_SRC_POP]     */
                  /* neutrons */

                  SampleDelnu();

                  /* Re-open buffer for writing (was closed in sampledelnu) */

                  WDB[DATA_BUF_REDUCED] = (double)NO;

                  /* Handle decay of pointwise precursors over interval */

                  DecayPointPrecDet();

                  /* Normalize source */
                  /* Should not do population control with delayed neutrons */
                  /* as the number of neutrons is already correct */

                  NormalizeDynSrc();

                  /* Get number of particles to simulate */

                  tosimulate = ReDistributeQues();

#ifdef DNPRINT
          fprintf(outp, "Moving to transport neutrons (printing from transportcycle.c)\n");
#endif

                  /* Start parallel timer */

                  StartTimer(TIMER_OMP_PARA);

                  /* Loop until source is empty */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
                  {
                    /* Get Open MP thread id */

                    id = OMP_THREAD_NUM;

                    /* Loop over source (first generation) */

                    while(FromSrc(id) > VALID_PTR)
                      Tracking(id);

                  }

                  /* Parallel loop over histories */
                  /* Track neutrons */

                  /* Even out ques for second generation */

                  tosimulate = ReDistributeQues();

                  while (tosimulate > 0)
                    {

                      /* Track current generation */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
                      {
                        /* Get Open MP thread id */

                        id = OMP_THREAD_NUM;

                        /* Track neutrons */

                        Tracking(id);
                      }

                      /* Even out ques for next generation */

                      tosimulate = ReDistributeQues();

                    }

                  /* Collect data from interval */

                  CollectDynData();

                  /* Get banked precursors */

                  GetBankedPrecursors();

                  /* Stop parallel timer */

                  StopTimer(TIMER_OMP_PARA);
                }
            }

          /* Get end time */

          t0 = TimerVal(TIMER_TRANSPORT) - t0;
          c0 = TimerCPUVal(TIMER_TRANSPORT) - c0;

          /* Score time */

          ptr = (long)RDB[RES_CYCLE_RUNTIME];
          AddStat(t0, ptr, 0);
          AddStat(c0, ptr, 1);

          /* Collect and clear buffered results */

          CollectResults();
          CollectDet();
          CollectPrecDet();
          PoisonEq();
          SolveOTFBurn();
          CalcMicroGroupXS();
          ClearBuf();

          /* Write savesrc data if requested */

          WriteDynSrc();

          /* Flush precursor source */

          FlushPrecSource();

          /* Flush bank */

          FlushBank();

          /* Stop cycle-wise transport timer */

          StopTimer(TIMER_TRANSPORT_CYCLE);

          if (TimerVal(TIMER_TRANSPORT_CYCLE) > 0.0)
            t0 = TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
              TimerVal(TIMER_TRANSPORT_CYCLE);
          else
            t0 = 0.0;

          /* CPU usage */

          ptr = (long)RDB[RES_CPU_USAGE];
          AddStat(t0, ptr, 0);

          /* Print cycle-wise output */

          PrintCycleOutput();

          /* Reset and restart cycle-wise transport timer */

          ResetTimer(TIMER_TRANSPORT_CYCLE);
          StartTimer(TIMER_TRANSPORT_CYCLE);

          /* Sort lists */

          SortAll();

          /* Print results */

          if (!((nb + 1) % (long)RDB[DATA_PRINT_INTERVAL]))
            {
              MatlabOutput();
              MicroDepOutput();
              DetectorOutput();
              MeshPlotter();
              PrintCoreDistr();
              PrintHistoryOutput();
              PrintPBData();
              PrintInterfaceOutput();
              PrintFinix();
              PrintPrecDet();
              RROutput();
              /*
              GeometryPlotter(NO);
              */
              FissMtxOutput();
              MORAOutput();
              WriteICMData();
            }
        }

      /* Update number of neutron histories run */

      WDB[DATA_NHIST_TOT] = RDB[DATA_NHIST_TOT] +
        RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

      /* Stop cycle-wise transport timer */

      StopTimer(TIMER_TRANSPORT_CYCLE);

      /***********************************************************************/
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
    {
      /***********************************************************************/

      /***** Dynamic external source simulation ******************************/

      /* Reset skip cycles */

      WDB[DATA_CRIT_SKIP] = 0.0;

      /* Set number of batches */

#ifdef MPI_MODE1

      maxb = (long)RDB[DATA_SRC_BATCHES];
      WDB[DATA_SRC_POP] = RDB[DATA_SRC_POP]/(double)mpitasks;

#else

      maxb = (long)(RDB[DATA_SRC_BATCHES]/(double)mpitasks);

#endif

      /* Set cycle-wise batch size */

      WDB[DATA_CYCLE_BATCH_SIZE] = RDB[DATA_SRC_POP];

      /* Set number of source points in batch */

      WDB[DATA_NHIST_CYCLE] = RDB[DATA_SRC_POP];

      /* Get number of bins */

      maxt = (long)RDB[DATA_DYN_NB];
      CheckValue(FUNCTION_NAME, "maxt", "", maxt, 1, 1000000000);

      /* Get pointer to time bin structure */

      tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);

      /* Get pointer to bins */

      tme = (long)RDB[tme + TME_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);

      /* Set time cut-off for first cycle */

      WDB[DATA_TIME_CUT_TMIN] = RDB[tme];
      WDB[DATA_TIME_CUT_TMAX] = RDB[tme + 1];

      /* Reset batch counters */

      WDB[DATA_BATCH_COUNT] = 0.0;
      WDB[DATA_MICRO_CALC_BATCH_COUNT] = 0.0;
      WDB[DATA_MICRO_CALC_BATCH_NUM] = 0.0;

      /* Check batching interval */

      if ((long)RDB[DATA_SRC_BATCHES] % (long)RDB[DATA_BATCH_INTERVAL])
        Error(0,
            "Total number of batches %ld is not a multiple of interval %ld",
            (long)RDB[DATA_SRC_BATCHES], (long)RDB[DATA_BATCH_INTERVAL]);

      /* Reset time bin index */

      WDB[DATA_DYN_TB] = 0.0;

      /* Plot geometry for first time interval */

      GeometryPlotter(NO);

      /* Stop tracks at outer boundary if time cut-off is used */

      if (RDB[DATA_TIME_CUT_TMAX] < INFTY)
        WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

      /* Reset cycle k-eff and put starting weight */

      WDB[DATA_CYCLE_KEFF] = 1.0;
      WDB[DATA_DYN_WGT0] = RDB[DATA_SRC_POP];

      /* Reset counter for solution relaxation */

      WDB[DATA_SOL_REL_ITER] = 0.0;

      WDB[DATA_SOL_REL_NTOT] = 0.0;
      WDB[DATA_SOL_REL_NCUR] = RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];
      WDB[DATA_SOL_REL_N1] = RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

      /* Start simulation */

      fprintf(outp, "Starting time dependent simulation...\n\n");

      /* Start active transport timer */

      StartTimer(TIMER_TRANSPORT_ACTIVE);

      /* Start cycle-wise transport timer */

      ResetTimer(TIMER_TRANSPORT_CYCLE);
      StartTimer(TIMER_TRANSPORT_CYCLE);

      do
        {
          /* Prepare coupled calculation iteration */

          PrepareCCIter();

          /* Loop over batches */

          for (nb = 0; nb < maxb ; nb++)
            {

              /* Put cycle index */

              WDB[DATA_CYCLE_IDX] = (double)nb;

              /* Initialize values in precursor detectors */
              /* This initializes the values in buffers for the first */
              /* time bin, currently this has to be done separately for */
              /* each batch */

              InitPrecDet();

              /* Set simulation normalization and calculate weights */
              /* of live neutrons and precursors */

              NormalizePrecDet();

              /* Do population control for initial source */

              PrecursorPopControl();

              /* Handle decay of precursor during interval */
              /* and calculate weight to emit */

              DecayMeshPrecDet();

              /* Stop tracks at outer boundary if time cut-off is used */

              if (RDB[DATA_TIME_CUT_TMAX] < INFTY)
                WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

              /* Sample delayed neutrons for the interval */

              SampleDelnu();

              /* Handle decay of pointwise precursors over interval */

              DecayPointPrecDet();

              /* Re-open buffer for writing */
              /* Was closed in sampledelnu  */

              WDB[DATA_BUF_REDUCED] = (double)NO;

              /* Get number of source neutrons */

              nsrc = (long)RDB[DATA_SRC_POP];

              /* Try to override with number of live neutrons */

              if ((ptr = (long)RDB[DATA_PTR_PREC_DET]) > VALID_PTR)
                {
                  nsrc = (long)RDB[ptr + PRECDET_N_LIVE];
                  sign = -1;
                }

              /* Get beginning time */

              t0 = TimerVal(TIMER_TRANSPORT);
              c0 = TimerCPUVal(TIMER_TRANSPORT);

              /* Start parallel timer */

              StartTimer(TIMER_OMP_PARA);

              /* Parallel loop for sampling source */

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

                for (nn = 0; nn < nsrc; nn++)
                  {
                    /* Calculate particle index */

                    idx = (long)RDB[DATA_NHIST_TOT];
                    idx = idx + (long)RDB[DATA_SRC_POP]*nb + nn;

                    /* Sample source point */

                    SampleSrcPoint(id, sign*nn, idx);
                  }
              }

              /* Parallel loop over histories */
              /* Track neutrons */

              tosimulate = 1;

              while (tosimulate > 0)
                {

                  /* Track current generation */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
                  {
                    /* Get Open MP thread id */

                    id = OMP_THREAD_NUM;

                    /* Track neutrons */

                    Tracking(id);
                  }

                  /* Even out ques for next generation */

                  tosimulate = ReDistributeQues();

                }

              /* Stop parallel timer */

              StopTimer(TIMER_OMP_PARA);

              /* Get banked precursors */

              GetBankedPrecursors();

              /* Store precursors at the end of the interval */
              /* Do we have to do something for the mesh precursors? */

              PrecursorsToStore();

              /* Add to batch counters */

              WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_COUNT] + 1.0;

              WDB[DATA_MICRO_CALC_BATCH_COUNT] =
                RDB[DATA_MICRO_CALC_BATCH_COUNT] + 1.0;

              /* Check batch interval */

              if ((long)RDB[DATA_BATCH_COUNT] == (long)RDB[DATA_BATCH_INTERVAL])
                {

                  /* Collect and clear buffered results */

                  CollectDynData();
                  CollectResults();
                  CollectDet();

                  /* Clear buffers after each batch */

                  ClearBuf();

                  /* Reset batch counter */

                  WDB[DATA_BATCH_COUNT] = 0.0;
                }

              /* Store neutrons at the end of the interval */

              BanksToStore();

              /* Clear buffered results */

              ClearBuf();

              /* Flush neutrons from bank before next batch */
              /* Should be empty already */

              FlushBank();

              /* Stop cycle-wise transport timer */

              StopTimer(TIMER_TRANSPORT_CYCLE);

              if (TimerVal(TIMER_TRANSPORT_CYCLE) > 0.0)
                t0 = TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
                  TimerVal(TIMER_TRANSPORT_CYCLE);
              else
                t0 = 0.0;

              /* CPU usage */

              ptr = (long)RDB[RES_CPU_USAGE];
              AddStat(t0, ptr, 0);

              /* Print cycle output */

              PrintCycleOutput();

              /* Sort lists */

              SortAll();

              /* Reset and restart cycle-wise transport timer */

              ResetTimer(TIMER_TRANSPORT_CYCLE);
              StartTimer(TIMER_TRANSPORT_CYCLE);

            }

          /************************************/
          /* Post-batches processing          */
          /************************************/


          /* Iterate coupled calculation */

          IterateCC();

          /* Update number of neutron histories run */
          /* For RNG-seed */

          WDB[DATA_NHIST_TOT] = RDB[DATA_NHIST_TOT] +
            RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

        }
      while (RDB[DATA_ITERATE] == (double)YES);

      /* Print FINIX output from first step */

      PrintFinix();

      /* Print results */

      if (!(1 % (long)RDB[DATA_PRINT_INTERVAL]))
        {
          MatlabOutput();
          MicroDepOutput();
          DetectorOutput();
          MeshPlotter();
          PrintCoreDistr();
          PrintHistoryOutput();
          PrintPBData();
          PrintInterfaceOutput();
          RROutput();
          /*
            GeometryPlotter(NO);
          */
          FissMtxOutput();
          MORAOutput();
          WriteICMData();
        }

      /* Moving on to remaining timesteps */

      /* Set normalization coefficient as the mean of */
      /* first time interval normalization coefficients */

      /* We would like to include only the normalization coefficients */
      /* from the last iteration (OK, since we clear stats each iter) */

      ptr = (long)RDB[RES_NORM_COEF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      WDB[DATA_NORM_COEF_N] = Mean(ptr,0);

      /* Check for remaining time intervals */

      if (maxt > 1)
        {

          /* Loop over remaining time intervals */

          for (nt = 1; nt < maxt; nt++)
            {

              /* Move neutrons from EOI stores to BOI stores */

              MoveStore();

              MoveIFC();

              /* Reset counter for solution relaxation */

              WDB[DATA_SOL_REL_ITER] = 0.0;

              WDB[DATA_SOL_REL_NTOT] = 0.0;
              WDB[DATA_SOL_REL_NCUR] = RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];
              WDB[DATA_SOL_REL_N1] = RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

              /* Set time bin index */

              WDB[DATA_DYN_TB] = (double)nt;

              /* Set new time cut-off */

              WDB[DATA_TIME_CUT_TMIN] = RDB[tme + nt];
              WDB[DATA_TIME_CUT_TMAX] = RDB[tme + nt + 1];

              /* Plot geometry for subsequent time interval */

              GeometryPlotter(NO);

              if (mpiid == 0)
                {
                  /* Iterate FINIX solvers */

                  IterateFinix();

                  /* Switch this to iteratecc (or?) */

                  if (RDB[DATA_CC_SIG_MODE] != (double)SIG_MODE_NONE)
                    {
                      /* Signal external program about moving to next timestep */

                      SignalExternal(SIGUSR2);

                      /* Read updated interfaces */

                      ptr = (long)RDB[DATA_PTR_IFC0];

                      while (ptr > VALID_PTR)
                        {
                          /* Update interface */

                          ReadInterface(ptr, YES);

                          /* Next interface */

                          ptr = NextItem(ptr);

                        }

                      /* Process updated interfaces */

                      ProcessInterface(YES);
                    }
                }

              /* Broadcast updated interfaces from 0 to other tasks */
#ifdef MPI

              /* Loop over interfaces */

              loc0 = (long)RDB[DATA_PTR_IFC0];

              while (loc0 > VALID_PTR)
                {
                  sz = (long)RDB[loc0 + IFC_MEM_SIZE];

                  /* Synchronise */

                  MPI_Barrier(my_comm);

                  /* Broadcast data from task 0 to other tasks */

                  MPITransfer(&WDB[loc0], NULL, sz, 0, MPI_METH_BC);

                  /* Synchronise */

                  MPI_Barrier(my_comm);

                  /* Next interface */

                  loc0 = NextItem(loc0);

                }

#endif

              /* Iteration loop for next timestep */

              if (RDB[DATA_RUN_CC] == (double)YES)
                WDB[DATA_ITERATE] = (double)YES;

              do
                {
                  /* Prepare coupled calculation iteration */

                  PrepareCCIter();

                  /* Reset batch counter */

                  WDB[DATA_BATCH_COUNT] = 0.0;

                  for (nb = 0; nb < maxb ; nb++)
                    {
                      /* Re-open buffer for writing */

                      WDB[DATA_BUF_REDUCED] = (double)NO;

                      /* Put cycle index */

                      WDB[DATA_CYCLE_IDX] = (double)nb;

                      /* Get precursors from store to source */
                      /* and neutrons from store to banks    */

                      ParticlesFromStore();

                      /* Do population control on precursors */

                      PrecursorPopControl();

                      /* Count live neutrons and their weight */

                      CountDynSrc();

                      /* Handle decay of precursor during interval */
                      /* and calculate weight to emit from mesh    */

                      DecayMeshPrecDet();

                      /* Resize the number of live neutrons */
                      /* This is the new population control */

                      ResizeDynSrc();

                      /* Sample delayed neutron source for time interval */
                      /* After this we should have RDB[DATA_SRC_POP]     */
                      /* neutrons */

                      SampleDelnu();

                      /* Re-open buffer for writing (was closed in sampledelnu) */

                      WDB[DATA_BUF_REDUCED] = (double)NO;

                      /* Handle decay of pointwise precursors over interval */

                      DecayPointPrecDet();

                      /* Normalize source */
                      /* Should not do population control with delayed neutrons */
                      /* as the number of neutrons is already correct */

                      NormalizeDynSrc();

                      /* Get number of particles to simulate */

                      tosimulate = ReDistributeQues();

#ifdef DNPRINT
                      fprintf(outp, "Moving to transport neutrons (printing from transportcycle.c)\n");
#endif

                      /* Start parallel timer */

                      StartTimer(TIMER_OMP_PARA);

                      /* Loop until source is empty */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
                      {
                        /* Get Open MP thread id */

                        id = OMP_THREAD_NUM;

                        /* Loop over source (first generation) */

                        while(FromSrc(id) > VALID_PTR)
                          Tracking(id);

                      }

                      /* Parallel loop over histories */
                      /* Track neutrons */

                      /* Even out ques for second generation */

                      tosimulate = ReDistributeQues();

                      while (tosimulate > 0)
                        {

                          /* Track current generation */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
                          {
                            /* Get Open MP thread id */

                            id = OMP_THREAD_NUM;

                            /* Track neutrons */

                            Tracking(id);
                          }

                          /* Even out ques for next generation */

                          tosimulate = ReDistributeQues();

                        }

                      /* Stop parallel timer */

                      StopTimer(TIMER_OMP_PARA);

                      /************************************/
                      /* Post-batch processing            */
                      /************************************/

                      /* Get banked precursors */

                      GetBankedPrecursors();

                      /* Store precursors at the end of the interval */
                      /* Do we have to do something for the mesh precursors? */

                      PrecursorsToStore();

                      /* Add to batch counters */

                      WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_COUNT] + 1.0;

                      WDB[DATA_MICRO_CALC_BATCH_COUNT] =
                        RDB[DATA_MICRO_CALC_BATCH_COUNT] + 1.0;

                      /* Check batch interval */

                      if ((long)RDB[DATA_BATCH_COUNT] == (long)RDB[DATA_BATCH_INTERVAL])
                        {

                          /* Collect and clear buffered results */

                          CollectDynData();
                          CollectResults();
                          CollectDet();
                          PoisonEq();
                          SolveOTFBurn();
                          CalcMicroGroupXS();

                          /* Clear buffers after each batch */

                          ClearBuf();

                          /* Reset batch counter */

                          WDB[DATA_BATCH_COUNT] = 0.0;
                        }

                      /* Store neutrons at the end of the interval */

                      BanksToStore();

                      /* Flush neutrons from bank before next batch */

                      FlushBank();

                      /* Stop cycle-wise transport timer */

                      StopTimer(TIMER_TRANSPORT_CYCLE);

                      if (TimerVal(TIMER_TRANSPORT_CYCLE) > 0.0)
                        t0 = TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
                          TimerVal(TIMER_TRANSPORT_CYCLE);
                      else
                        t0 = 0.0;

                      /* CPU usage */

                      ptr = (long)RDB[RES_CPU_USAGE];
                      AddStat(t0, ptr, 0);

                      /* Print cycle output */

                      PrintCycleOutput();

                      /* Sort lists */

                      SortAll();

                      /* Reset and restart cycle-wise transport timer */

                      ResetTimer(TIMER_TRANSPORT_CYCLE);
                      StartTimer(TIMER_TRANSPORT_CYCLE);
                    }

                  /* After all batches have been calculated */

                  /* Iterate coupled codes */

                  IterateCC();

                  /* Check iteration flag */

                }
              while (RDB[DATA_ITERATE] == (double)YES);

              /* Print FINIX output on each time step */

              PrintFinix();

              /* Print results */

              if (!((nt + 1) % (long)RDB[DATA_PRINT_INTERVAL]))
                {
                  MatlabOutput();
                  MicroDepOutput();
                  DetectorOutput();
                  MeshPlotter();
                  PrintCoreDistr();
                  PrintHistoryOutput();
                  PrintPBData();
                  PrintInterfaceOutput();
                  RROutput();
                  /*
                    GeometryPlotter(NO);
                  */
                  FissMtxOutput();
                  MORAOutput();
                  WriteICMData();
                }
            }
        }

      /* Signal externally coupled program to end calculation */

      SignalExternal(SIGTERM);

      /* Stop cycle-wise transport timer */

      StopTimer(TIMER_TRANSPORT_CYCLE);

      /* Collect final mesh-based precursor detector values */

      FinalizeCCPrecDet();

      /* Move EOI-stores to BOI-stores */

      MoveStore();

      /* Loop over batches to get live neutrons and precursor points */

      for (nb = 0; nb < maxb ; nb++)
        {

          /* Put cycle index */

          WDB[DATA_CYCLE_IDX] = (double)nb;

          /* Clear scoring buffers */

          ClearBuf();

          /* Get particles of this batch from store */

          ParticlesFromStore();

          /* Reduce scoring buffers */

          ReduceBuffer();

          /* Write dynamic source at the end of the simulation */

          WriteDynSrc();

          /* Flush banks and precursor source since particlesfromstore filled */
          /* them */

          FlushBank();
          FlushPrecSource();
        }

      /* Remove store files */

      RemoveStoreFiles();

      /***********************************************************************/
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /***********************************************************************/

      /***** Neutron criticality source simulation ***************************/

      /* Reset simulated batch size */

      WDB[DATA_SIMUL_BATCH_SIZE] = 0.0;

      /* Reset number of histories in gpop mode */

      WDB[DATA_GROW_POP_NHIST] = 0.0;

      /* Reset time bin index */

      WDB[DATA_DYN_TB] = 0.0;

      if ((RDB[DATA_USE_FSP] == (double)NO) ||
          ((RDB[DATA_BURN_STEP] + RDB[DATA_SOL_REL_ITER] == 0.0)
           && (RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)))
        {
          /* Check if response-matrix based acceleration is used */

          if ((long)RDB[DATA_RMX_CONVG_ACC] == YES)
            InitialCritSrc();

          fprintf(outp, "Sampling initial source...\n");

#ifdef OPEN_MP
#pragma omp parallel private(id, idx)
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

                idx = (long)RDB[DATA_NHIST_TOT] + nn;

                /* Sample source point */

                SampleSrcPoint(id, nn, idx);
              }
          }

          fprintf(outp, "OK.\n\n");

          /* Set number of inactive batches */

          skip = (long)RDB[DATA_CRIT_SKIP];
          nb0 = 0;
        }
      else
        {
          fprintf(outp, "Continuing from previous fission source\n\n");

          ResizeFissionSrc();

          /* Set number of inactive batches                          */
          /* by setting the initial inactive batch number            */
          /* This should make sure that the output and processing    */
          /* is done after the same number of live batches as before */

          nb0 = (long)RDB[DATA_CRIT_SKIP] - (long)RDB[DATA_FSP_CRIT_SKIP];

          skip = (long)RDB[DATA_CRIT_SKIP];
        }

      /* Start cycle-wise transport timer */

      ResetTimer(TIMER_TRANSPORT_CYCLE);
      StartTimer(TIMER_TRANSPORT_CYCLE);

      /* Set number of batches */

#ifdef MPI_MODE1

      maxb = (long)RDB[DATA_CRIT_CYCLES];

#else

      maxb = (long)(RDB[DATA_CRIT_CYCLES]/mpitasks);

#endif

      /* Check batching interval */

      if (maxb % (long)RDB[DATA_BATCH_INTERVAL])
        Error(0,
            "Total number of batches %ld is not a multiple of interval %ld",
            maxb, (long)RDB[DATA_BATCH_INTERVAL]);

      /* Loop over batches */

      for (nb = nb0; nb < maxb + skip; nb++)
        {
          /* Check number of skip cycles */

          if (nb == skip)
            {
              /* Clear statistics */

              ClearStat(-1);

              /* Start active transport timer */

              StartTimer(TIMER_TRANSPORT_ACTIVE);

              /* Reset number of active neutron histories */

              WDB[DATA_NHIST_CYCLE] = 0.0;

              /* Reset batch counters */

              WDB[DATA_BATCH_COUNT] = 0.0;
              WDB[DATA_MICRO_CALC_BATCH_COUNT] = 0.0;
              WDB[DATA_MICRO_CALC_BATCH_NUM] = 0.0;
            }

          /* Check reset for user iteration */

          if (((long)RDB[DATA_ITER_MODE] == ITER_MODE_USER) &&
              ((nb == skip + (long)RDB[DATA_ITER_NCYC]) ||
               (nb == skip + (long)(2.0*RDB[DATA_ITER_NCYC])) ||
               (nb == skip + (long)(3.0*RDB[DATA_ITER_NCYC]))))
            {
              /* Clear statistics */

              ClearStat(-1);

              /* Reset batch counters */

              WDB[DATA_BATCH_COUNT] = 0.0;
              WDB[DATA_MICRO_CALC_BATCH_COUNT] = 0.0;
              WDB[DATA_MICRO_CALC_BATCH_NUM] = 0.0;
            }

          /* Put cycle index */

          WDB[DATA_CYCLE_IDX] = (double)nb;

          /* Normalize source */

          NormalizeCritSrc();

          /* Clear events */

          ProcessEvents();

          /* Clear sensitivity blocks */

          ProcessSensEBlocks();

          /* Get beginning time */

          t0 = TimerVal(TIMER_TRANSPORT);
          c0 = TimerCPUVal(TIMER_TRANSPORT);

          /* MGa: Reset the DD communications  */

          ResetDDComms();

          /* Start parallel timer */

          StartTimer(TIMER_OMP_PARA);

          /* Loop until source is empty */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
          {
            /* Get Open MP thread id */

            id = OMP_THREAD_NUM;

            /* Loop over source */

            while(FromSrc(id) > VALID_PTR)
              Tracking(id);
          }

          /* Stop parallel timer */

          StopTimer(TIMER_OMP_PARA);

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

              /* Start parallel timer */

              StartTimer(TIMER_OMP_PARA);

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

              /* Stop parallel timer */

              StopTimer(TIMER_OMP_PARA);
            }

          /* Get end time */

          t0 = TimerVal(TIMER_TRANSPORT) - t0;
          c0 = TimerCPUVal(TIMER_TRANSPORT) - c0;

          /* Score time */

          ptr = (long)RDB[RES_CYCLE_RUNTIME];
          AddStat(t0, ptr, 0);
          AddStat(c0, ptr, 1);

          /* Reset normalization coefficients */

          WDB[DATA_NORM_COEF_N] = -1.0;
          WDB[DATA_NORM_COEF_G] = -1.0;

          /* K-eff iteration */

          IterateKeff();

          /* Add to batch counters */

          WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_COUNT] + 1.0;

          WDB[DATA_MICRO_CALC_BATCH_COUNT] =
            RDB[DATA_MICRO_CALC_BATCH_COUNT] + 1.0;

          /* Check batch interval */

          if ((long)RDB[DATA_BATCH_COUNT] == (long)RDB[DATA_BATCH_INTERVAL])
            {
              /* Collect and clear buffered results */

              CollectResults();
              CalcMicroGroupXS();
              CollectPrecDet();
              CollectDet();
              CollectSensResults();
              CollectUncResults();
              PoisonEq();
              SolveOTFBurn();
              ClearBuf();

              /* Reset batch counter */

              WDB[DATA_BATCH_COUNT] = 0.0;
            }

          /* Stop cycle-wise transport timer */

          StopTimer(TIMER_TRANSPORT_CYCLE);

          if (TimerVal(TIMER_TRANSPORT_CYCLE) > 0.0)
            t0 = TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
              TimerVal(TIMER_TRANSPORT_CYCLE);
          else
            t0 = 0.0;

          /* CPU usage */

          ptr = (long)RDB[RES_CPU_USAGE];
          AddStat(t0, ptr, 0);

          /* Print cycle-wise output */

          PrintCycleOutput();

          /* Reset and restart cycle-wise transport timer */

          ResetTimer(TIMER_TRANSPORT_CYCLE);
          StartTimer(TIMER_TRANSPORT_CYCLE);

          /* Sort lists */

          SortAll();

          /* Print results */

          if ((nb > skip) && (!((nb - skip + 1) %
                                (long)RDB[DATA_PRINT_INTERVAL])))
            {
              MatlabOutput();
              MicroDepOutput();
              DetectorOutput();
              SensitivityOutput();
              MeshPlotter();
              PrintCoreDistr();
              PrintHistoryOutput();
              PrintPBData();
              PrintFinix();
              PrintPrecDet();
              PrintInterfaceOutput();
              RROutput();
              /*
              GeometryPlotter(NO);
              */
              FissMtxOutput();
               MORAOutput();
              WriteICMData();
            }

          /* Add to time vector */

          ptr = (long)RDB[DATA_BTCH_TIME];
          WDB[ptr + nb] = TimerVal(TIMER_RUNTIME);

          /* Stopping criterion for gpop */

          if ((nb > skip) && ((long)RDB[DATA_GROW_POP_SIM] == YES) &&
              ((long)RDB[DATA_NHIST_CYCLE] > (long)RDB[DATA_GROW_POP_MAX_HIS]))
            break;
        }

      /* Flush bank if not passing it to next step */

      if(RDB[DATA_USE_FSP] == (double)NO)
        FlushBank();

      /* Stop cycle-wise transport timer */

      StopTimer(TIMER_TRANSPORT_CYCLE);

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid mode");

  /***************************************************************************/

  /***** Transport cycle completed *******************************************/

  /* Collect results from MPI tasks */

  CollectParallelData();

  /* Call importance map solver */

  ImportanceSolver();

  /* Put completed flag */

  WDB[DATA_SIMULATION_COMPLETED] = (double)YES;

  /* Dump buffered source points into file */

  WriteSourceFile(-1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1);

  /* Clear scoring buffer */

  ClearBuf();

  /* Stop transport timers */

  StopTimer(TIMER_TRANSPORT);
  StopTimer(TIMER_TRANSPORT_TOTAL);

  if (((long)RDB[DATA_CRIT_CYCLES] > 0) || ((long)RDB[DATA_SRC_BATCHES] > 0))
    StopTimer(TIMER_TRANSPORT_ACTIVE);

  /* Remember previous value */

  if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
    WDB[DATA_PRED_TRANSPORT_TIME] = TimerVal(TIMER_TRANSPORT);
  else
    WDB[DATA_CORR_TRANSPORT_TIME] = TimerVal(TIMER_TRANSPORT);

  /* Set estimate for coefficient calculation */

  WDB[DATA_COEF_TRANSPORT_TIME] = TimerVal(TIMER_TRANSPORT);

  /* Put keff for burnup iteration */

  if ((long)RDB[DATA_OPTI_IMPLICIT_RR] == YES)
    ptr = (long)RDB[RES_IMP_KEFF];
  else
    ptr = (long)RDB[RES_COL_KEFF];

  WDB[DATA_BURN_PREV_KEFF] = Mean(ptr, 0);
  WDB[DATA_BURN_PREV_DKEFF] = StdDev(ptr, 0);

  /* Put poison concentrations */

  PutPoisonConc();

  /* Put on-the-fly burnup concentrations */

  PutOTFBurnConc();

  /* Calculate activities (poison concentrations may be updated) */

  CalculateActivities();

  /* Print output */

  MatlabOutput();
  MicroDepOutput();
  DetectorOutput();
  SensitivityOutput();
  MeshPlotter();
  PrintCoreDistr();
  PrintHistoryOutput();
  PrintPBData();
  PrintFinix();
  PrintPrecDet();
  PrintInterfaceOutput();
  RROutput();

  /* This is called from VRCycle() in VR iteration mode */

  if ((long)RDB[DATA_RUN_VR_ITER] == NO)
    GeometryPlotter(NO);

  FissMtxOutput();
  MORAOutput();
  WriteICMData();
  /*
  ARESOutput();
  */

  /* Statistical tests */

  StatTests();

  /* Print cell statistics */

  PrintCellStats();

  /***************************************************************************/
}

/*****************************************************************************/
