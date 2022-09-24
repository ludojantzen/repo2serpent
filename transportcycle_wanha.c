/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : transportcycle.c                               */
/*                                                                           */
/* Created:       2011/05/23 (JLe)                                           */
/* Last modified: 2015/09/15 (JLe)                                           */
/* Version:       2.1.25                                                     */
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
  long nb, nb0, maxb, skip, nn, nt, maxt, id, idx, ptr, tme;
  double t0, c0;

  /***************************************************************************/

  /***** Main transport cycle ************************************************/

  /* Reset total and active timer */

  ResetTimer(TIMER_TRANSPORT);
  ResetTimer(TIMER_TRANSPORT_ACTIVE);

  /* Start transport timer */
      
  StartTimer(TIMER_TRANSPORT);
  StartTimer(TIMER_TRANSPORT_TOTAL);

  /* Reset completed flag and sort counter */

  WDB[DATA_SIMULATION_COMPLETED] = (double)NO;
  WDB[DATA_SORT_COUNT] = 1.0;

  /* Check mode */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC)
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

      /* Check batching interval */

      if ((long)RDB[DATA_SRC_BATCHES] % (long)RDB[DATA_BATCH_INTERVAL])
	Error(0, 
	    "Total number of batches %ld is not a multiple of interval %ld",
	    (long)RDB[DATA_SRC_BATCHES], (long)RDB[DATA_BATCH_INTERVAL]);

      /* Start simulation */

      fprintf(out, "Starting external source simulation...\n\n");

      /* Loop over batches */
      
      for (nb = 0; nb < (long)RDB[DATA_SRC_BATCHES]; nb++)
	{
	  /*******************************************************************/
	  
	  /***** Loop over source batches ************************************/

	  /* Reset time bin index */

	  WDB[DATA_DYN_TB] = 0.0;

	  /* Set time cut-off for first cycle */
	  
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
		idx = idx + (long)(nb*RDB[DATA_SRC_POP]) + nn;

		/* Sample source point */

		SampleSrcPoint(id, nn, idx);

		/* Track */

		Tracking(id);
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

		  /* Set new time cut-off */

		  WDB[DATA_TIME_CUT_TMIN] = RDB[tme + nt];
		  WDB[DATA_TIME_CUT_TMAX] = RDB[tme + nt + 1];
		  
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
		    
		    while(FromSrc(id) > VALID_PTR)
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

	  /******************************************************************/
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
	       
      /* Start active transport timer */

      StartTimer(TIMER_TRANSPORT_ACTIVE);

      /* Start cycle-wise transport timer */

      ResetTimer(TIMER_TRANSPORT_CYCLE);
      StartTimer(TIMER_TRANSPORT_CYCLE);

      /* Start simulation */

      fprintf(out, "Starting time dependent simulation...\n\n");
      
      /* Reset time bin index */

      WDB[DATA_DYN_TB] = 0.0;

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

      do
	{
	  /* Prepare coupled calculation iteration */

	  ClearInterfaceStat();
	  
	  /* Re-open buffer for writing */
	  
	  WDB[DATA_BUF_REDUCED] = (double)NO;
	  
	  /* Loop over batches */
	  
	  for (nb = 0; nb < (long)RDB[DATA_SRC_BATCHES]; nb++)
	    {
	      
	      /* Re-open buffer for writing */
	      
	      WDB[DATA_BUF_REDUCED] = (double)NO;
	      
	      /* Set batch counters */
	      
	      WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_INTERVAL];
	      
	      WDB[DATA_MICRO_CALC_BATCH_COUNT] = 
		RDB[DATA_MICRO_CALC_BATCH_SIZE];
	      
	      /* Put cycle index */
	      
	      WDB[DATA_CYCLE_IDX] = (double)nb;
	      
	      /* Print cycle output */
	      
	      PrintCycleOutput();
	    
	      /* Get beginning time */
	      
	      t0 = TimerVal(TIMER_TRANSPORT);
	      c0 = TimerCPUVal(TIMER_TRANSPORT);
	      
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
		    idx = idx + (long)(nb*RDB[DATA_SRC_POP]) + nn;
		    
		    /* Sample source point */
		    
		    SampleSrcPoint(id, nn, idx);
		    
		    /* Track */
		    
		    Tracking(id);
		  }
	      }
	      
	      /* Store neutrons at the end of the interval */
	      
	      BanksToStore();
	      
	      /* Flush neutrons from bank before next batch */ 
	      
	      FlushBank();
	      
	      /* Collect detectors */
	      
	      CollectDet();
	      
	      /* Stop parallel timer */
	      
	      StopTimer(TIMER_OMP_PARA);	    
	      
	      /* Collect data from interval */
	      
	      CollectDynData();
	      
	      /* Reset normalization coefficients */
	      /* Separate coefficients for separate batches */
	      /* Will be averaged for subsequent time intervals */
	      
	      WDB[DATA_NORM_COEF_N] = -1.0;
	      WDB[DATA_NORM_COEF_G] = -1.0;
	      
	      NormCoef(PARTICLE_TYPE_GAMMA);
	      
	      /* Collect results */
	      
	      CollectResults();
	      
	      /* Clear buffers after each batch */
	      
	      ClearBuf();
	      
	    }
	  
	  /************************************/
	  /* Post-batch processing            */
	  /************************************/
	  
	  /* Iterate coupled calculation */
	  
	  IterateCC();
	  
	  /* Check iteration flag */
	  
	}
      while(RDB[DATA_ITERATE] == (double)YES);

      /* Moving on to remaining timesteps */

      /* Set normalization coefficient as the mean of */
      /* first time interval normalization coefficients */
      
      ptr = (long)RDB[RES_NORM_COEF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      WDB[DATA_NORM_COEF_N] = Mean(ptr,0);
      
      /* Check for remaining time intervals */

      if (maxt > 1)
	{
	  /* Re-open buffer for writing */

	  WDB[DATA_BUF_REDUCED] = (double)NO;
	      
	  /* Loop over remaining time intervals */
	      
	  for (nt = 1; nt < maxt; nt++)
	    {
	      /* Move neutrons from EOI stores to BOI stores */

	      MoveStore();

	      /* Reset counter for solution relaxation */

	      WDB[DATA_SOL_REL_ITER] = 0.0;

	      WDB[DATA_SOL_REL_NTOT] = 0.0;
	      WDB[DATA_SOL_REL_NCUR] = RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];
	      WDB[DATA_SOL_REL_N1] = RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

	      if (RDB[DATA_CC_SIG_MODE] != (double)SIG_MODE_NONE)
		{
		  /* Signal external program about moving to next timestep */

		  SignalExternal(SIGUSR2);
	      
		  /* Read updated interfaces */

		  ptr = (long)RDB[DATA_PTR_IFC0];

		  while(ptr > VALID_PTR)
		    {
		      /* Update interface */

		      ReadInterface(ptr, YES);

		      /* Next interface */

		      ptr = NextItem(ptr);

		    }

		  /* Process updated interfaces */

		  ProcessInterface(YES);
		}

	      /* Iteration loop for next timestep */
	      
	      do
		{
		  /* Prepare coupled calculation iteration */

		  ClearInterfaceStat();
		  
		  /* Re-open buffer for writing */
		  
		  WDB[DATA_BUF_REDUCED] = (double)NO;
		  
		  /* Set time bin index */
		  
		  WDB[DATA_DYN_TB] = (double)nt;
		  
		  /* Set new time cut-off */
		  
		  WDB[DATA_TIME_CUT_TMIN] = RDB[tme + nt];
		  WDB[DATA_TIME_CUT_TMAX] = RDB[tme + nt + 1];
		  
		  for (nb = 0; nb < (long)RDB[DATA_SRC_BATCHES]; nb++)
		    {
		      /* Re-open buffer for writing */
		      
		      WDB[DATA_BUF_REDUCED] = (double)NO;
		      
		      /* Set batch counter */
		      
		      WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_INTERVAL];
		      
		      /* Put cycle index */
		      
		      WDB[DATA_CYCLE_IDX] = (double)nb;
		      
		      /* Print cycle output */

		      PrintCycleOutput();

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
			
			while(FromSrc(id) > VALID_PTR)
			  Tracking(id);      
		      }
		      
		      /* Stop parallel timer */
		      
		      StopTimer(TIMER_OMP_PARA);
		      
		      /************************************/
		      /* Post-batch processing            */
		      /************************************/		
		      
		      BanksToStore();
		      
		      FlushBank();
		      
		      /* Might break something */
		      
		      CollectResults();
		      
		      /* Collect detectors */
		      
		      CollectDet();
		      
		      /* Collect data from interval */
		      
		      CollectDynData();
		      
		      /* Clear buffered results */
		      
		      ClearBuf();
		    }

		  IterateCC();
		  
		  /* Check iteration flag */
		  
		}
	      while(RDB[DATA_ITERATE] == (double)YES);

	      if((nt+1)%50 == 0)
		PrintFinix();
	    }
	}
    
      /* Signal external program about ending the calculation */
      
      SignalExternal(SIGTERM);
      
      MatlabOutput();
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

      return;

      /* Katso noi allaolevat vielä, mitä niistä pitää tehdä */
	      
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
      PoisonEq();
      CalcMicroGroupXS();
      ClearBuf();

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
       

      /* Update number of neutron histories run */

      WDB[DATA_NHIST_TOT] = RDB[DATA_NHIST_TOT] +
	RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

      /* Stop cycle-wise transport timer */

      StopTimer(TIMER_TRANSPORT_CYCLE);

      /***********************************************************************/
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /***********************************************************************/
      
      /***** Neutron criticality source simulation ***************************/

      /* Reset simulated batch size */

      WDB[DATA_SIMUL_BATCH_SIZE] = 0.0;

      /* Reset time bin index */

      WDB[DATA_DYN_TB] = 0.0;

      /* Generate initial source */

      if((RDB[DATA_USE_FSP] == (double)NO) || 
	 ((RDB[DATA_BURN_STEP] + RDB[DATA_SOL_REL_ITER] == 0.0)
         && (RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)))
	{

	  fprintf(out, "Sampling initial source...\n");

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

	  fprintf(out, "OK.\n\n");

          /* Set number of inactive batches */

	  skip = (long)RDB[DATA_CRIT_SKIP];
	  nb0 = 0;
	}
      else
	{
	  fprintf(out, "Continuing from previous fission source\n\n");	

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
	    }

	  /* Put cycle index */

	  WDB[DATA_CYCLE_IDX] = (double)nb;

	  /* Normalize source */

	  NormalizeCritSrc();

	  /* Clear events */

	  ProcessEvents();

	  /* Get beginning time */

	  t0 = TimerVal(TIMER_TRANSPORT);
	  c0 = TimerCPUVal(TIMER_TRANSPORT);

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
	      CollectDet();
	      PoisonEq();
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
	  
	  if (!((nb + 1) % (long)RDB[DATA_PRINT_INTERVAL]))
	    {
	      MatlabOutput();
	      DetectorOutput();
	      MeshPlotter();
	      PrintCoreDistr();
	      PrintHistoryOutput();
	      PrintPBData();
	      PrintFinix();
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

  /* Calculate activities (poison concentrations may be updated) */

  CalculateActivities();

  /* Store simulator output data */

  StoreSimData();

  /* Print output */

  MatlabOutput();
  DetectorOutput();
  MeshPlotter();
  PrintCoreDistr();
  PrintHistoryOutput();
  PrintPBData();
  PrintFinix();
  PrintInterfaceOutput();
  RROutput();

  GeometryPlotter(NO);

  FissMtxOutput();
  MORAOutput();
  WriteICMData();
  CAXOutput();
  /*
  ARESOutput();
  */

  /* Statistical tests */

  StatTests();

  /***************************************************************************/
}

/*****************************************************************************/
