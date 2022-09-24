/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : main.c                                         */
/*                                                                           */
/* Created:       2010/11/22 (JLe)                                           */
/* Last modified: 2020/06/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Main program file                                            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Main:"

/*****************************************************************************/

int main(int argc, char** argv)
{
  long ptr, idx[10000], ncoef, nhisv, more, cyc;
  char str[MAX_STR];
  double t;

  /* Reset coefficient and history variation indexes and transport time */

  ncoef = -1;
  nhisv = -1;

  cyc = 0;
  t = 0.0;

  /* Initialise MPI */

  InitMPI(argc, argv);

  /* Loop over branches */

  do
    {
      /* Start runtime and init timers */

      StartTimer(TIMER_RUNTIME);
      ResetTimer(TIMER_INIT);
      StartTimer(TIMER_INIT);
      StartTimer(TIMER_INIT_TOTAL);

      /* Init main data block */

      InitData();

      /* Initialise signal handler */

      InitSignal();

      /***********************************************************************/

      /***** Initial processing before MPI parallelization *******************/

      /* Check MPI id number */

      if (mpiid == 0)
        {
          /* Get system stat */

          SystemStat();

          /* Process command line */

          if (ParseCommandLine(argc, argv) < 0)
            return 0;

          /* Check coefficient calculation mode */

          if (((long)RDB[DATA_COEF_CALC_SPECIAL_MODE] ==
               SPECIAL_COEF_MODE_COE_ONLY) && (ncoef == -1))
            ncoef = 0;

          /* Init OpenMP related stuff */

          InitOMP();

          /* Check particle disperser mode */

          if ((long)RDB[DATA_PARTICLE_DISPERSER_MODE] == YES)
            {
              /* Run disperse subroutine */

              Disperse();

              /* Free memory */

              FreeMem();

              /* Exit */

              exit(-1);
            }

          /* Remove warning message file */

          sprintf(str, "%s.wrn", GetText(DATA_PTR_INPUT_FNAME));
          remove(str);

          /* Check for coefficient calculation */

          if ((nhisv == -1) &&
              ((ncoef == -1) || (((long)RDB[DATA_COEF_CALC_SPECIAL_MODE] ==
                                  SPECIAL_COEF_MODE_COE_ONLY) && (ncoef == 0))))
            {
              /* print title */

              PrintTitle();

              /* Time stamp */

              fprintf(outp, "Begin calculation on %s\n\n",
                      GetText(DATA_PTR_DATE));
            }

          /* Read input */

          ReadInput(GetText(DATA_PTR_INPUT_FNAME));
          fprintf(outp, "\n");

          /* Initialise socket based communication */

          InitSocket();

          /* Process surfaces */

          ProcessSurfaces();

          /* Reconfigure complement cells */

          ProcessComplementCells();

          /* Check for re-deplete */

          if ((long)RDB[DATA_PARTICLE_REDEPLETE_MODE] == YES)
            {
              /* Process inventory list */

              ProcessInventory();

              /* Print depletion output */

              PrintDepOutput();

              /* Free memory */

              FreeMem();

              /* Exit */

              exit(-1);
            }

          /* Check that the mode is right for group constant generation */

          if ((long)RDB[DATA_PTR_GCU0] > 0)
            if ((long)RDB[DATA_OPTI_MODE] != 4)
              Note(0, "Optimization mode 4 shoud be used for GC generation");

          /* Check coefficient, history variation and case matrix stuff */

          CheckCoefCalc();

          /* Coefficient calculations */

          ncoef = SetCoefCalc(ncoef);

          /* History variations calculations */

          nhisv = SetHisvCalc(nhisv, cyc);

          /* Reset burnup mode if no burnable materials are defined, or */
          /* set if restart file is read. (miks tää ei voi olla tuolla  */
          /* setoptimization.c:ssä?) */

          if ((long)RDB[DATA_BURN_MATERIALS_FLAG] == NO)
            WDB[DATA_BURNUP_CALCULATION_MODE] = (double)NO;
          else if ((long)RDB[DATA_READ_RESTART_FILE] == YES)
            WDB[DATA_BURNUP_CALCULATION_MODE] = (double)YES;

          /* Check for coefficient calculation */

          if ((long)RDB[DATA_PTR_COEF0] > VALID_PTR)
            {
              /* Check run index */

              if ((long)RDB[DATA_COEF_CALC_IDX] < 0)
                {
                  /* Original run, set restart file writing */

                  WDB[DATA_WRITE_RESTART_FILE] = (double)YES;

                  /* Tän nimen nollaaminen on ehkä aina tarpeetonta, */
                  /* mutta case matrixien kanssa se johtaa erroriin. */

                  if ((long)RDB[DATA_CASEMTX_RUN] == NO)
                    WDB[DATA_RESTART_WRITE_PTR_FNAME] = NULLPTR;

                  /* Switch group constant calculation off */

                  WDB[DATA_PTR_GCU0] = NULLPTR;
                  WDB[DATA_PTR_ADF0] = NULLPTR;
                  WDB[DATA_PTR_PPW0] = NULLPTR;
                  WDB[DATA_PTR_TRC0] = NULLPTR;
                  WDB[DATA_PTR_MDEP0] = NULLPTR;
                  WDB[DATA_OPTI_GC_CALC] = (double)NO;
                  WDB[DATA_B1_CALC] = (double)NO;
                  WDB[DATA_OPTI_POISON_CALC] = (double)NO;
                }
              else
                {
                  /* Set or reset burnup mode (set on earlier when restart */
                  /* file flag is checked, but should be off if branches   */
                  /* are run without burnup). */

                  if ((long)RDB[DATA_BURN_PTR_DEP] > VALID_PTR)
                    WDB[DATA_BURNUP_CALCULATION_MODE] = (double)YES;
                  else
                    WDB[DATA_BURNUP_CALCULATION_MODE] = (double)NO;

                  /* Reset restart file writing and depletion history */

                  WDB[DATA_WRITE_RESTART_FILE] = (double)NO;

                  /* Tän nimen nollaaminen on ehkä aina tarpeetonta, */
                  /* mutta case matrixien kanssa se johtaa erroriin. */

                  if ((long)RDB[DATA_CASEMTX_RUN] == NO)
                    WDB[DATA_RESTART_READ_PTR_FNAME] = NULLPTR;

                  WDB[DATA_BURN_PTR_DEP] = NULLPTR;
                  WDB[DATA_BURN_TOT_STEPS] = 0.0;

                  /* Switch equilibrium poison calculations off */
                  /* (NOTE: the output value will be zero) */

                  WDB[DATA_XENON_EQUILIBRIUM_MODE] = -1.0;
                  WDB[DATA_PTR_XENON_MAT_LIST] = NULLPTR;

                  WDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] = -1.0;
                  WDB[DATA_PTR_SAMARIUM_MAT_LIST] = NULLPTR;
                }
            }

          /* Check restart file (NOTE: Ei-palamamoodi toimii versioss 2.1.24 */
          /* pelkästään fotonilaskussa. Neutronidatamuotoiset ACE-nuklidit   */
          /* korvataan tässä ensin hajoamisdatalla, joka sitten myöhemmin    */
          /* muutetaan alkuainekohtaisiksi koostumuksiksi ja fotonidataksi.  */

          if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)
            ReadRestartFile(RESTART_CHECK);
          else
            ReadRestartFile(RESTART_REPLACE);

          /* Process time binnings */

          ProcessTimeBins();

          /* Remove void cells */

          RemoveVoidCells();

          /* Read STL geometries */

          ReadSTLGeometry();

          /* Check duplicate input definitions */

          CheckDuplicates();

          /* Process FINIX definitions */

          ProcessFinix();

          /* Read pebble bed geometries */

          ReadPBGeometry();

          /* Read unstructured mesh based geometry */

          ReadUMSHGeometry();

          /* Process depletion history (voidaan kutsua jo tässä) */

          ProcessDepHis();

          /* Set optimization */

          SetOptimization();

          /* Initialize secondary RNG */

          srand48(parent_seed);

          /* Update memory size */

          WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] +
            (double)MemCount();

          /* Process burnup material divisors */

          ProcessDivisors();

          /* Processing for MSR calculations */

          ProcessMSR();

          /* Divide burnable zones */

          DivideBurnMat();

          /* Update memory size */

          WDB[DATA_TOT_MAT_BYTES] = RDB[DATA_TOT_MAT_BYTES] +
            (double)MemCount();

          /* Process stochastic geometries */

          ProcessPBGeometry();

          /* Process unstructured mesh based geometries */

          ProcessUMSHGeometry();

          /* This is used for testing and debugging only (terminates run) */

          if (1 == 2)
            WriteUMSHtoSTL();

          /* Create universes in geometry */

          CreateGeometry();

          /* Process reprocessors */

          ProcessReprocessors();

          /* Count number of zones */

          ZoneCount(-1, -1, 0);

          /* Create super-imposed search meshes */

          ProcessCellMesh();

          /* Process universe transformations */

          ProcessTransformations();

          /* Process universe symmetries */

          ProcessSymmetries();

          /* Check and remove unused stuff */

          CheckUnused();

          /* Process cells */

          ProcessCells();

          /* Process nests */

          ProcessNests();

          /* Process STL geometries */

          ProcessSTLGeometry();

          /* Set material pointers */

          FindMaterialPointers();

          /* Process boundary conditions (must be called after symmetries */
          /* and transformations are processed) */

          ProcessBC();

          /* Update memory size */

          WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] +
            (double)MemCount();

          /* Make depletion zones */

          MakeDepletionZones(-1, -1, 0, 0, 0, idx, 0.0, 0.0, 0.0, 0);

          /* Update memory size */

          WDB[DATA_TOT_MAT_BYTES] = RDB[DATA_TOT_MAT_BYTES] +
            (double)MemCount();

          /* Process sources, energy grids and detectors (must be done before */
          /* unused cells and materials are removed) */

          InitPrecDetSource();
          ProcessSources();
          ProcessUserEGrids();

          /* Update memory size */

          WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] +
            (double)MemCount();

          /* Add live and file detectors for precursors */

          AllocPrecDet();

          /* Process detectors */

          ProcessDetectors();

          /* Update memory size */

          WDB[DATA_TOT_RES_BYTES] = RDB[DATA_TOT_RES_BYTES] +
            (double)MemCount();

          /* Remove unused materials */

          ptr = (long)RDB[DATA_PTR_M0];
          RemoveFlaggedItems(ptr, MATERIAL_OPTIONS, OPT_USED, NO);

          /* Process lattices */

          ProcessLattices();

          /* Find universe boundaries */

          UniverseBoundaries();

          /* Calculate nest volumes */

          NestVolumes();

          /* Calculate cell volumes */

          CellVolumes();

          /* Count number of cells */

          CellCount(-1, -1, 0, 1);

          /* Calculate material volumes */

          MaterialVolumes();

          /* MGa: Set MPI id's for decomposed materials */

          SetDDIDSimple();

          /* Print geometry data to output file */

          PrintGeometryData();

          /* Update memory size */

          WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] +
            (double)MemCount();

          /* Process core power distributions */

          ProcessCPD();

          /* Process statistics */

          ProcessStats();

          /* Process entropy stuff */

          ProcessEntropy();

          /* Process GC stuff */

          ProcessGC();

          /* Process ICM stuff */

          ProcessICM();

          /* Close surface list */

          if ((ptr = (long)RDB[DATA_PTR_S0]) > VALID_PTR)
            CloseList(ptr);

          /* Close universe list */

          if ((ptr = (long)RDB[DATA_PTR_U0]) > VALID_PTR)
            CloseList(ptr);

          /* Close cell list */

          if ((ptr = (long)RDB[DATA_PTR_C0]) > VALID_PTR)
            CloseList(ptr);

          /* Update memory size */

          WDB[DATA_TOT_RES_BYTES] = RDB[DATA_TOT_RES_BYTES] +
            (double)MemCount();

          /* Process multi-physics interfaces */

          ProcessInterface((long)NO);

          /* Process user-defined interations */

          UserIter(1);

          /* Initialize internally coupled codes */

          InitInternal();

          /* Prepare adaptive cell search meshes */

          PrepareCellSearchMesh();

          /* Pre-sort cell lists */

          PreSort();

          /* Find materials at given positions */

          MatPos();

          /* Monte Carlo volume calculator */

          VolumesMC();

          /* Experimental version of the disperser routine */

          Disperse2();

          /* Test STL geometries */

          TestSTLGeometry();

          /* Break if command-line volume MC mode */

          if ((long)RDB[DATA_VOLUME_CALCULATION_MODE] == YES)
            return -1;

          /* Process variance reduction and response matrix stuff */
          /* (pitää kutsua ennen ExpandPrivateArrays():tä */

          ProcessVR();
          ProcessRMX();

          /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel */
          /* calculation */

          ExpandPrivateArrays();

          /* Used for reverse-engineering STL geometries */

          if (1 == 2)
            STLMatFinder();

          /* Interactive command-line plotter */

          CommandLinePlotter();

          /* Plot geometry */

          GeometryPlotter(YES);

          /* Sample interface data if requested */

          SampleIFCData(YES);

          /* Update memory size */

          WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] +
            (double)MemCount();

          /* Process burnup interface */

          ProcessBUInterface();

          /* Process nuclides */

          ProcessNuclides();

          /* Update memory size */

          WDB[DATA_TOT_XS_BYTES] = RDB[DATA_TOT_XS_BYTES] +
            (double)MemCount();

          /* Process data-interfaces */

          ProcessDataInterfaces();

          /* Read data-interfaces */

          ReadDataInterfaces();

          /* Update memory size */

          WDB[DATA_TOT_IFC_BYTES] = RDB[DATA_TOT_IFC_BYTES] +
            (double)MemCount();

          /* Process inventory list */

          ProcessInventory();

          /* Update memory size */

          WDB[DATA_TOT_RES_BYTES] = RDB[DATA_TOT_RES_BYTES] +
            (double)MemCount();
        }

      /***********************************************************************/

      /**** MPI parallel part ************************************************/

      /* Distribute data to parallel MPI tasks */

      ShareInputData();

      /* Update memory size */

      WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] + (double)MemCount();

      /* Process energy grids */

      UnionizeGrid();

      /* Process XS data */

      ProcessXSData();

      /* Generate cache-optimized block */

      CacheXS();

      /* Update memory size */

      WDB[DATA_TOT_XS_BYTES] = RDB[DATA_TOT_XS_BYTES] + (double)MemCount();

      /* Process data for sensitivity calculation */

      ProcessSensitivities();

      /* Process the covariance data that was read in processnuclides */

      ProcessCovarianceData();

      /* Process statistics for sensitivity calculation */

      ProcessSensStats();

      /* Update memory size */

      WDB[DATA_TOT_RES_BYTES] = RDB[DATA_TOT_RES_BYTES] + (double)MemCount();

      /* Allocate memory for precursor statistics. Loops over nuclide list */
      /* Needs precursor group lists to be set (set at ProcessXSData) */

      ProcessPrecDet();

      /* Process mesh plots */

      ProcessMeshPlots();

      /* Update memory size */

      WDB[DATA_TOT_RES_BYTES] = RDB[DATA_TOT_RES_BYTES] + (double)MemCount();

      /* Update memory size */

      WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] + (double)MemCount();

      /* Process materials */

      ProcessMaterials();

      /* Update memory size */

      WDB[DATA_TOT_MAT_BYTES] = RDB[DATA_TOT_MAT_BYTES] + (double)MemCount();

      /* This is used for testing and debugging only (terminates run) */

      if (1 == 2)
        WriteTetMeshtoGeo();

      /* Link reactions to sources and detectors */

      LinkReactions();

      /* Allocate memory for interface statistics */

      AllocInterfaceStat();

      /* Update memory size */

      WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] + (double)MemCount();

      /* Process fission matrixes */

      ProcessFissMtx();

      /* Update memory size */

      WDB[DATA_TOT_RES_BYTES] = RDB[DATA_TOT_RES_BYTES] + (double)MemCount();

      /* Init particle structures */

      InitHistories();

      /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */

      ExpandPrivateArrays();

     /* MGa: Initialize the DD communications */

      InitDDComms();

      /* Disallow memory allocation from here on*/

      Mem(MEM_DENY);

      /* Check if calculation should proceed */

      if (((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES) ||
          ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES) ||
          ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES))
        {
          /* Stop init timer */

          StopTimer(TIMER_INIT);
          StopTimer(TIMER_INIT_TOTAL);

          /* Set coefficient transport time */

          WDB[DATA_COEF_TRANSPORT_TIME] = t;

          /* Check mode and start calculation */

          if ((long)RDB[DATA_RUN_VR_ITER] == YES)
            {
              /* Variance reduction iterations */

              VRCycle();
            }
          else if ((long)RDB[DATA_COM_SOCKET] > 0)
            {
              /* Supervised calculation */

              SupervisedCycle();

              /* Deinitialize the socket communicator */

              DeinitSocket();
            }
          else if ((long)RDB[DATA_COEF_CALC_IDX] > 0)
            {
              /* Coefficient calculations */

              CoefCycle();
            }
          else if ((long)RDB[DATA_PTR_RIA0] > VALID_PTR)
            {
              /* RIA simulation */

              RIACycle();
            }
          else if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES) &&
                   ((long)RDB[DATA_BURN_PTR_DEP] > VALID_PTR))
            {
              /* Internal burnup mode, go to burnup cycle */

              BurnupCycle();
            }
          else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
            {
              /* Dynamic calculation mode, iteration is done inside */
              /* transportcycle */

              PrepareTransportCycle();

              /* Run transport cycle */

              TransportCycle();
            }
          else
            {
              /* Single transport cycle */

              PrepareTransportCycle();

              /* Run transportcycle(s) */

              do
                {
                  /* Prepare coupled calculation iteration if needed */

                  PrepareCCIter();

                  /* Run transport cycle */

                  TransportCycle();

                  /* Iterate coupled calculation routines */

                  IterateCC();

                  /* Repeat if needed */
                }
              while(RDB[DATA_ITERATE] == (double)YES);

              /* Signal externally coupled program to end calculation */

              SignalExternal(SIGTERM);
            }

          /* Remember transport time (needed for estimating total in */
          /* coefficent calculation */

          t = RDB[DATA_COEF_TRANSPORT_TIME];
        }

      /* Remember burnup cycle index */

      cyc = (long)RDB[DATA_BURN_STEP];

      /* Get flag for more coefficient calculations (must be done */
      /* before the memory is free'd) */

      if ((long)RDB[DATA_COEF_CALC_SPECIAL_MODE] == SPECIAL_COEF_MODE_HIS_ONLY)
        more = NO;
      else
        more = (long)RDB[DATA_MORE_COEF_CALC];

      /* MGa: Free the DD communications */

      FreeDDComms();

      /* Free memory */

      FreeMem();

      /* Stop runtime timer */

      StopTimer(TIMER_RUNTIME);

      /* Update coefficient index */

      ncoef++;
    }
  while (more == YES);

  /* Finalize MPI */

  FinalizeMPI();

  /* Exit */

  return 0;
}

/*****************************************************************************/
