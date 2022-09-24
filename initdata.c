/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : initdata.c                                     */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2020/05/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Inits values in main data block                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitData:"

/*****************************************************************************/

void InitData()
{
  long ptr;
  double val;
  char *path, *strnomp, *seed, tmpstr[MAX_STR];
  struct timeb tb;
  clock_t cput0;

  /* Set error pointer */

  errp = stdout;

  /* Check that system is 64 bit */

  if (sizeof(long) != 8)
    {
      /* Print error */

      fprintf(errp, "\nSerpent 2 must be run in a 64-bit system. ");

      if (sizeof(long) == 4)
        fprintf(errp, "Your system appears to be 32-bit.");

      fprintf(errp, "\n\n");

      /* Abort */

      exit(-1);
    }

  /* Reset pointers */

  RDB = NULL;
  WDB = NULL;
  ACE = NULL;
  ASCII = NULL;
  PRIVA = NULL;
  BUF = NULL;
  RES1 = NULL;
  RES2 = NULL;
  RES3 = NULL;
  SEED = NULL;
  SEED0 = NULL;

  /* Allocate memory for random number seed vectors */

  SEED = (unsigned long *)Mem(MEM_ALLOC, MAX_OMP_THREADS*RNG_SZ,
                              sizeof(unsigned long));
  SEED0 = (unsigned long *)Mem(MEM_ALLOC, MAX_OMP_THREADS*RNG_SZ,
                              sizeof(unsigned long));

  /* Set pointer to standard output */

  if (mpiid == 0)
    outp = stdout;
  else
    {
      outp = fopen("/dev/null", "w");
      return;
    }

  /* Set line-buffering for stdout */

  setlinebuf(outp);

  /* Allocate memory for main data block */

  ReallocMem(DATA_ARRAY, DATA_FIXED_BLOCK_SIZE);

  /* Allow memory allocation */

  Mem(MEM_ALLOW);

  /* Allocate data at the beginning of array to simplify error testing */

  ReallocMem(RES1_ARRAY, DATA_FIXED_BLOCK_SIZE);

  /* Allocate memory for ASCII data block to simplify error testing */

  ASCII = (char *)Mem(MEM_ALLOC, (long)(VALID_PTR + 1), sizeof(char));
  WDB[DATA_ASCII_DATA_SIZE] = (double)(VALID_PTR + 1);

  /* Reset status of private memory blocks */

  WDB[DATA_PRIVA_MEM_READY] = (double)NO;

  /* Init list pointers */

  WDB[DATA_PTR_C0] = NULLPTR;
  WDB[DATA_PTR_S0] = NULLPTR;
  WDB[DATA_PTR_M0] = NULLPTR;
  WDB[DATA_PTR_L0] = NULLPTR;
  WDB[DATA_PTR_T0] = NULLPTR;
  WDB[DATA_PTR_NST0] = NULLPTR;
  WDB[DATA_PTR_GPL0] = NULLPTR;
  WDB[DATA_PTR_TR0] = NULLPTR;
  WDB[DATA_PTR_PB0] = NULLPTR;
  WDB[DATA_PTR_IFC0] = NULLPTR;
  WDB[DATA_PTR_DIV0] = NULLPTR;

  /* Init file names */

  WDB[DATA_PTR_ACEDATA_FNAME_LIST] = NULLPTR;
  WDB[DATA_PTR_DECDATA_FNAME_LIST] = NULLPTR;
  WDB[DATA_PTR_NFYDATA_FNAME_LIST] = NULLPTR;
  WDB[DATA_PTR_SFYDATA_FNAME_LIST] = NULLPTR;

 /* Maximum fraction of system memory to use */
  /* Try to get fraction from environmental variable */

  if ((path = getenv("SERPENT_MEM_FRAC")) != NULL)
    if (sscanf(path, "%lf", &val) == 1)
      if ((val > 0.0) && (val < 1.0))
        WDB[DATA_CPU_MEM_FRAC] = val;

  /* Set to default value if not set in environmental variable */

  if (RDB[DATA_CPU_MEM_FRAC] == 0.0)
    WDB[DATA_CPU_MEM_FRAC] = 0.8;

  /* Few-group constant generation */

  WDB[DATA_ERG_FG_PTR_PREDEF] = (double)PutText("default2");

  /* Thinning tolerance and energy boundaries (tarkista noi fotonirajat) */

  WDB[DATA_ERG_TOL] = -1.0;

  WDB[DATA_NEUTRON_EMIN] = 1E-11;
  WDB[DATA_NEUTRON_EMAX] = 20.0;

  WDB[DATA_PHOTON_EMIN] = 1E-3;
  WDB[DATA_PHOTON_EMAX] = 100.0;

  /* Set boundary conditions and albedos (NOTE: Toi STOP_AT_BOUNDARY */
  /* käännettiin pois juuri ennen 2.1.25 jakelua) */

  WDB[DATA_STOP_AT_BOUNDARY] = (double)NO;
  WDB[DATA_GEOM_BC0] = (double)BC_BLACK;
  WDB[DATA_GEOM_BC1] = (double)BC_BLACK;
  WDB[DATA_GEOM_BC2] = (double)BC_BLACK;
  WDB[DATA_GEOM_BC3] = (double)BC_BLACK;
  WDB[DATA_GEOM_ALBEDO1] = 1.0;
  WDB[DATA_GEOM_ALBEDO2] = 1.0;
  WDB[DATA_GEOM_ALBEDO3] = 1.0;

  /* K-eff iteration */

  WDB[DATA_ITER_MODE] = (double)ITER_MODE_NONE;
  WDB[DATA_ITER_KEFF] = 1.0;
  WDB[DATA_ITER_VAL] = -1.0;
  WDB[DATA_ITER_NCYC] = 50.0;
  WDB[DATA_ITER_FIX] = (double)NO;

  WDB[DATA_ITER_ALB_F1] = 1.0;
  WDB[DATA_ITER_ALB_F2] = 1.0;
  WDB[DATA_ITER_ALB_F3] = 1.0;

  /* Set initial date */

  WDB[DATA_PTR_DATE] = (double)PutText(TimeStamp());

  /* Get initial cpu time */

  cput0 = clock();
  WDB[DATA_CPU_T0] = (double)cput0;

  /* Reset calculation modes */

  WDB[DATA_NEUTRON_TRANSPORT_MODE] = (double)NO;
  WDB[DATA_PHOTON_TRANSPORT_MODE] = (double)NO;
  WDB[DATA_BURNUP_CALCULATION_MODE] = (double)NO;
  WDB[DATA_VOLUME_CALCULATION_MODE] = (double)NO;
  WDB[DATA_PARTICLE_DISPERSER_MODE] = (double)NO;
  WDB[DATA_PARTICLE_REDEPLETE_MODE] = (double)NO;
  WDB[DATA_QUICK_PLOT_MODE] = (double)NO;
  WDB[DATA_INTERACTIVE_PLOT_MODE] = (double)NO;
  WDB[DATA_COPY_INPUTS] = (double)NO;
  WDB[DATA_PHOTON_PRODUCTION] = (double)NO;
  WDB[DATA_NGAMMA_SRC_SIM] = (double)NO;
  WDB[DATA_MULTI_PARTICLE_TRANSPORT] = (double)NO;
  WDB[DATA_MATPOS_N_PTS] = 0.0;

  /* Implicit photon production */

  WDB[DATA_PHOTON_IMPL_WMIN] = 0.1;
  WDB[DATA_PHOTON_IMPL_NMAX] = 10.0;

  /* Normalisation */

  WDB[DATA_PTR_NORM] = NULLPTR;
  WDB[DATA_NORM_U235_FISSE] = U235_FISSE;

  /* Replay mode */

  WDB[DATA_OPTI_REPLAY] = (double)NO;

  /* Set default random number seed */

  ftime(&tb);
  parent_seed = (unsigned long)(tb.time*1000 + tb.millitm);

  /* Override with environment variable */

  if ((seed = getenv("SERPENT_RNG_SEED")) != NULL)
    {
      parent_seed = (unsigned long)atoi(seed);
      WDB[DATA_OPTI_REPLAY] = (double)YES;
    }

  /* XS data plotter parameters */

  WDB[DATA_XSPLOT_NE] = -1.0;
  WDB[DATA_XSPLOT_EMIN] = -1.0;
  WDB[DATA_XSPLOT_EMAX] = -1.0;

  /* Print compositions */

  WDB[DATA_BURN_PRINT_COMP] = (double)NO;
  WDB[DATA_BURN_PRINT_COMP_LIM] = 1.0;
  WDB[DATA_BURN_MAT_OUTPUT] = (double)BURN_OUT_MAT_PARENT;
  WDB[DATA_BURN_PRINT_STEP] = (double)BURN_OUT_PRINT_ALL;
  WDB[DATA_BURN_PRINT_OUTPUT] = (double)YES;

  /* Decay only mode */

  WDB[DATA_BURN_DECAY_CALC] = (double)NO;

  /* Print histories */

  WDB[DATA_OPTI_PRINT_HIS] = (double)NO;

  /* Depletion cut-offs */

  WDB[DATA_DEP_HALF_LIFE_CUTOFF] = INFTY;
  WDB[DATA_DEP_TTA_CUTOFF]       = 1E-28;

  /* Time cut-off for transport simulation */

  WDB[DATA_TIME_CUT_TMIN] = 0.0;
  WDB[DATA_TIME_CUT_TMAX] = INFTY;

  /* Dynamic simulation mode with intervals */

  WDB[DATA_DYN_DT] = -INFTY;
  WDB[DATA_DYN_WMIN] = -INFTY;

  /* Reset generation cut-off */

  WDB[DATA_GEN_CUT] = (double)MAX_GENERATIONS;
  WDB[DATA_MAX_PROMPT_CHAIN_LENGTH] = -1.0;

  /* Fission source entropy */

  WDB[DATA_OPTI_ENTROPY_CALC] = (double)NO;

  WDB[DATA_ENTROPY_NX] = 5.0;
  WDB[DATA_ENTROPY_NY] = 5.0;
  WDB[DATA_ENTROPY_NZ] = 5.0;

  WDB[DATA_ENTROPY_XMIN] = -INFTY;
  WDB[DATA_ENTROPY_XMAX] =  INFTY;

  WDB[DATA_ENTROPY_YMIN] = -INFTY;
  WDB[DATA_ENTROPY_YMAX] =  INFTY;

  WDB[DATA_ENTROPY_ZMIN] = -INFTY;
  WDB[DATA_ENTROPY_ZMAX] =  INFTY;

  /* Source point animation */

  WDB[DATA_SOURCE_PT_ANIM] = (double)NO;
  WDB[DATA_SOURCE_PT_ANIM_F] = 1.0;
  WDB[DATA_SOURCE_PT_ANIM_PALETTE] = (double)PALETTE_HOT;

  /* Microscopic partial total cross section limit for including nuclide */
  /* in reaction lists. */

  WDB[DATA_MIN_TOTXS] = 1E-5;

  /* Use unresolved resonance probability tables (turned off by default) */

  WDB[DATA_USE_URES] = (double)NO;
  WDB[DATA_URES_PTR_USE_LIST] = NULLPTR;

  /* Plotter */

  WDB[DATA_STOP_AFTER_PLOT] = (double)STOP_AFTER_PLOT_NONE;
  WDB[DATA_PLOTTER_MODE] = (double)NO;
  WDB[DATA_IGNORE_GEOM_PLOTS] = (double)NO;

  /* Printing mixtures */

  WDB[DATA_DECOMPOSE_MIXTURES] = (double)NO;

  /* Optimization and memory/data options */

  WDB[DATA_OPTI_MODE] = 4.0;
  WDB[DATA_OPTI_UNIONIZE_GRID] = -1.0;
  WDB[DATA_OPTI_RECONSTRUCT_MICROXS] = -1.0;
  WDB[DATA_OPTI_RECONSTRUCT_MACROXS] = -1.0;
  WDB[DATA_OPTI_INCLUDE_SPECIALS] = (double)NO;
  WDB[DATA_OPTI_DIX] = -1.0;

  /* Delta-tracking */

  WDB[DATA_OPT_USE_DT] = (double)YES;
  WDB[DATA_DT_NTHRESH] = 0.1;
  WDB[DATA_DT_PTHRESH] = 0.1;

  /* Maximum number of cells to create search lists and presort options */

  WDB[DATA_MAX_CELL_SEARCH_LIST] = 2000.0;
  WDB[DATA_PRESORT_NP] = -1.0;
  WDB[DATA_PRESORT_NB] = 10.0;

  /* Implicit Monte Carlo */

  WDB[DATA_OPT_IMPL_FISS] = -1.0;
  WDB[DATA_OPT_IMPL_CAPT] = -1.0;
  WDB[DATA_OPT_IMPL_NXN] = -1.0;

  /* Set roulette off by default */

  WDB[DATA_OPT_ROULETTE_W0] = -1.0;
  WDB[DATA_OPT_ROULETTE_P0] = 0.5;

  /* Init cycle-keff */

  WDB[DATA_CYCLE_KEFF] = 1.0;
  WDB[DATA_N_POP_EIG] = 1.0;

  /* Simulation mode and batching interval */

  WDB[DATA_SIMULATION_MODE] = -1.0;
  WDB[DATA_BATCH_INTERVAL] = 1.0;

  /* OpenMP stuff */

  if ((strnomp = getenv("SERPENT_OMP_NUM_THREADS")) != NULL)
    WDB[DATA_OMP_MAX_THREADS] = (double)atoi(strnomp);
  else if ((strnomp = getenv("OMP_NUM_THREADS")) != NULL)
    WDB[DATA_OMP_MAX_THREADS] = (double)atoi(strnomp);
  else
    WDB[DATA_OMP_MAX_THREADS] = 1.0;

  /* Allocate empty memory ACE array to deal with zero pointers */

  WDB[DATA_PTR_ACE0] = NULLPTR;
  ReallocMem(ACE_ARRAY, DATA_FIXED_BLOCK_SIZE);

  WDB[DATA_PTR_ACE_NFY_DATA] = NULLPTR;
  WDB[DATA_PTR_ACE_SFY_DATA] = NULLPTR;

  /* Root universe */

  WDB[DATA_PTR_ROOT_UNIVERSE] = (double)PutText("0");

  /* Event bank size and maximum number of generations */

  WDB[DATA_EVENT_BANK_SZ] = 1000.0;
  WDB[DATA_EVENT_MAX_GEN] = 20.0;
  WDB[DATA_EVENT_RECORD_FLAGS] = 0.0;

#ifndef OLD_IFP

  SetOption(DATA_EVENT_RECORD_FLAGS, RECORD_EVENT_IFP);

#endif

  /* Monte Carlo volume calculation */

  WDB[DATA_VOLUME_MC] = NO;
  WDB[DATA_VOLUME_MC_NMAX] = 0.0;
  WDB[DATA_VOLUME_MC_TMAX] = -1.0;
  WDB[DATA_VOLUME_MC_EMAX] = -1.0;

  /* Default confidence levels for TMP majorant generation */

#ifndef TRADMAJO

  /* Revisited majorant */

  WDB[DATA_QPARAM_TMS] = 2.0E-5;
  WDB[DATA_QPARAM_DBRC] = 2.0E-5;

#else

  /* "f-parameter" of traditional majorant */

  WDB[DATA_QPARAM_TMS] = 3.0;
  WDB[DATA_QPARAM_DBRC] = 3.0;

#endif

  /* DBRC */

  WDB[DATA_USE_DBRC] = (double)NO;
  WDB[DATA_PTR_DBRC] = NULLPTR;
  WDB[DATA_DBRC_EMIN] = 0.4E-6;
  WDB[DATA_DBRC_EMAX] = 210E-6;

  /* Energy cut-off for transmutation reactions */

  WDB[DATA_BURN_ENECUT] = 10.0;

  /* Normalization */

  WDB[DATA_NORM_BURN] = BURN_NORM_ALL;

  /* Burnup mode */

  WDB[DATA_BURN_BUMODE] = (double)BUMODE_CRAM;
  WDB[DATA_BU_SPECTRUM_COLLAPSE] = -1.0;
  WDB[DATA_BURN_CALC_NSF] = (double)NO;

  /* CRAM:n asteluku */

  WDB[DATA_BURN_CRAM_K] = 14.0;

  /* Predictor-corrector calculation */

  WDB[DATA_BURN_PRED_TYPE] = (double)PRED_TYPE_CONSTANT;
  WDB[DATA_BURN_PRED_NSS] = 1.0;
  WDB[DATA_BURN_CORR_TYPE] = (double)CORR_TYPE_LINEAR;
  WDB[DATA_BURN_CORR_NSS] = 1.0;
  WDB[DATA_BURN_SIE] = (double)NO;

  /* Corrector iteration */

  WDB[DATA_BURN_CI_TYPE] = (double)CI_TYPE_OUTER;
  WDB[DATA_BURN_CI_MAXI] = 1.0;
  WDB[DATA_BURN_CI_LAST] = (double)NO;

  WDB[DATA_BURN_CI_NBATCH] = -1.0;
  WDB[DATA_BURN_CI_CYCLES] = -1.0;
  WDB[DATA_BURN_CI_SKIP]   = -1.0;

  /* Core power distribution */

  WDB[DATA_CORE_PDE_DEPTH] = -1.0;
  WDB[DATA_CORE_PDE_L1] = -1.0;
  WDB[DATA_CORE_PDE_L2] = -1.0;

  /* Uniform fission source method */

  WDB[DATA_UFS_PTR_SRC_MESH] = NULLPTR;
  WDB[DATA_UFS_ORDER] = 1.0;
  WDB[DATA_UFS_MIN] = 0.001;
  WDB[DATA_UFS_MAX] = 10.0;

  /* Energy grid for coarse multi-group xs */

  WDB[DATA_COARSE_MG_NE] = 4000.0;
  WDB[DATA_OPTI_MG_MODE] = -1.0;

  /* Reset minimum and maximum energies of cross section data */

  WDB[DATA_NEUTRON_XS_EMIN] = INFTY;
  WDB[DATA_NEUTRON_XS_EMAX] = -INFTY;

  WDB[DATA_PHOTON_XS_EMIN] = INFTY;
  WDB[DATA_PHOTON_XS_EMAX] = -INFTY;

  /* Delayed nubar flag */

  WDB[DATA_USE_DELNU] = -1.0;

  /* Doppler broadening mode */

  WDB[DATA_TMS_MODE] = TMS_MODE_NONE;
  WDB[DATA_USE_DENSITY_FACTOR] = (double)NO;
  WDB[DATA_USE_DOPPLER_PREPROCESSOR] = (double)NO;

  /* Ures energy boundaries */

  WDB[DATA_URES_EMIN] = INFTY;
  WDB[DATA_URES_EMAX] = -INFTY;

  /* Critical spectrum calculation */

  WDB[DATA_B1_CALC] = (double)NO;
  WDB[DATA_B1_MODE] = (double)CRIT_SPECTRUM_OLD;
  WDB[DATA_B1_FM_DIFF] = (double)FM_DIFF_INF;
  WDB[DATA_B1_BURNUP_CORR] = (double)NO;
  WDB[DATA_B1_ERR_LIMIT] = 1E-7;
  WDB[DATA_B1_KEFF_TGT] = 1.0;
  WDB[DATA_B1_MAX_ITER] = 25.0;
  WDB[DATA_B1_B2_INIT] = 1E-6;

  /* Default micro-group structure */

  WDB[DATA_MICRO_PTR_EGRID] = (double)PutText("defaultmg");

  /* Batching for micro-group calculation (-1 = set equal to */
  /* number of cycles) */

  WDB[DATA_MICRO_CALC_BATCH_SIZE] = 20;

  /* Shared scoring buffer (NOTE: Tän pitää olla defaulttina NO, tai */
  /* sitä ei voi muuttaa readinput.c:stä) */

  WDB[DATA_OPTI_SHARED_BUF] = (double)NO;

  /* Shared RES2 array (NOTE: Tälle tehdään viritys initomp.c:ssä, */
  /* jotta arvoa voi muuttaa readinput.c:ssä). */

  WDB[DATA_OPTI_SHARED_RES2] = (double)YES;

  /* Reproducibility in OpenMP and MPI modes */

  WDB[DATA_OPTI_OMP_REPRODUCIBILITY] = (double)YES;
  WDB[DATA_OPTI_MPI_REPRODUCIBILITY] = (double)NO;

  /* Include scattering production in removal xs */

  WDB[DATA_GC_REMXS_MULT] = (double)YES;

  /* Solution relaxation factor for coupled calculation */

  WDB[DATA_SOL_REL_FACT] = 1.0;

  /* Solution relaxation maximum population for coupled calculation */

  WDB[DATA_SOL_REL_MAX_POP] = INFTY/1e6;

  /* Analog reaction rate calculation */

  WDB[DATA_ANA_RR_NCALC] = (double)ARR_MODE_NONE;
  WDB[DATA_ANA_RR_PCALC] = (double)ARR_MODE_NONE;

  /* Terminate on die */

  WDB[DATA_TERMINATE_ON_DIE] = (double)YES;

  /* Recorded tracks for plots */

  WDB[DATA_TRACK_PLOTTER_HIS] = -1.0;
  WDB[DATA_TRACK_PLOTTER_FILE] = (double)NO;

  /* Track plots */

  WDB[DATA_TRACK_PLOT_TMIN] = 0.0;
  WDB[DATA_TRACK_PLOT_TMAX] = INFTY;
  WDB[DATA_TRACK_PLOT_FRAMES] = 1.0;
  WDB[DATA_TRACK_PLOT_NHIS] = 100.0;

  /* UFS mode */

  WDB[DATA_UFS_MODE] = (double)UFS_MODE_NONE;

  /* Default cross section data library */

  if ((path = getenv("SERPENT_ACELIB")) != NULL)
    {
      ptr = ReallocMem(DATA_ARRAY, 2);
      WDB[DATA_PTR_ACEDATA_FNAME_LIST] = (double)ptr;
      WDB[ptr++] = (double)PutText(path);
      WDB[ptr] = NULLPTR;
    }

  /* Reset neutron, gamma and precursor buffer factors */

  WDB[DATA_PART_NBUF_FACTOR] = -1.0;
  WDB[DATA_PART_GBUF_FACTOR] = -1.0;
  WDB[DATA_PART_PBUF_FACTOR] = -1.0;

  /* Reaction sampling */

  WDB[DATA_NPHYS_SAMPLE_FISS] = (double)YES;
  WDB[DATA_NPHYS_SAMPLE_CAPT] = (double)YES;
  WDB[DATA_NPHYS_SAMPLE_SCATT] = (double)YES;

  /* Number of time bins */

  WDB[DATA_DYN_NB] = 1.0;

  /* Alpha eigenvalue */

  WDB[DATA_ALPHA_EIG] = 0.0; /* 1.431E+04 = EPR-laskun kriittinen alpha */

  /* Nubar scaling factor in external source simulation */

  WDB[DATA_EXT_SRC_NUBAR_F] = 1.0;

  /* Minimum xs for CFE */

  WDB[DATA_CFE_N_MIN_L] = 20.0;
  WDB[DATA_CFE_N_MIN_T] = -1.0;
  WDB[DATA_CFE_G_MIN_L] = 20.0;
  WDB[DATA_CFE_G_MIN_T] = -1.0;

  /* Energy boundary for cache-optimized xs block */

  WDB[DATA_CACHE_OPTI_EMAX] = 1E-5;

  /* Ures infinite dilution cut-off */

  WDB[DATA_URES_DILU_CUT] = 1E-9;

  /* Poison calculation */

  WDB[DATA_OPTI_POISON_CALC] = (double)NO;
  WDB[DATA_OPTI_POISON_CALC_XE135M] = (double)NO;

  /* Eddington factors */

  WDB[DATA_OPTI_EDDINGTON_CALC] = (double)NO;

  /* Equilibrium poison calculation */

  WDB[DATA_XENON_EQUILIBRIUM_MODE] = -1.0;
  WDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] = -1.0;

  /* MPI batch size */

  WDB[DATA_OPTI_MPI_BATCH_SIZE] = 10000.0;

  /* Actinide limits for burnup calculation */

  WDB[DATA_BU_ACT_MIN_Z] = 90.0;
  WDB[DATA_BU_ACT_MAX_Z] = 96.0;

  /* Number of progenies for beta-eff and prompt lifetime calculation */

  WDB[DATA_IFP_OPT_PRINT_ALL] = (double)NO;
  WDB[DATA_IFP_CHAIN_LENGTH] = 15.0;
  WDB[DATA_PERT_VAR_A] = 1E-8;
  WDB[DATA_PERT_VAR_C] = 1E-6;
  WDB[DATA_PERT_N_BATCH] = 9.0;

  /* ICM stuff */

  WDB[DATA_ICM_CALC] = (double)NO;
  WDB[DATA_ICM_NSEG] = -1.0;
  WDB[DATA_ICM_NSUB] = 1.0;
  WDB[DATA_ICM_NMU0] = 1.0;
  WDB[DATA_ICM_NMU1] = 1.0;
  WDB[DATA_ICM_NMU2] = 1.0;

  /* Group constant generation at multiple levels */

  WDB[DATA_MULTI_LEVEL_GCU] = (double)NO;

  /* Remove void cells and ignore undefined cells */

  WDB[DATA_REMOVE_VOID_CELLS] = (double)NO;
  WDB[DATA_IGNORE_UNDEFINED_CELLS] = (double)YES;

  /* Maximum number of undefined positions */

  WDB[DATA_MAX_UNDEF_POS] = 0.0;

  /* Print interval */

  WDB[DATA_PRINT_INTERVAL] = 50.0;

  /* Write/write restart file */

  WDB[DATA_WRITE_RESTART_FILE] = (double)NO;
  WDB[DATA_READ_RESTART_FILE] = (double)NO;
  WDB[DATA_RESTART_READ_IDX] = -1.0;
  WDB[DATA_RESTART_READ_ZERO_XE] = (double)NO;
  WDB[DATA_RESTART_READ_ZERO_SM] = (double)NO;
  WDB[DATA_RESTART_READ_CONTINUE] = (double)NO;
  WDB[DATA_RESTART_START_POINT] = -1.0;

  /* Statistical tests on group constants */

  WDB[DATA_GC_STAT_TESTS] = (double)NO;
  WDB[DATA_RUN_STAT_TESTS] = (double)NO;

  /* STL geometry stuff */

  WDB[DATA_STL_TEMP_ARRAY_SIZE] = 100.0;
  WDB[DATA_STL_TEST_N_PTS] = 0.0;
  WDB[DATA_STL_TEST_N_DIR] = 0.0;
  WDB[DATA_STL_GEOM_TEST_MODE] = (double)NO;
  WDB[DATA_STL_MESH_FILES] = (double)NO;

  /* NOTE: Jos tää on suurempi niin se voi aiheuttaa infinite loopin */
  /* ST-moodissa (esim pupu testikeissi). Pienempi arvo puolestaan   */
  /* voi aiheuttaa geometriaerroreita kun säteet osuu liian lähelle  */
  /* kolmion reunaa. */

  WDB[DATA_STL_FACET_EXD] = 1E-4;

  /* Setting this to YES is a fall-back option to enfore the use of */
  /* delta-tracking in STL search mesh cells with facets. */

  WDB[DATA_STL_ENFORCE_DT] = (double)NO;

  /* Maximum tracking loop */

  WDB[DATA_NEUTRON_MAX_TRACK_LOOP] = 1000000.0;
  WDB[DATA_PHOTON_MAX_TRACK_LOOP] = 1000000.0;

  /* Error flags */

  WDB[DATA_NEUTRON_MAX_TRACK_LOOP_ERR] = (double)YES;
  WDB[DATA_PHOTON_MAX_TRACK_LOOP_ERR] = (double)YES;

  /* Reset parameters for Wielandt method */

  WDB[DATA_WIELANDT_MODE] = (double)WIELANDT_MODE_NONE;
  WDB[DATA_WIELANDT_KEFF] = INFTY;
  WDB[DATA_WIELANDT_KP] = 0.0;
  WDB[DATA_WIELANDT_P] = 0.0;

  /* Response matrix calculation */

  WDB[DATA_RMTX_CALC] = NO;
  WDB[DATA_RMTX_MFP_CALC] = NO;

  /* No response-matrix based acceleration */

  WDB[DATA_RMX_CONVG_ACC] = (double)NO;

  /* Multi-response stuff for RMX solver */

  WDB[DATA_RMX_MULT_CONVG_LIM] = 0.01;
  WDB[DATA_RMX_MULT_CONVG_PASS] = 0.99;
  WDB[DATA_RMX_MULT_MAX_ITER] = 50.0;

  /* RMX test modes */

  WDB[DATA_RMX_TEST_MODE] = (double)NO;

  /* Reset coefficient calculation index and error mode */

  WDB[DATA_COEF_CALC_IDX] = -1.0;
  WDB[DATA_COEF_CALC_INCLUDE_ERRORS] = (double)NO;

  /* Reset flag for more coefficient calculations */

  WDB[DATA_MORE_COEF_CALC] = (double)NO;
  WDB[DATA_COEF_CALC_SPECIAL_MODE] = 0.0;

  /* History variation break point */

  WDB[DATA_HISV_BREAK_POINT] = 0.0;

  /* Case matrix calculation */

  WDB[DATA_CASEMTX_RUN] = (double)NO;

  /* Homogeneous flux solution */

  WDB[DATA_HOMOFLUX_SOLVER] = 1.0;
  WDB[DATA_ADF_TRAPZ_PT] = 100.0;
  WDB[DATA_HOMOFLUX_DIFFCOEF] = (double)HOMOFLUX_DIFF_INF;

  /* Klein-Nishina threshold energy */

  WDB[DATA_PHOTON_EKN] = INFTY;

  /* Use Doppler broadening of Compton photons */

  WDB[DATA_PHOTON_USE_DOPPLER] = (double)YES;

  /* Use TTB-approximation */

  WDB[DATA_PHOTON_USE_TTB] = (double)YES;

  /* Separate TTB data for positrons */

  WDB[DATA_PHOTON_TTBPM] = (double)YES;

  /* Energy threshold for linear interpolation of TTB number yield */

  WDB[DATA_PHOTON_TTBLINE] = 0.1;

  /* Detailed angular distribution for Compton electrons */

  WDB[DATA_PHOTON_COMP_EANG] = (double)NO;

  /* Default photon data library */

  if ((path = getenv("SERPENT_DATA")) != NULL) {
    sprintf(tmpstr, "%s/%s/", path, "photon_data");
    WDB[DATA_PHOTON_DATA_DIR] = (double)PutText(tmpstr);
  }

  /* Fallback option for Compton profiles (0: use 2.1.30, 1: use 2.1.29) */
  WDB[DATA_PHOTON_CP_FALLBACK_2129] = (double)NO;

  /* Photon data file names */

  WDB[DATA_PHOTON_COH_FNAME] = (double)PutText("cohff.dat");
  WDB[DATA_PHOTON_INCOH_FNAME] =(double) PutText("incohsf.dat");
  WDB[DATA_PHOTON_RELAX_FNAME] = (double)PutText("relax.dat");
  WDB[DATA_PHOTON_PEXS_FNAME] = (double)PutText("xspe.dat");
  WDB[DATA_PHOTON_PESS_FNAME] = (double)PutText("xspess.dat");
  WDB[DATA_PHOTON_PETOT_FNAME] = (double)PutText("xspetot.dat");
  WDB[DATA_PHOTON_CP_FNAME] = (double)PutText("ComptonProfiles.dat");

  /* Size of stopping power and TTB arrays */

  WDB[DATA_ELECTRON_SP_N] = 200;

  /* Stopping power debugging energy grid */

  WDB[DATA_ELECTRON_SP_E] = NULLPTR;

  /* Flag for printing electron and positron stopping powers */

  WDB[DATA_ELECTRON_PRINT_SP] = (double)NO;

  /* Calculate radiative stopping powers (instead of using data) */

  WDB[DATA_ELECTRON_SP_RAD_CALC] = (double)YES;

  /* Electron data file names */

  WDB[DATA_ELECTRON_SP_FNAME] = (double)PutText("el_stopping_power.dat");
  WDB[DATA_ELECTRON_BR_FNAME] = (double)PutText("pdebr.dat");
  WDB[DATA_ELECTRON_GSCONFIG_FNAME] = (double)PutText("gsconfig.dat");

  /* Decay source flag */

  WDB[DATA_USE_DECAY_SRC] = (double)NO;

  /* Weight windows */

  WDB[DATA_USE_WEIGHT_WINDOWS] = (double)NO;

  WDB[DATA_WWD_LOWER_BOUND] = -1.0;
  WDB[DATA_WWD_UPPER_BOUND] = -1.0;
  WDB[DATA_WWD_MAX_SPLIT] = 1E4;
  WDB[DATA_WWD_MIN_ROULETTE] = 1E-18;
  WDB[DATA_WWD_SURVIVAL_F] = 1.0;

  /* Heat, gamma and light nuclide production cross sections */

  WDB[DATA_INCLUDE_HEAT_PROD_XS] = (double)NO;
  WDB[DATA_INCLUDE_PHOT_PROD_XS] = (double)NO;
  WDB[DATA_INCLUDE_PROT_PROD_XS] = (double)NO;
  WDB[DATA_INCLUDE_DEUT_PROD_XS] = (double)NO;
  WDB[DATA_INCLUDE_TRIT_PROD_XS] = (double)NO;
  WDB[DATA_INCLUDE_HE3_PROD_XS] = (double)NO;
  WDB[DATA_INCLUDE_HE4_PROD_XS] = (double)NO;

  /* Factor for reading f*RDB[SRC_POP] precursors points for initial source  */
  /* in transient simulations (some of the precursor points might be removed */
  /* in population control) */

  WDB[DATA_PREC_SRC_FACT] = 10.0;

  /* Threshold for storing precursors generated by implicit estimator     */
  /* during transient simulation, negative value is based on prec. weight */
  /* Positive based on prec. emission during typical time interval */

  WDB[DATA_PREC_STORE_TRESH] = 1.0;

  /* Additional normalization factor for transient simulations. Normalization */
  /* is based on the criticality source simulation, but this can be used to   */
  /* multiply the normalization of the criticality source simulation */

  WDB[DATA_TRANS_NORM_FACT] = 1.0;

  /* File name flags for file-based particle stores to use in coupled */
  /* dynamic calculations, 0 = A, 1 = B */

  WDB[DATA_BOI_STORE_NAME] = 0.0;
  WDB[DATA_EOI_STORE_NAME] = 1.0;

  /* CMM method for diffusion coefficients */

  WDB[DATA_CMM_CALC] = (double)CMM_REMXS_IMPL;

  /* Response matrix calculation */

  WDB[DATA_RMTX_CALC] = (double)NO;
  WDB[DATA_RMTX_SOLVE_FORWARD] = (double)YES;
  WDB[DATA_RMTX_SOLVE_ADJOINT] = (double)YES;

  /* Variance reduction iterations */

  WDB[DATA_RUN_VR_ITER] = (double)NO;

  /* Geometry-based importances */

  WDB[DATA_USE_GEOM_IMP] = (double)NO;

  /* Exponential decay heat fit */

  WDB[DATA_EXPO_DEC_FIT_NF] = -1.0;
  WDB[DATA_EXPO_DEC_FIT_NT] = 1000.0;
  WDB[DATA_EXPO_DEC_FIT_TMIN] = 1E-3;
  WDB[DATA_EXPO_DEC_FIT_TMAX] = 6000;

  /* Particle speeds (used with track plot animation) */

  WDB[DATA_NEUTRON_SPD] = -1.0;
  WDB[DATA_PHOTON_SPD] = -1.0;

  /* Energy-dependent fission yields */

  WDB[DATA_FISSY_ENE_DEP] = (double)YES;

  /* Confidential flag */

  WDB[DATA_CONFIDENTIAL] = (double)NO;

  /* Energy cut-offs off */

  WDB[DATA_NEUTRON_ECUT] = -INFTY;
  WDB[DATA_PHOTON_ECUT] = RDB[DATA_PHOTON_EMIN];

  /* Photon mean free path cut-off */

  WDB[DATA_PHOTON_MFPCUT] = -INFTY;

  /* Source file buffer size */

  WDB[DATA_SRC_FILE_BUF_SIZE] = 1000000.0;

  /* Source importance calculation */

  WDB[DATA_SRC_IMP_CALC] = (double)NO;

  /* Maximum running time */

  WDB[DATA_MAX_TRANSPORT_RUNTIME] = -1.0;
  WDB[DATA_STOP_AFTER_PROCESSING] = (double)NO;

  /* Options for sensitivity calculations */

  WDB[DATA_SENS_MODE] = SENS_MODE_NONE;
  WDB[DATA_SENS_LAST_GEN] = 15.0;
  WDB[DATA_SENS_SCORE_TYPE] = (double)SENS_SCORE_TYPE_EVENT;

  /* Elemental decomposition */

  WDB[DATA_ELEM_DECOMP] = (double)NO;
  WDB[DATA_ELEM_DECOMP_PTR_LIST] = NULLPTR;

  /* Options for energy deposition */

  WDB[DATA_EDEP_MODE] = (double)EDEP_MODE_CONSTANT;
  WDB[DATA_EDEP_CAPT_E] = 0.0;
  WDB[DATA_EDEP_DELAYED] = (double)YES;
  WDB[DATA_EDEP_KEFF_CORR] = (double)YES;
  WDB[DATA_EDEP_LOCAL_EGD] = (double)NO;

  /* Domain decomposition */

  WDB[DATA_DD_DECOMPOSE] = (double)NO;
  WDB[DATA_DD_BUFF_NUM_PARTS] = 128;
  WDB[DATA_DD_SECT0] = 0.0;
  WDB[DATA_DD_ORIG_X0] = 0.0;
  WDB[DATA_DD_ORIG_Y0] = 0.0;
  WDB[DATA_DD_ORIG_Z0] = 0.0;

  /* On-the-fly burnup mode */

  WDB[DATA_OTF_BURN_MODE] = (double)NO;

  /* Global density factor */

  WDB[DATA_GLOBAL_DF] = 1.0;

  /* Common particle ques in OpenMP mode */

  WDB[DATA_COMMON_QUE_LIM] = 0.0;

  /* Default parameters for growing population simulation */

  WDB[DATA_GROW_POP_SIM] = (double)NO;
  WDB[DATA_GROW_POP_C] = 1.0;
  WDB[DATA_GROW_POP_MAX_POP] = 1E6;

  /* Randomize decay and fission yield data */

  WDB[DATA_BURN_RANDOMIZE_DEC] = (double)NO;
  WDB[DATA_BURN_RANDOMIZE_FY] = (double)NO;

  /* Energy-dependent importances */

  WDB[DATA_ENE_IMP_CALC] = (double)NO;
  WDB[DATA_ENE_IMP_POW] = -0.5;

  /***************************************************************************/
}

/*****************************************************************************/
