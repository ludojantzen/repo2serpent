/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setoptimization.c                              */
/*                                                                           */
/* Created:       2011/07/17 (JLe)                                           */
/* Last modified: 2019/09/06 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: - Set various options based on optimization                  */
/*                                                                           */
/* Comments: - Tää aliohjelma tekee oikeastaan kaikkea muutakin kuin         */
/*             asettaa optimoinnin, mm. korjaa käyttäjän tekemiä valintoja   */
/*             jne... Nimen voisi muuttaa.                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetOptimization:"

/*****************************************************************************/

void SetOptimization()
{
  long n, ptr;
  char tmpstr[MAX_STR];
  long maxg;

  /***************************************************************************/

  /***** Set options etc. ****************************************************/

  /* Track plotter */

  if ((long)RDB[DATA_STOP_AFTER_PLOT] == STOP_AFTER_PLOT_TRACKS)
    {
      /* Check number of tracks */

      if (RDB[DATA_TRACK_PLOTTER_HIS] < 1.0)
        Die(FUNCTION_NAME, "Invalid number of tracks");

      /* Reset burnup mode */

      WDB[DATA_BURNUP_CALCULATION_MODE] = (double)NO;

      /* No importance solver */

      WDB[DATA_RMTX_CALC] = (double)NO;

      /* Set number of histories */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        {
          WDB[DATA_CRIT_POP] = RDB[DATA_TRACK_PLOTTER_HIS];
          WDB[DATA_CRIT_CYCLES] = 0.0;
          WDB[DATA_CRIT_SKIP] = 1.0;
        }
      else
        {
          WDB[DATA_SRC_POP] = RDB[DATA_TRACK_PLOTTER_HIS];
          WDB[DATA_SRC_BATCHES] = 1.0;
          WDB[DATA_CRIT_SKIP] = 0.0;
        }

      /* Set flag */

      SetOption(DATA_EVENT_RECORD_FLAGS, RECORD_EVENT_PLOTTER);

      /* Reset minimum collision frequency */

      WDB[DATA_CFE_N_MIN_L] = INFTY;
      WDB[DATA_CFE_G_MIN_L] = INFTY;
    }

  /* RIA Simulation */

  if ((long)RDB[DATA_PTR_RIA0] > VALID_PTR)
    {
      /* Set simulation mode to criticality */

      WDB[DATA_SIMULATION_MODE] = (double)SIMULATION_MODE_CRIT;

      /* Check parameters */

      if ((long)RDB[DATA_SRC_POP] < 1)
        Error(0, "Source batch size must be defined");
      if ((long)RDB[DATA_SRC_BATCHES] < 1)
        Error(0, "Number of source batches cycles must be defined");

      /* Write source to file */

      sprintf(tmpstr, "%s.src", GetText(DATA_PTR_INPUT_FNAME));
      WDB[DATA_PTR_CRIT_SRC_DET] = (double)PutText(tmpstr);
    }

  /* Add k-eff iteration cycles to skip cycles */

  if ((long)RDB[DATA_ITER_MODE] != ITER_MODE_NONE)
    {
      if ((long)RDB[DATA_ITER_FIX] == YES)
        WDB[DATA_CRIT_SKIP] = 2.0*RDB[DATA_CRIT_SKIP] + RDB[DATA_ITER_NCYC];
      else
        WDB[DATA_CRIT_SKIP] = RDB[DATA_CRIT_SKIP] + RDB[DATA_ITER_NCYC];
    }

  /* Additional cycles in user-defined mode */

  if ((long)RDB[DATA_ITER_MODE] == ITER_MODE_USER)
    {
      /* Check MPI mode */

#ifdef MPI_MODE2

      Die(FUNCTION_NAME, "MPI mode 2 not supported");

#endif

      /* Add active cycles */

      WDB[DATA_CRIT_CYCLES] = RDB[DATA_CRIT_CYCLES] + 3.0*RDB[DATA_ITER_NCYC];
    }

  /* Check that simulation mode is set */

  if ((long)RDB[DATA_SIMULATION_MODE] < 0)
    Error(0, "Simulation mode must be set using \"pop\", \"nps\" or \"dyn\"");
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /* Check parameters */

      if ((long)RDB[DATA_CRIT_POP] < 1)
        Error(0, "Population size must be defined");
      if ((long)RDB[DATA_CRIT_CYCLES] + (long)RDB[DATA_CRIT_SKIP] < 1)
        Error(0, "Number of criticality cycles must be defined");
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC)
    {
      /* Check parameters */

      if ((long)RDB[DATA_SRC_POP] < 1)
        Error(0, "Source batch size must be defined");
      if ((long)RDB[DATA_SRC_BATCHES] < 1)
        Error(0, "Number of source batches cycles must be defined");

      /* Reset number of skip cycles */

      WDB[DATA_CRIT_SKIP] = 0.0;

      /* External source mode + precursor tracking = DELDYN mode */

      if ((long)RDB[DATA_PTR_PREC_DET] > VALID_PTR)
        WDB[DATA_SIMULATION_MODE] = (double)SIMULATION_MODE_DELDYN;

      /* Time dependent simulation + coupled calculation = DYN mode */

      if (RDB[DATA_RUN_CC] == (double)YES)
        WDB[DATA_SIMULATION_MODE] = (double)SIMULATION_MODE_DYN;
    }
  else
    Die(FUNCTION_NAME, "Invalid simulation mode");

  /* Set output interval to 1 to print output after each time interval */
  /* in coupled time dependent (transient) simulations since results   */
  /* are cleared after each time interval */

  if (RDB[DATA_SIMULATION_MODE] == (double)SIMULATION_MODE_DYN)
    WDB[DATA_PRINT_INTERVAL] = 1.0;

  /* Include total list in mode 0 */

  WDB[DATA_OPTI_MODE0_INCLUDE_TOTAL] = 1.0;

  /* Set unionized energy grid thinning tolerance */

  if (RDB[DATA_ERG_TOL] < 0.0)
    {
      if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
        WDB[DATA_ERG_TOL] = 0.0;
      else
        WDB[DATA_ERG_TOL] = 5E-5;
    }

  /* Check optimization mode */

  switch ((long)RDB[DATA_OPTI_MODE])
    {
    case 1:
      {
        WDB[DATA_OPTI_RECONSTRUCT_MICROXS] = (double)NO;
        WDB[DATA_OPTI_RECONSTRUCT_MACROXS] = (double)NO;
        WDB[DATA_OPTI_IMPLICIT_RR] = (double)NO;
        WDB[DATA_OPTI_GC_CALC] = (double)NO;

        if ((long)RDB[DATA_OPTI_MG_MODE] < 0)
          WDB[DATA_OPTI_MG_MODE] = (double)NO;

        if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES)
          Note(0, "Option 'set xscalc 2' ignored in optimization mode 1");

        if ((long)RDB[DATA_OPTI_DIX] < 0)
          WDB[DATA_OPTI_DIX] = (double)NO;

        WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)NO;

        break;
      }
    case 2:
      {
        WDB[DATA_OPTI_RECONSTRUCT_MICROXS] = (double)YES;
        WDB[DATA_OPTI_RECONSTRUCT_MACROXS] = (double)NO;
        WDB[DATA_OPTI_IMPLICIT_RR] = (double)NO;
        WDB[DATA_OPTI_GC_CALC] = (double)NO;

        if ((long)RDB[DATA_OPTI_MG_MODE] < 0)
          WDB[DATA_OPTI_MG_MODE] = (double)YES;

        if ((long)RDB[DATA_OPTI_DIX] < 0)
          WDB[DATA_OPTI_DIX] = (double)NO;

        if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES)
          Note(0, "Option 'set xscalc 2' ignored in optimization mode 2");

        WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)NO;

        break;
      }
    case 3:
      {
        WDB[DATA_OPTI_RECONSTRUCT_MICROXS] = (double)NO;
        WDB[DATA_OPTI_RECONSTRUCT_MACROXS] = (double)YES;
        WDB[DATA_OPTI_IMPLICIT_RR] = (double)YES;
        WDB[DATA_OPTI_GC_CALC] = (double)YES;

        if ((long)RDB[DATA_OPTI_MG_MODE] < 0)
          WDB[DATA_OPTI_MG_MODE] = (double)NO;

        if ((long)RDB[DATA_OPTI_DIX] < 0)
          WDB[DATA_OPTI_DIX] = (double)YES;

        if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] < 0)
          WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)YES;

        break;
      }
    case 4:
      {
        WDB[DATA_OPTI_RECONSTRUCT_MICROXS] = (double)YES;
        WDB[DATA_OPTI_RECONSTRUCT_MACROXS] = (double)YES;
        WDB[DATA_OPTI_IMPLICIT_RR] = (double)YES;
        WDB[DATA_OPTI_GC_CALC] = (double)YES;

        if ((long)RDB[DATA_OPTI_MG_MODE] < 0)
          WDB[DATA_OPTI_MG_MODE] = (double)NO;

        if ((long)RDB[DATA_OPTI_DIX] < 0)
          WDB[DATA_OPTI_DIX] = (double)NO;

        if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] < 0)
          WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)YES;

        break;
      }
    default:
      Error(0, "Invalid optimization mode %ld", (long)RDB[DATA_OPTI_MODE]);
    }

  /* Switch to multi-group TMS mode if mg mode */

  if (((long)RDB[DATA_OPTI_MG_MODE] == YES) &&
      ((long)RDB[DATA_TMS_MODE] == TMS_MODE_CE))
    WDB[DATA_TMS_MODE] = (double)TMS_MODE_MG;

  /* Set unionization flag */

  if (((long)RDB[DATA_OPT_USE_DT] == YES) ||
      ((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == YES) ||
      ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == YES) ||
      ((long)RDB[DATA_OPTI_DIX] == YES))
    WDB[DATA_OPTI_UNIONIZE_GRID] = (double)YES;
  else
    WDB[DATA_OPTI_UNIONIZE_GRID] = (double)NO;

  /* No reconstruction of microscopic xs if double-indexing is on */

  if ((long)RDB[DATA_OPTI_DIX] == YES)
    WDB[DATA_OPTI_RECONSTRUCT_MICROXS] = (double)NO;

  /* Set group constant generation on if B1 is used */

  if ((long)RDB[DATA_B1_CALC] == YES)
    WDB[DATA_OPTI_GC_CALC] = (double)YES;

  /* Set group constant calculation on if universe is given and off if */
  /* set to null */

  if ((long)RDB[DATA_PTR_GCU0] > 0)
    WDB[DATA_OPTI_GC_CALC] = (double)YES;
  else if ((long)RDB[DATA_PTR_MORA0] > 0)
    WDB[DATA_OPTI_GC_CALC] = (double)YES;
  else if ((long)RDB[DATA_PTR_GCU0] < 0)
    WDB[DATA_OPTI_GC_CALC] = (double)NO;

  /* Set group constant calculation off if not neutron transport mode */

  if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == NO)
    {
      WDB[DATA_OPTI_GC_CALC] = (double)NO;
      WDB[DATA_B1_CALC] = (double)NO;
      WDB[DATA_PTR_GCU0] = -1.0;
    }

  /* Switch CMM calculation off if GC calculation is not on */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    WDB[DATA_CMM_CALC] = (double)NO;

  /* Use implicit reaction rates if group constants are generated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == YES)
    WDB[DATA_OPTI_IMPLICIT_RR] = (double)YES;

  /* Stop particle tracks at outer geometry boundary if wwgen or wwin is on */

  if (((long)RDB[DATA_USE_WEIGHT_WINDOWS] == YES) ||
      ((long)RDB[DATA_RMTX_CALC] == YES))
    WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

  /* Check that group constant calculation is used with fum */

  if (((long)RDB[DATA_OPTI_GC_CALC] == NO) && ((long)RDB[DATA_B1_CALC] == YES))
    Error(0, "Group constant calculation needed for B1 mode");

  /* Check that group constant calculation is used with mdep */

  if (((long)RDB[DATA_OPTI_GC_CALC] == NO) &&
      ((long)RDB[DATA_PTR_MDEP0] > VALID_PTR))
    Error(0, "Group constant generation must be on with mdep");

  /* Check CFE collision distance when GC calculation is on */

  if (((long)RDB[DATA_OPTI_GC_CALC] == YES) &&
      (RDB[DATA_CFE_N_MIN_L] > 100.0) &&
      ((long)RDB[DATA_STOP_AFTER_PLOT] == STOP_AFTER_PLOT_NONE))
    Error(0,
          "Insufficient CFE collision distance for group constant generation");

  /* Check on-the-fly mode */

  if ((long)RDB[DATA_TMS_MODE] != TMS_MODE_NONE)
    {
      /* Spectrum-collapse off */

      WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)NO;

      /* Grid thinning off */

      WDB[DATA_ERG_TOL] = 0.0;
    }

  /* Set delayed nubar flag if not set in input */

  if ((long)RDB[DATA_USE_DELNU] == -1)
    {
      /* Check simulation mode */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        WDB[DATA_USE_DELNU] = (double)YES;
      else
        WDB[DATA_USE_DELNU] = (double)NO;
    }

  /* Reset replay option if no reproducibility */

  if ((mpitasks > 1) && ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO))
    WDB[DATA_OPTI_REPLAY] = (double)NO;

  if (((long)RDB[DATA_OMP_MAX_THREADS] > 1) &&
      ((long)RDB[DATA_OPTI_OMP_REPRODUCIBILITY] == NO))
    WDB[DATA_OPTI_REPLAY] = (double)NO;

  /* Disable grid thinning if problems */

  if (((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == NO) ||
      (((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == NO) &&
       ((long)WDB[DATA_OPTI_MG_MODE] == NO)))
    WDB[DATA_ERG_TOL] = 0.0;

  /* Set defaults for implicit reactions */

  if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == NO)
    {
      /* No weight windows, use implcit (n,xn) by default */

      if (RDB[DATA_OPT_IMPL_FISS] == -1.0)
        WDB[DATA_OPT_IMPL_FISS] = (double)NO;

      if (RDB[DATA_OPT_IMPL_CAPT] == -1.0)
        WDB[DATA_OPT_IMPL_CAPT] = (double)NO;

      if (RDB[DATA_OPT_IMPL_NXN] == -1.0)
        WDB[DATA_OPT_IMPL_NXN] = (double)YES;
    }
  else
    {
      /* Weight windows, set all to analog by default */

      if (RDB[DATA_OPT_IMPL_FISS] == -1.0)
        WDB[DATA_OPT_IMPL_FISS] = (double)NO;

      if (RDB[DATA_OPT_IMPL_CAPT] == -1.0)
        WDB[DATA_OPT_IMPL_CAPT] = (double)NO;

      if (RDB[DATA_OPT_IMPL_NXN] == -1.0)
        WDB[DATA_OPT_IMPL_NXN] = (double)NO;

      /* Switch group constant generation off unless universe is specified */

      if ((long)RDB[DATA_PTR_GCU0] < VALID_PTR)
        {
          WDB[DATA_OPTI_GC_CALC] = (double)NO;
          WDB[DATA_B1_CALC] = (double)NO;
        }
    }

  /* Disable shared results array if source biasing is in use (ei   */
  /* ole aavistustakaan miksi tämä tehdään, mutta se kaataa RES2-   */
  /* muistinvarauksen jos moodi on initdata.c:ssä YES, ja inputissa */
  /* on mesh plotteja. */

  /*
  if ((long)RDB[DATA_UFS_MODE] != UFS_MODE_NONE)
    WDB[DATA_OPTI_SHARED_RES2] = (double)NO;
  */
  /* Set explicit fission if criticality source simulation */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    WDB[DATA_OPT_IMPL_FISS] = (double)NO;

  /* Set neutron buffer factor */

  if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
    {
      /* Check value */

      if ((long)RDB[DATA_PART_NBUF_FACTOR] < 100)
        WDB[DATA_PART_NBUF_FACTOR] = 100.0;
    }
  else if ((long)RDB[DATA_PART_NBUF_FACTOR] < 0.0)
    {
      /* Set default */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        WDB[DATA_PART_NBUF_FACTOR] = 5.0;
      else
        WDB[DATA_PART_NBUF_FACTOR] = 1000.0;

      /* Check if too low */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        if (RDB[DATA_CRIT_POP]/RDB[DATA_OMP_MAX_THREADS] < 1000.0)
          WDB[DATA_PART_NBUF_FACTOR] = RDB[DATA_PART_NBUF_FACTOR]*
            1000.0/RDB[DATA_CRIT_POP]*RDB[DATA_OMP_MAX_THREADS];
    }

  /* Set gamma buffer factor */

  if (RDB[DATA_PART_GBUF_FACTOR] < 0.0)
    WDB[DATA_PART_GBUF_FACTOR] = 1000;

  /* Override micro-group batching for group constant generation */

  if ((long)RDB[DATA_MICRO_CALC_BATCH_SIZE] == -1)
    {
      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        WDB[DATA_MICRO_CALC_BATCH_SIZE] = RDB[DATA_CRIT_CYCLES];
      else
        WDB[DATA_MICRO_CALC_BATCH_SIZE] = RDB[DATA_SRC_BATCHES];
    }

  if ((long)RDB[DATA_BATCH_INTERVAL] > 1)
    WDB[DATA_MICRO_CALC_BATCH_SIZE] = RDB[DATA_BATCH_INTERVAL];

  /* Check implicit capture and multi-group mode */

  if (((long)RDB[DATA_OPTI_MG_MODE] == YES) &&
      ((long)RDB[DATA_OPT_IMPL_CAPT] == YES))
    Error(0, "Implicit capture does not work in this optimization mode");

  /* Wieland shift */

  if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
    {
      /* Sanity check of input */

      if ((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT)
        Die(FUNCTION_NAME,
            "Wielandt shift can only be used with criticality source mode");

      /* Check multiple eigenvalue simulations */

      if ((long)RDB[DATA_N_POP_EIG] > 1)
        Error(0, "Wielandt method not allowed with multiple eig. simulations");

      /* Compute initial guess of shifted eigenvalue */

      if ((long)RDB[DATA_WIELANDT_MODE] == WIELANDT_MODE_FIX_K)
        {
          WDB[DATA_CYCLE_KEFF] =
            (RDB[DATA_CYCLE_KEFF]*RDB[DATA_WIELANDT_KEFF])/
            (RDB[DATA_WIELANDT_KEFF] - RDB[DATA_CYCLE_KEFF]);

          /* Check */

          if (RDB[DATA_CYCLE_KEFF] < 0.3)
            WDB[DATA_CYCLE_KEFF] = 0.3;
        }
    }

  /* Check geometry importances and DT */

  if ((long)RDB[DATA_USE_GEOM_IMP] == YES)
    {
      /* Check for neutrons */

      if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
        if (RDB[DATA_DT_NTHRESH] < 1.0)
          Error(0, "Geometry importances cannot be used with delta-tracking");

      /* Check for photons */

      if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
        if (RDB[DATA_DT_PTHRESH] < 1.0)
          Error(0, "Geometry importances cannot be used with delta-tracking");
    }

  /* Check analog (n,xn) with group constant generation */

  if (((long)RDB[DATA_OPTI_GC_CALC] == YES) &&
      ((long)RDB[DATA_OPT_IMPL_NXN] == NO))
    Error(0, "Analog (n,xn) does not work with group constant generation");

  /* Set weight cut-off if implicit capture is on */

  if ((long)RDB[DATA_OPT_IMPL_CAPT] == YES)
    if (RDB[DATA_OPT_ROULETTE_W0] < 0.0)
      WDB[DATA_OPT_ROULETTE_W0] = 0.001;

  /* No tracking loop errors if DT is enforced in STL geometries */

  if ((long)RDB[DATA_STL_ENFORCE_DT] == YES)
    {
      /* Error flags */

      WDB[DATA_NEUTRON_MAX_TRACK_LOOP_ERR] = (double)NO;
      WDB[DATA_PHOTON_MAX_TRACK_LOOP_ERR] = (double)NO;
    }

  /* Check source point importance calculation */

   if ((long)RDB[DATA_SRC_IMP_CALC] == YES)
     {
       /* Reset size */

       WDB[DATA_CRIT_POP] = 1.0;
       WDB[DATA_CRIT_CYCLES] = 0.0;
       WDB[DATA_SRC_POP] = 1.0;
       WDB[DATA_SRC_BATCHES] = 1.0;
     }

   /* Set multiple particle transport mode */

   n = (long)(RDB[DATA_NEUTRON_TRANSPORT_MODE] +
              RDB[DATA_PHOTON_TRANSPORT_MODE]);

   if (n == 0)
     {
       /* Set default mode */

       WDB[DATA_NEUTRON_TRANSPORT_MODE] = (double)YES;
     }
   else if (n > 1)
     {
       /* Set multi-particle flag */

       WDB[DATA_MULTI_PARTICLE_TRANSPORT] = (double)YES;
     }

   /* Check if weight windows are used without normalization */

   if (((long)RDB[DATA_USE_WEIGHT_WINDOWS] == YES) &&
       ((long)RDB[DATA_PTR_NORM] < VALID_PTR) &&
       ((long)RDB[DATA_USE_DECAY_SRC] == NO))
     Note(0, "Normalization should be defined when using weight windows");

   /* Figure out the number of event generations to store */

   if (!((long)RDB[DATA_EVENT_RECORD_FLAGS] & RECORD_EVENT_PLOTTER) &&
       ((long)RDB[DATA_EVENT_RECORD_FLAGS] & RECORD_EVENT_SENS))
     {
       /* Recording events for Sens but not for plotter */

       /* Maximum generation is the maximum of Sens and IFP */

       maxg = (long)RDB[DATA_IFP_CHAIN_LENGTH];

       if ((long)RDB[DATA_SENS_LAST_GEN] > maxg)
         maxg = (long)RDB[DATA_SENS_LAST_GEN];

       /* Put maximum generation */

       WDB[DATA_EVENT_MAX_GEN] = (double)maxg;
     }

   /* Check number of normalizations */

   if ((ptr = (long)RDB[DATA_BURN_PTR_DEP]) < VALID_PTR)
     {
       /* No burnup intervals, check number of normalizations */

       if ((ptr = (long)RDB[DATA_PTR_NORM]) > VALID_PTR)
         if (NextItem(ptr) > VALID_PTR)
           Error(ptr, "Multiple normalizations");
     }
   else
     {
       /* Get number of burnup intervals */

       n = ListSize(ptr);

       /* Check number of normalizations */

       if ((ptr = (long)RDB[DATA_PTR_NORM]) > VALID_PTR)
         if (ListSize(ptr) > n)
           Error(ptr,
                 "Number of normalizations exceeds number of burnup intervals");
     }

   /* Check that restart files are not used with domain decomposition */

   if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
       (((long)RDB[DATA_WRITE_RESTART_FILE] == YES) ||
        ((long)RDB[DATA_READ_RESTART_FILE] == YES)))
     Error(0, "Domain decomposition does not yet work with restart files");

   /* Switch DD off if only 1 mpi task */

   if (((long)RDB[DATA_DD_DECOMPOSE] == YES) && (mpitasks == 1))
     WDB[DATA_DD_DECOMPOSE] = (double)NO;

   /* Set IFP chain length in DD */

   if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
     WDB[DATA_IFP_CHAIN_LENGTH] = -1.0;

   /* Check that materials are divided by div cards in dd */

   if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
       ((long)RDB[DATA_PTR_DIV0] < VALID_PTR))
     Error(0, "Domain decomposition without material division");

   /* Check that divided compositions are not written with DD */

   if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
       ((long)RDB[DATA_BURN_MAT_OUTPUT] != BURN_OUT_MAT_PARENT))
     Error(0, "Parameter \"set depout\" must be 2 with domain decomposition");

   /* Check OTF burnup mode and optimization mode 2 (syy miksi   */
   /* toi antaa vääriä tuloksia toistaiseksi tuntematon, ehkä    */
   /* noi reaktiot tulee huomioitua kahteen kertaan, tai ainakin */
   /* keff:issä ja ryhmävakioissa näkyy isoja eroja, tosin xenon */
   /* iteraatio näyttää toimivan Jle / 2.1.31 / 26.3.2018). */

   if (((long)RDB[DATA_OTF_BURN_MODE] == YES) &&
       ((long)RDB[DATA_OPTI_MG_MODE] == YES))
     Error(0, "OTF burnup solver does not work in optimization mode 2");

   /* Remember setting for implicit reaction rates */

   WDB[DATA_OPTI_IMPLICIT_RR0] = RDB[DATA_OPTI_IMPLICIT_RR];

   /* Group constant generation of if track plotter */

   if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
     {
       WDB[DATA_OPTI_GC_CALC] = (double)NO;
       WDB[DATA_PTR_GCU0] = 0.0;
     }

   /* Check if coefficient calculation is on but gc generation is not */

   if (((long)RDB[DATA_PTR_COEF0] > VALID_PTR) &&
       ((long)RDB[DATA_PTR_GCU0] < VALID_PTR))
     Note(0, "Group constant generation is not switched on");

   /* Test ures + implicit capture + no fissions */

   if (((long)RDB[DATA_USE_URES] == YES) &&
       ((long)RDB[DATA_OPT_IMPL_CAPT] == YES) &&
       ((long)RDB[DATA_NPHYS_SAMPLE_FISS] == NO))
     Error(0, "Implicit capture must be off when fissions are not sampled");

  /***************************************************************************/
}

/*****************************************************************************/
