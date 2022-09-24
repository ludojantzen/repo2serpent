/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : preparetransportcycle.c                        */
/*                                                                           */
/* Created:       2011/05/23 (JLe)                                           */
/* Last modified: 2019/11/10 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calls all subroutines that need to be run before the         */
/*              transport cycle                                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrepareTransportCycle:"

/*****************************************************************************/

void PrepareTransportCycle()
{
  double mem;
  long ptr, type, trans, id;

  /* Reset cycle-wise MPI timer */

  ResetTimer(TIMER_MPI_OVERHEAD);

  /* Start process timers */

  ResetTimer(TIMER_PROCESS);
  StartTimer(TIMER_PROCESS);
  StartTimer(TIMER_PROCESS_TOTAL);

  /* Memory size before processing */

  mem = RDB[DATA_TOTAL_BYTES];

  /* Get burnup step type (not set in transport mode) */

  type = (long)RDB[DATA_BURN_STEP_TYPE];

  /* Check type and set transport flag */

  if ((type == DEP_STEP_DEC_STEP) || (type == DEP_STEP_DEC_TOT) ||
      (type == DEP_STEP_ACT_STEP) || (type == DEP_STEP_ACT_TOT))
    trans = NO;
  else
    trans = YES;

  /* Check VR iteration index */

  if ((long)RDB[DATA_VR_ITER_IDX] > 0)
    trans = NO;

  /* Reset poison concenctration and OTF burn data */

  if (trans == YES)
    {
      ResetPoisonConc();
      ResetOTFBurn();
    }

  /* Normalize compositions */

  NormalizeCompositions();

  /* Update average densities for micro depletion */

  UpdateMicroDens();

  /* Check step type */

  if (trans == YES)
    {
      /* Generate reaction lists */

      ProcessReactionLists();

      /* Calculate total cross sections */

      MaterialTotals();

      /* Distribute material data from to parallel MPI tasks (burnable */
      /* material compositions are located in different position and   */
      /* distributed in CollectBurnData(), called from BurnupCycle().  */

      DistributeMaterialData();

      /* Print for debugging */

      PrintReactionLists();

      /* Energy groups for spectrum-collapse method */

      ProcessBurnupEGroups();

      /* Calculate majorants for delta-tracking */

      CalculateDTMajorants();

      /* Plot XS data */

      XSPlotter();
    }

  /* Calculate masses */

  CalculateMasses();

  /* Calculate activities (toi if-lause on siltä varalta että rikotaan  */
  /* jotain, aktiivisuudet lasketaan muuten transportcycle.c:n lopussa. */
  /* Kutsua tarvitaan RadGammaSrc():n takia. Lisätty 16.5.2015 / JLE).  */

  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
    CalculateActivities();

  /* Process photon energy cut-offs */

  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
    ProcessPhotonEcut();

  fprintf(outp, "Clearing results and statistics...\n");

  /* Clear temporary data */

  ClearValuePairs();

  /* Reset statistics */

  ClearStat(-1);

  /* Reset transmutation reaction rates etc. */

  ClearTransmuXS();

  /* Clear micro-group data */

  ClearMicroGroupXS();

  /* Check memory size */

  if (((long)RDB[DATA_BURN_STEP] > 0) && (mem != RDB[DATA_TOTAL_BYTES]))
    Die(FUNCTION_NAME, "Change in total memory size");

  /* Clear some counters */

  ptr = (long)RDB[DATA_PTR_OMP_HISTORY_COUNT];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  ClearPrivateData(ptr);

  /* Reset B1 iteration counters */

  WDB[DATA_B1_REPEATED] = 0.0;
  WDB[DATA_B1_CONVERGED] = 0.0;

  /* Reset collision counters */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    PutPrivateData(ptr, 0, id);

  /* Reset relaxed transmutation reaction rates */

  ClearRelTransmuXS();

  /* Reset lost particle counter */

  WDB[DATA_UNDEF_POS_COUNT] = 0.0;

  /* Reset counter for solution relaxation */

  if (RDB[DATA_BURN_SIE] == (double)NO)
    {
      /* Not following the SIE iteration flow */

      WDB[DATA_SOL_REL_ITER] = 0.0;
      WDB[DATA_SOL_REL_NTOT] = 0.0;
      WDB[DATA_SOL_REL_NCUR] = RDB[DATA_CRIT_CYCLES]*RDB[DATA_CRIT_POP];
      WDB[DATA_SOL_REL_N1] = RDB[DATA_CRIT_CYCLES]*RDB[DATA_CRIT_POP];

    }
  else
    {
      /* Following SIE iteration where the coupled calculation power */
      /* is relaxed between corrector iterations */

      /* These are only reset if SIE moves to next step */

      if (RDB[DATA_BURN_CI_I] == 0)
       {
         WDB[DATA_SOL_REL_ITER] = 0.0;
          WDB[DATA_SOL_REL_NTOT] = 0.0;
          WDB[DATA_SOL_REL_NCUR] = RDB[DATA_CRIT_CYCLES]*RDB[DATA_CRIT_POP];
          WDB[DATA_SOL_REL_N1] = RDB[DATA_CRIT_CYCLES]*RDB[DATA_CRIT_POP];
       }
    }

  if (RDB[DATA_RUN_CC] == (double)YES)
    {
      /* Set iteration flag on */

      WDB[DATA_ITERATE] = (double)YES;
    }

  fprintf(outp, "OK.\n\n");

  /* Stop process timers */

  StopTimer(TIMER_PROCESS);
  StopTimer(TIMER_PROCESS_TOTAL);
}

/*****************************************************************************/
