/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : preparecciter.c                                */
/*                                                                           */
/* Created:       2014/10/15 (VVa)                                           */
/* Last modified: 2019/02/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calls all subroutines that need to be run before starting    */
/*              coupled calculation iteration of a transport cycle           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrepareCCIter:"

/*****************************************************************************/

void PrepareCCIter()
{

  long ptr, id;

  /* Return if not running coupled calculation */

  if(RDB[DATA_RUN_CC] == (double)NO)
    return;

  /* Start process timers */

  ResetTimer(TIMER_PROCESS);
  StartTimer(TIMER_PROCESS);
  StartTimer(TIMER_PROCESS_TOTAL);

  fprintf(outp, "Clearing results and statistics...\n");

  /* Reset statistics */

  ClearStat(-1);

  /* Reset transmutation reaction rates etc. */

  ClearTransmuXS();

  /* Clear micro-group data */

  ClearMicroGroupXS();

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

  fprintf(outp, "OK.\n\n");

  /* Stop process timers */

  StopTimer(TIMER_PROCESS);
  StopTimer(TIMER_PROCESS_TOTAL);

}

/*****************************************************************************/
