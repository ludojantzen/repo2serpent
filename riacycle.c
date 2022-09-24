/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : riacycle.c                                     */
/*                                                                           */
/* Created:       2012/10/24 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Runs a criticality source calculation, followed by time-     */
/*              dependent source simulation.                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RIACycle:"

/*****************************************************************************/

void RIACycle()
{
  long ptr;

  fprintf(outp, "Running RIA simulation...\n\n");

  /***************************************************************************/

  /***** First simulation (criticality) **************************************/

  /* Switch group constant generation off */

  WDB[DATA_OPTI_GC_CALC] = (double)NO;

  /* Switch detectors off */

  WDB[DATA_PTR_DET0] = -WDB[DATA_PTR_DET0];

  /* Prepare transport cycle */

  PrepareTransportCycle();

  /* Run criticality source simulation to get source distribution */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /* Run transportcycle(s) */

      do
        {

          PrepareCCIter();

          TransportCycle();

          IterateCC();

        }
      while(RDB[DATA_ITERATE] == (double)YES);


      /* Signal externally coupled program to end calculation */

      SignalExternal(SIGUSR2);

    }
  else
    Die(FUNCTION_NAME, "Simulation mode should be crit");

  ptr = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
  printf("%E\n", Mean(ptr, 0));

  /***************************************************************************/

  /***** Second simulation (source) ******************************************/

  /* Switch mode */

  WDB[DATA_SIMULATION_MODE] = SIMULATION_MODE_SRC;

  /* Switch detectors back on */

  WDB[DATA_PTR_DET0] = -WDB[DATA_PTR_DET0];

  /* Reset number of skip cycles */

  WDB[DATA_CRIT_SKIP] = 0.0;

  /* Set source */

  WDB[DATA_PTR_SRC0] = RDB[DATA_PTR_RIA_SRC];

  /* Prepare transport cycle */

  PrepareTransportCycle();

  /* Run simulation */

  TransportCycle();

  /***************************************************************************/

}

/*****************************************************************************/
