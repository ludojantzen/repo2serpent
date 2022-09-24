/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : supervisedcycle.c                              */
/*                                                                           */
/* Created:       2019/02/03 (VVa)                                           */
/* Last modified: 2019/03/29 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Handles the calculation in the case of a supervisor code.    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "krakensignals.h"

#define FUNCTION_NAME "SupervisedCycle:"

/*****************************************************************************/

void SupervisedCycle()
{
  long socket, signal;

  /* Return if not running coupled calculation */

  if(RDB[DATA_RUN_CC] == (double)NO)
    return;

  /* Return if not communicating via sockets */

  if ((socket = (long)RDB[DATA_COM_SOCKET]) == 0)
    return;

  fprintf(outp, "Starting supervised coupled calculation.\n\n");

  signal = -1;

  while (signal != KRAKEN_TERMINATE)
    {

    /* Get next signal */

    signal = GetSignalFromSupervisor();

    switch (signal)
      {
      case KRAKEN_TERMINATE:
        fprintf(outp, "Received signal to terminate calculation.\n");

        break;
      case KRAKEN_GIVE_INPUT_FIELD_TEMPLATES:
        fprintf(outp, "Received signal to give input field templates.\n");

        /* Now simply send zero as the number of accepted fields */

        SendIFCInputTemplates();

        break;
      case KRAKEN_GIVE_OUTPUT_FIELD_TEMPLATES:
        fprintf(outp, "Received signal to give output field templates.\n");

        /* Now simply send zero as the number of provided fields */

        SendIFCOutputTemplates();

        break;
      case KRAKEN_SOLVE_STEADY_STATE:

        fprintf(outp, "Received signal to solve steady state.\n");

        /* This should be called only before the first iteration */

        if (RDB[DATA_SOL_REL_ITER] == 0)
          PrepareTransportCycle();

        /* Prepare coupled calculation iteration */

        PrepareCCIter();

        /* Run transport cycle */

        TransportCycle();

        /* Collect data and relax results */

        RelaxCCResults();

        /* Increment iteration number */

        WDB[DATA_SOL_REL_ITER] = RDB[DATA_SOL_REL_ITER] + 1.0;

        break;
      case KRAKEN_GIVE_FIELD_DATA:

        fprintf(outp, "Received signal to send output data.\n");

        SendIFCOutputData();

        break;
      case KRAKEN_TAKE_FIELD_DATA:

        fprintf(outp, "Received signal to receive input data.\n");

        ReceiveIFCInputData();

        /* Sample interface data if requested */

        SampleIFCData(NO);

        break;
      default:
        fprintf(outp, "Signal %ld not implemented for supervised cycle.\n", signal);
      }
    }
 }

/*****************************************************************************/
