/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : signalhandler.c                                */
/*                                                                           */
/* Created:       2013/04/12 (VVa)                                           */
/* Last modified: 2014/11/05 (VVa)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Catches signals and handles them accordingly                 */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SignalHandler:"

/*****************************************************************************/

void SignalHandler(int signo)
{

  /* Change waiting flag if SIGUSR1 is received */
  /* = external program has updated interfaces  */

  if (signo == SIGUSR1)
    WDB[DATA_WAITING] = (double)NO;

  /* Change iteration flag if SIGUSR2 is received               */
  /* = external program has deemed convergence to be sufficient */
  
  if (signo == SIGUSR2)
    {
      WDB[DATA_WAITING] = (double)NO;

      WDB[DATA_ITERATE] = (double)NO;
    }

}

/*****************************************************************************/
