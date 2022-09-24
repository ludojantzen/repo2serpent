/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : initsignal.c                                   */
/*                                                                           */
/* Created:       2013/04/12 (VVa)                                           */
/* Last modified: 2014/11/21 (VVa)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Initializes the signal handler and signal data               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitSignal:"

/*****************************************************************************/

void InitSignal()
{
  struct sigaction wakeup;
     
  /* Wakeup signals will be handled by SignalHandler function */

  wakeup.sa_handler = SignalHandler;

  /* Initialize the wakeup sigaction */

  sigemptyset(&wakeup.sa_mask);

  /* No special flags*/

  wakeup.sa_flags = 0;
  
  /* Associate SIGUSR1 with the sigaction wakeup */

  sigaction (SIGUSR1, &wakeup, NULL);

  /* Associate SIGUSR2 with the sigaction wakeup  */
  /* Will be handled differently in SignalHandler */

  sigaction (SIGUSR2, &wakeup, NULL);

}

/*****************************************************************************/
