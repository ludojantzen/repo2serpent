/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getsignalfromsupervisor.c                      */
/*                                                                           */
/* Created:       2019/02/03 (VVa)                                           */
/* Last modified: 2019/02/04 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Waits for a signal from a supervisor code (e.g. Cerberus)    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetSignalFromSupervisor:"

/*****************************************************************************/

long GetSignalFromSupervisor()
{
  long socket, value;

  /* Return if not running coupled calculation */

  if(RDB[DATA_RUN_CC] == (double)NO)
    return 0;

  /* Return if not communicating via sockets */

  if ((socket = (long)RDB[DATA_COM_SOCKET]) == 0)
    return 0;

  /* Receive a single integer */

  SocketReceiveLong(&value);

  return value;
 }

/*****************************************************************************/
