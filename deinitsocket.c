/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : deinitsocket.c                                 */
/*                                                                           */
/* Created:       2019/02/04 (VVa)                                           */
/* Last modified: 2019/02/04 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Deinitializes the socket used for communications.            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DeinitSocket:"

void DeinitSocket()
{
  long socket;

  /* Check that the connection is still on */

  if ((socket = (long)RDB[DATA_COM_SOCKET]) <= 0)
    return;

  /* Close the socket */

  close(socket);

  /* Zero-out the socket address */

  WDB[DATA_COM_SOCKET] = 0;
}
