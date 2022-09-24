/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : socketsendstring.c                             */
/*                                                                           */
/* Created:       2019/02/04 (VVa)                                           */
/* Last modified: 2019/03/27 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sends a string of unknown length through socket.             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SocketSendString:"

/*****************************************************************************/

void SocketSendString(const char *buffer)
{
  long str_length;

  /* Count string length */

  str_length = strlen(buffer);

  /* First send the length of the string */

  SocketSendLong(str_length);

  /* Then send the data itself */

  SocketSend(buffer, str_length*sizeof(char));
}
