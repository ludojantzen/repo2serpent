/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : socketreceivestring.c                          */
/*                                                                           */
/* Created:       2019/02/04 (VVa)                                           */
/* Last modified: 2019/03/27 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Receives a string of unknown length through socket.          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SocketReceiveString:"

/*****************************************************************************/

void SocketReceiveString(char *buffer)
{
  int64_t str_length;

  /* Get the length of the string to read */

  SocketReceiveLong(&str_length);

  /* Check that the string will fit into our buffer */

  if (str_length >= MAX_STR)
    Die("recvString", "The string to be received would be %ld characters, which "
        "exceeds the buffer of %ld characters.", str_length, MAX_STR);

  /* Read the string */

  SocketReceive(buffer, str_length*sizeof(char));

  /* Null terminate it */

  buffer[str_length] = '\0';
}
