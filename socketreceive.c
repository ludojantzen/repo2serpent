/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : socketreceive.c                                */
/*                                                                           */
/* Created:       2019/02/04 (VVa)                                           */
/* Last modified: 2019/02/04 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Receives a known amount of bytes through socket.             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SocketReceive:"

/*****************************************************************************/

void SocketReceive(char* buffer, long length)
{
  long socket;
  long pos = 0;
  long count;

  /* Check that the connection is still on */

  if ((socket = (long)RDB[DATA_COM_SOCKET]) == 0)
    Die(FUNCTION_NAME, "Trying to receive from a closed socket.");

  /* Get data until the requested length has been obtained */

  while (pos < length)
    {
      /* Read data into current position at buffer */

      count = recv(socket, buffer+pos, length-pos, 0);

      /* Check that no error occurred */

      if (count <= 0)
        Die(FUNCTION_NAME, "Error (%d) when trying to receive.", errno);

      /* Update position in buffer */

      pos += count;
    }
}
