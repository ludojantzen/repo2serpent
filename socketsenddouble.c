/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : socketsenddouble.c                             */
/*                                                                           */
/* Created:       2019/02/04 (VVa)                                           */
/* Last modified: 2019/03/27 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sends an eight byte floating point number through socket.    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SocketSendDouble:"

/*****************************************************************************/

void SocketSendDouble(double value)
{
  SocketSend((char*)&value, sizeof(double));
}
