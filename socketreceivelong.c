/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : socketreceivelong.c                            */
/*                                                                           */
/* Created:       2019/02/04 (VVa)                                           */
/* Last modified: 2019/03/27 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Receives an eight byte integer through socket.               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SocketReceiveLong:"

/*****************************************************************************/

void SocketReceiveLong(long *value)
{
  SocketReceive((char*)value, sizeof(int64_t));
}
