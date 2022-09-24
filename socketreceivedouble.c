/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : socketreceivedouble.c                          */
/*                                                                           */
/* Created:       2019/02/04 (VVa)                                           */
/* Last modified: 2019/03/27 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Receives an eight byte floating point number through socket. */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SocketReceiveDouble:"

/*****************************************************************************/

void SocketReceiveDouble(double *value)
{
  SocketReceive((char*)value, sizeof(double));
}
