/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sendifcmesh.c                                  */
/*                                                                           */
/* Created:       2019/02/05 (VVa)                                           */
/* Last modified: 2019/03/29 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sends interface mesh to Cerberus     .                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "krakenmeshtypes.h"

#define FUNCTION_NAME "SendIFCMesh:"

void SendIFCMesh(long msh)
{
  long mshtype, nx, ny, nz;
  double minx, miny, minz;
  double maxx, maxy, maxz;

  /* Get mesh type */

  mshtype = RDB[msh + MESH_TYPE];

  if (mshtype == MESH_TYPE_CARTESIAN)
    {
      /* Get mesh data */

      nx = RDB[msh + MESH_N0];
      ny = RDB[msh + MESH_N1];
      nz = RDB[msh + MESH_N2];

      minx = 0.01*RDB[msh + MESH_MIN0];
      maxx = 0.01*RDB[msh + MESH_MAX0];
      miny = 0.01*RDB[msh + MESH_MIN1];
      maxy = 0.01*RDB[msh + MESH_MAX1];
      minz = 0.01*RDB[msh + MESH_MIN2];
      maxz = 0.01*RDB[msh + MESH_MAX2];

      /* Send type */

      SocketSendLong(KRAKEN_MESH_TYPE_STRUCTURED_CARTESIAN);

      /* Send data */

      SocketSendLong(nx);
      SocketSendLong(ny);
      SocketSendLong(nz);

      SocketSendDouble(minx);
      SocketSendDouble(maxx);
      SocketSendDouble(miny);
      SocketSendDouble(maxy);
      SocketSendDouble(minz);
      SocketSendDouble(maxz);
    }
  else
    {
      Die(FUNCTION_NAME,
          "IFC output mesh sender not implemented for regular mesh type %ld",
          mshtype);
    }
}
