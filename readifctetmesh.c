/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifctetmesh.c                               */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2018/11/08 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads tet-mesh multi-physics interfaces                      */
/*                                                                           */
/* Comments:   -Split from readinterface.c for 2.1.22                        */
/*             -Deprecated in favor of OpenFOAM interfaces.                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCTetMesh:"

/*****************************************************************************/

void ReadIFCTetMesh(long loc0, long update)
{
  Die(FUNCTION_NAME, "Interface type %ld is deprecated", (long)RDB[loc0 + IFC_TYPE]);
}

/*****************************************************************************/
