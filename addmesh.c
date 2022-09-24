/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addmesh.c                                      */
/*                                                                           */
/* Created:       2011/05/14 (JLe)                                           */
/* Last modified: 2017/03/22 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description:  Adds value in mesh structure                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddMesh:"

/*****************************************************************************/

void AddMesh(long msh, double val, double x, double y, double z, long id)
{
  long ptr, idx;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Check content */

  if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_RES)
    Die(FUNCTION_NAME, "Invalid content type");

  /* Get mesh index */

  if ((idx = MeshIndex(msh, x, y, z, -1.0)) < 0)
    return;

  /* Get pointer to data */
  
  ptr = (long)RDB[msh + MESH_PTR_RES2];
  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

  /* Add value */

  AddPrivateRes(ptr + idx + 1, val, id);

  /* Add to total */

  AddPrivateRes(ptr, val, id);
}

/*****************************************************************************/
