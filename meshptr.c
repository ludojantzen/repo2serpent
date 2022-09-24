/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshptr.c                                      */
/*                                                                           */
/* Created:       2012/01/17 (JLe)                                           */
/* Last modified: 2017/03/22 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description:  Returns pointer in mesh by coordinates                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshPtr:"

/*****************************************************************************/

long MeshPtr(long msh, double x, double y, double z)
{
  long ptr, loc0, idx;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Check content */

  if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_PTR)
    Die(FUNCTION_NAME, "Invalid content type");

  /* Get mesh index */

  if ((idx = MeshIndex(msh, x, y, z, -1.0)) < 0)
    return -1;

  /* Get direct pointer to data */
  
  ptr = (long)RDB[msh + MESH_PTR_PTR] + idx;
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check type */

  if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_ADAPTIVE)
    {
      /* Adaptive mesh, check pointer */

      if ((loc0 = (long)RDB[ptr]) == NULLPTR)
        {
          /* Tää voi olla ongelma */

          return -1;
        }
      else if (loc0 < -VALID_PTR)
        {
          /* Pointer to new mesh, call recursively */

          ptr = MeshPtr(-loc0, x, y, z);
        }
    }

  /* Return pointer */
  
  return ptr;
}

/*****************************************************************************/
