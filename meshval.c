/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshval.c                                      */
/*                                                                           */
/* Created:       2011/05/14 (JLe)                                           */
/* Last modified: 2017/03/22 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description:  Returns value in mesh by coordinates                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshVal:"

/*****************************************************************************/

double MeshVal(long msh, double x, double y, double z)
{
  long ptr, idx;
  double val;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Avoid compiler warning */

  val = -1.0;

  /* Get mesh index */

  if ((idx = MeshIndex(msh, x, y, z, -1.0)) < 0)
    return 0.0;

  /* Check type */

  if ((long)RDB[msh + MESH_CONTENT] == MESH_CONTENT_RES)
    {
      /* Get pointer to data */
      
      ptr = (long)RDB[msh + MESH_PTR_RES2];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      
      /* Get value */
      
      val = SumPrivateRes(ptr + idx + 1);
    }
  else if ((long)RDB[msh + MESH_CONTENT] == MESH_CONTENT_DAT)
    {
      /* Get pointer to data */
      
      ptr = (long)RDB[msh + MESH_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Get value */
      
      val = RDB[ptr + idx + 1];
    }
  else
    Die(FUNCTION_NAME, "Invalid content type");

  /* Return value */
  
  return val;
}

/*****************************************************************************/
