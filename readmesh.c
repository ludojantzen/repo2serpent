/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readmesh.c                                     */
/*                                                                           */
/* Created:       2011/05/14 (JLe)                                           */
/* Last modified: 2015/10/08 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description:  Returns value from mesh structure                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadMesh:"

/*****************************************************************************/

double ReadMesh(long msh, long i, long j, long k)
{
  long n0, n1, ptr;
  double val;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
  
  /* Get sizes */
  
  n0 = (long)RDB[msh + MESH_N0];
  n1 = (long)RDB[msh + MESH_N1];

  /* Check indexes */
  
  CheckValue(FUNCTION_NAME, "i", "", i, 0, n0 - 1);
  CheckValue(FUNCTION_NAME, "j", "", j, 0, n1 - 1);

  /* Avoid compiler warning */

  val = -1.0;

  /* Check type */

  if ((long)RDB[msh + MESH_CONTENT] == MESH_CONTENT_RES)
    {  
      /* Get pointer to data */

      ptr = (long)RDB[msh + MESH_PTR_RES2] + i + j*n0 + k*n0*n1 + 1;
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      
      /* Get value */
      
      val = GetPrivateRes(ptr);
    }
  else if ((long)RDB[msh + MESH_CONTENT] == MESH_CONTENT_DAT)
    {  
      /* Get pointer to data */

      ptr = (long)RDB[msh + MESH_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Get value */
      
      val = RDB[ptr + i + j*n0 + k*n0*n1 + 1];
    }
  else
    Die(FUNCTION_NAME, "Invalid content type");

  /* Return value */
      
  return val;
}

/*****************************************************************************/
