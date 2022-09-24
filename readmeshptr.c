/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readmeshptr.c                                  */
/*                                                                           */
/* Created:       2012/01/17 (JLe)                                           */
/* Last modified: 2012/01/17 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description:  Returns pointer from mesh structure                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadMeshPtr:"

/*****************************************************************************/

long ReadMeshPtr(long msh, long i, long j, long k)
{
  long n0, n1, n2, ptr;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Check content */

  if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_PTR)
    Die(FUNCTION_NAME, "Invalid content type");
  
  /* Get sizes */
  
  n0 = (long)RDB[msh + MESH_N0];
  n1 = (long)RDB[msh + MESH_N1];
  n2 = (long)RDB[msh + MESH_N2];

  /* Check indexes */
  
  CheckValue(FUNCTION_NAME, "i", "", i, 0, n0 - 1);
  CheckValue(FUNCTION_NAME, "j", "", j, 0, n1 - 1);
  CheckValue(FUNCTION_NAME, "k", "", k, 0, n2 - 1);

  /* Get pointer to pointer */
  
  ptr = (long)RDB[msh + MESH_PTR_PTR] + i + j*n0 + k*n0*n1;
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Return pointer */
  
  return ptr;
}

/*****************************************************************************/
