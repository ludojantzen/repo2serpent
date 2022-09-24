/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addmeshidx.c                                   */
/*                                                                           */
/* Created:       2014/05/23 (JLe)                                           */
/* Last modified: 2015/10/04 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description:  Adds value in mesh structure using indexes                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddMeshIdx:"

/*****************************************************************************/

void AddMeshIdx(long msh, double val, long i, long j, long k, long id)
{
  long ptr, idx, n0, n1, n2;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Check content */

  if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_RES)
    Die(FUNCTION_NAME, "Invalid content type");

  /* Get sizes */
  
  n0 = (long)RDB[msh + MESH_N0];
  n1 = (long)RDB[msh + MESH_N1];
  n2 = (long)RDB[msh + MESH_N2];

  /* Check */

  if ((i < 0) || (i > n0 - 1))
    return;
  else if ((j < 0) || (j > n1 - 1))
    return;
  else if ((k < 0) || (k > n2 - 1))
    return;

  /* Calculate index */

  idx = i + j*n0 + k*n0*n1;

  /* Get pointer to data */
  
  ptr = (long)RDB[msh + MESH_PTR_RES2];
  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

  /* Add value */

  AddPrivateRes(ptr + idx + 1, val, id);

  /* Add to total */

  AddPrivateRes(ptr, val, id);
}

/*****************************************************************************/
