/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshtot.c                                      */
/*                                                                           */
/* Created:       2011/05/14 (JLe)                                           */
/* Last modified: 2015/10/04 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description:  Returns sum of all mesh values                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshTot:"

/*****************************************************************************/

double MeshTot(long msh)
{
  long ptr;
  double tot;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Avoid compiler warning */

  tot = -1.0;

  /* Check type */

  if ((long)RDB[msh + MESH_CONTENT] == MESH_CONTENT_RES)
    {
      /* Get pointer to data */
      
      ptr = (long)RDB[msh + MESH_PTR_RES2];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      
      /* Get value */

      tot = SumPrivateRes(ptr);
    }
  else if ((long)RDB[msh + MESH_CONTENT] == MESH_CONTENT_DAT)
    {
      /* Get pointer to data */
      
      ptr = (long)RDB[msh + MESH_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Get value */
      
      tot = RDB[ptr];
    }
  else
    Die(FUNCTION_NAME, "Invalid content type");

  /* Return value */

  return tot;
}

/*****************************************************************************/
