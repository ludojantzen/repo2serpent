/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fissmtxindex.c                                 */
/*                                                                           */
/* Created:       2013/01/17 (JLe)                                           */
/* Last modified: 2017/03/22 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Returns index to fission matrix                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FissMtxIndex:"

/*****************************************************************************/

long FissMtxIndex(long mat, long id)
{
  long type, idx, ptr, uni, lvl, fmx, msh;
  double x, y, z;

  /* Get matrix type */

  if ((type = (long)RDB[DATA_FMTX_TYPE]) < 1)
    return -1;

  /* Avoid compiler warning */

  idx = -1;

  /* Check type */

  if (type == FISSION_MATRIX_TYPE_MAT)
    {
      /* Check material pointer */

      if (mat < VALID_PTR)
        idx = -1;
      else
        idx = (long)RDB[mat + MATERIAL_FMTX_IDX];
    }
  else if (type == FISSION_MATRIX_TYPE_UNI)
    {
      /* Get collision universe */
      
      ptr = (long)RDB[DATA_PTR_COLLISION_UNI];
      uni = (long)GetPrivateData(ptr, id);
          
      /* Check pointer */
  
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);
      
      /* Get index */
          
      idx = (long)RDB[uni + UNIVERSE_FMTX_IDX];
    }
  else if (type == FISSION_MATRIX_TYPE_XYZ)
    {
      /* Pointer to first level */
      
      lvl = (long)RDB[DATA_PTR_LVL0];
      CheckPointer(FUNCTION_NAME, "(lvl)", DATA_ARRAY, lvl);

      /* Pointer to private data */

      lvl = (long)RDB[lvl + LVL_PTR_PRIVATE_DATA];
      CheckPointer(FUNCTION_NAME, "(lvl)", PRIVA_ARRAY, lvl);

      /* Get coordinates */

      x = GetPrivateData(lvl + LVL_PRIV_X, id);
      y = GetPrivateData(lvl + LVL_PRIV_Y, id);
      z = GetPrivateData(lvl + LVL_PRIV_Z, id);

      /* Get pointer to matrix */

      fmx = (long)RDB[DATA_PTR_FMTX];
      CheckPointer(FUNCTION_NAME, "(fmx)", DATA_ARRAY, fmx);

      /* Get pointer to mesh */

      msh = (long)RDB[fmx + FMTX_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get index */

      idx = MeshIndex(msh, x, y, z, -1.0);
    }
  else
    Die(FUNCTION_NAME, "Invalid matrix type");

  /* Return index */

  return idx;  
}

/*****************************************************************************/
