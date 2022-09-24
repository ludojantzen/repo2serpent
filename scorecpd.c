/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorecpd.c                                     */
/*                                                                           */
/* Created:       2011/05/15 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Scores core power distribution                               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreCPD:"

/*****************************************************************************/

void ScoreCPD(double fissE, double wgt, double z, long id)
{
  long ptr, lat, n0, n1, n2, nz, ncol;
  double zmin, zmax;

  /* Check depth */

  if (RDB[DATA_CORE_PDE_DEPTH] < 1.0)
    return;

  /* Reset indexes */

  n0 = -1;
  n1 = -1;
  n2 = -1;

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Core level */

  ptr = (long)RDB[DATA_CORE_PDE_PTR_CORE];
  while (ptr > VALID_PTR)
    {
      /* Pointer to lattice */

      lat = (long)RDB[ptr + CPD_PTR_LAT];

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(lat)", DATA_ARRAY, lat);

      /* Check for collision */

      if ((n0 = (long)TestValuePair(lat + LAT_COL_CELL_IDX, (double)ncol, 
                                    id)) > -1)
        break;
      
      /* Next */
      
      ptr = NextItem(ptr);
    }

  /* Assembly level */

  ptr = (long)RDB[DATA_CORE_PDE_PTR_ASS];
  while (ptr > VALID_PTR)
    {
      /* Pointer to lattice */

      lat = (long)RDB[ptr + CPD_PTR_LAT];

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(lat)", DATA_ARRAY, lat);

      /* Check for collision */

      if ((n1 = (long)TestValuePair(lat + LAT_COL_CELL_IDX, (double)ncol, 
                                    id)) > -1)
        break;
      
      /* Next */
      
      ptr = NextItem(ptr);
    }

  /* Axial bin */

  zmin = RDB[DATA_CORE_PDE_ZMIN];
  zmax = RDB[DATA_CORE_PDE_ZMAX];
  nz = (long)RDB[DATA_CORE_PDE_NZ];

  /* Check z-coordinate */
      
  if (nz == 1)
    n2 = -1;
  else if ((z >= zmin) && (z < zmax))
    {
      /* Calculate bin width */

      n2 = (long)((z - zmin)/(zmax - zmin)*((double)nz));
    }

  /* Score */

  if (n0 > -1)
    {
      CheckValue(FUNCTION_NAME, "n0", "", n0, 0, 
                 (long)RDB[DATA_CORE_PDE_N0] - 1);

      ptr = (long)RDB[DATA_CORE_PDE_PTR_RES0];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(fissE, wgt, ptr, id, n0);

      if (n1 > -1)
        {
          CheckValue(FUNCTION_NAME, "n1", "", n1, 0, 
                     (long)RDB[DATA_CORE_PDE_N1] - 1);
          
          ptr = (long)RDB[DATA_CORE_PDE_PTR_RES1];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(fissE, wgt, ptr, id, -1, n0, n1);
          
          if (n2 > -1)
            {
              CheckValue(FUNCTION_NAME, "n2", "", n2, 0, 
                         (long)RDB[DATA_CORE_PDE_N2] - 1);
              
              ptr = (long)RDB[DATA_CORE_PDE_PTR_RES2];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(fissE, wgt, ptr, id, -1, n0, n1, n2);
            }
        }
      else if (n2 > -1)
        {
          CheckValue(FUNCTION_NAME, "n2", "", n2, 0, 
                     (long)RDB[DATA_CORE_PDE_N2] - 1);
          
          ptr = (long)RDB[DATA_CORE_PDE_PTR_RES2];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(fissE, wgt, ptr, id, -1, n0, 0, n2);
        }
    }
}

/*****************************************************************************/
