/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addstlpoint.c                                  */
/*                                                                           */
/* Created:       2014/03/10 (JLe)                                           */
/* Last modified: 2017/03/22 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Tests for coincident points and adds point to STL structure  */
/*                                                                           */
/* Comments: - NOTE: Toi silmukka naapuricellien yli ei oo ihan              */
/*             optimaalinen.                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddSTLPoint()"

/*****************************************************************************/

long AddSTLPoint(long ***arr, long stl, long fle, long sld, 
                 double x, double y, double z)
{
  long loc0, ptr, i, j, k, N, i0, j0, k0;
  double dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax, dr;

  /* Get boundaries */

  xmin = RDB[stl + STL_XMIN];
  xmax = RDB[stl + STL_XMAX];
  ymin = RDB[stl + STL_YMIN];
  ymax = RDB[stl + STL_YMAX];
  zmin = RDB[stl + STL_ZMIN];
  zmax = RDB[stl + STL_ZMAX];

  /* Scale and shift origin */

  x = x*RDB[fle + STL_FILE_SCALING] - RDB[fle + STL_FILE_X0];
  y = y*RDB[fle + STL_FILE_SCALING] - RDB[fle + STL_FILE_Y0];
  z = z*RDB[fle + STL_FILE_SCALING] - RDB[fle + STL_FILE_Z0];

  /* Get temporary lattice size */

  N = (long)RDB[DATA_STL_TEMP_ARRAY_SIZE];

  /* Check points */

  CheckValue(FUNCTION_NAME, "x", "", x, xmin, xmax);
  CheckValue(FUNCTION_NAME, "y", "", y, ymin, ymax);
  CheckValue(FUNCTION_NAME, "z", "", z, zmin, zmax);
  
  /* Check if points are merged */

  if ((dr = RDB[stl + STL_MERGE_RAD]) <= 0.0)
    {
      /* Not merged, addjust boundaries */

      xmin = xmin - 1E-6;
      xmax = xmax + 1E-6;
      ymin = ymin - 1E-6;
      ymax = ymax + 1E-6;
      zmin = zmin - 1E-6;
      zmax = zmax + 1E-6;
    }
  else
    {
      /* Check that radius is reasonable */
      
      CheckValue(FUNCTION_NAME, "dr", "", dr, 0.0, 1.0);

      /* Adjust boundaries */

      xmin = xmin - dr - 1E-6;
      xmax = xmax + dr + 1E-6;
      ymin = ymin - dr - 1E-6;
      ymax = ymax + dr + 1E-6;
      zmin = zmin - dr - 1E-6;
      zmax = zmax + dr + 1E-6;

      /* Calculate indexes */
      
      i0 = (long)((double)N*(x - xmin)/(xmax - xmin));
      j0 = (long)((double)N*(y - ymin)/(ymax - ymin));
      k0 = (long)((double)N*(z - zmin)/(zmax - zmin));
      
      /* Loop over cells in the neighbourhood */
      
      for (i = i0 - 1; i < i0 + 2; i++)
        for (j = j0 - 1; j < j0 + 2; j++)
          for (k = k0 - 1; k < k0 + 2; k++)
            {
              /* Check boundary values */
              
              if ((i < 0) || (i > N - 1) || 
                  (j < 0) || (j > N - 1) ||
                  (k < 0) || (k > N - 1))
                continue;
              
              /* Get pointer and loop over contents */
              
              ptr = arr[i][j][k];
              while (ptr > VALID_PTR)
                {
                  /* Calculate distance and check */
                  
                  if ((dx = fabs(x - RDB[ptr + STL_POINT_X])) < dr)
                    if ((dy = fabs(y - RDB[ptr + STL_POINT_Y])) < dr)
                      if ((dz = fabs(z - RDB[ptr + STL_POINT_Z])) < dr)
                        if (dx*dx + dy*dy + dz*dz < dr*dr)
                          {
                            /* Combine point */
                            
                            return ptr;
                          }
                  
                  /* Next */
                  
                  ptr = NextItem(ptr);
                }
            }
    }

  /* Not merged with previous, allocate memory */

  ptr = ReallocMem(DATA_ARRAY, STL_POINT_BLOCK_SIZE);

  /* Put coordinates */

  WDB[ptr + STL_POINT_X] = x;
  WDB[ptr + STL_POINT_Y] = y;
  WDB[ptr + STL_POINT_Z] = z;

  /* Calculate index */

  i = (long)((double)N*(x - xmin)/(xmax - xmin));
  j = (long)((double)N*(y - ymin)/(ymax - ymin));
  k = (long)((double)N*(z - zmin)/(zmax - zmin));
  
  /* Update list */

  loc0 = arr[i][j][k];
  WDB[ptr + LIST_PTR_NEXT] = (double)loc0;
  arr[i][j][k] = ptr;

  /* Add counter */

  WDB[sld + STL_SOLID_N_POINTS] = RDB[sld + STL_SOLID_N_POINTS] + 1.0;

  /* Return pointer */

  return ptr;
}

/*****************************************************************************/
