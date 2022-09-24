/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stlraytest.c                                   */
/*                                                                           */
/* Created:       2014/11/24 (JLe)                                           */
/* Last modified: 2016/06/15 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description:  Performs ray test needed for cell search in STL geometry    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "STLRayTest:"

/*****************************************************************************/

long STLRayTest(long sld, long msh, double x, double y, double z, 
                double u, double v, double w, long mode, long id)
{
  long loc0, loc1, ptr, i, cross, fail;
  double xmin, xmax, ymin, ymax, zmin, zmax, l, d, b, u0, v0, w0, dir, dir0;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(sld)", DATA_ARRAY, sld);
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Check coordinates and direction cosines */

  CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);

  CheckValue(FUNCTION_NAME, "u", "", u, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "v", "", v, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "w", "", w, -1.0, 1.0);

  /* Check if ray is too parallel to search mesh cell boundaries */

  if ((fabs(u) < 0.001) || (fabs(v) < 0.001) || (fabs(w) < 0.001))
    return STL_RAY_TEST_FAIL_PARA;

  /* Get bounding box */

  xmin = RDB[sld + STL_SOLID_XMIN];
  xmax = RDB[sld + STL_SOLID_XMAX];
  ymin = RDB[sld + STL_SOLID_YMIN];
  ymax = RDB[sld + STL_SOLID_YMAX];
  zmin = RDB[sld + STL_SOLID_ZMIN];
  zmax = RDB[sld + STL_SOLID_ZMAX];

  /* Check values */
  
  CheckValue(FUNCTION_NAME, "xmin", "", xmin, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "xmax", "", xmax, xmin, INFTY);

  CheckValue(FUNCTION_NAME, "ymin", "", ymin, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "ymax", "", ymax, ymin, INFTY);

  CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);

  /* Reset number of crossings */

  cross = 0;

  /* Main loop */

  for (i = 0; i < 10000; i++)
    {  
      /* Check if point is outside the bounding box */

      if ((x < xmin) || (x > xmax) || 
          (y < ymin) || (y > ymax) || 
          (z < zmin) || (z > zmax))
        break;
        
      /* Distance to search mesh boundaries */

      b = NearestMeshBoundary(msh, x, y, z, u, v, w, &fail);
      CheckValue(FUNCTION_NAME, "b", "", b, -ZERO, INFTY);

      /* Check failure */

      if (fail == YES)
        return STL_RAY_TEST_FAIL_EDGE;

      /* Reset minimum distance to facets and direction */

      d = INFTY;
      dir0 = INFTY;

      /* Get pointer to search mesh cell */

      if ((loc0 = MeshPtr(msh, x, y, z)) < VALID_PTR)
        return NO;
      
      /* Pointer to content */

      if ((loc0 = (long)RDB[loc0]) > VALID_PTR)
        {
          /* Check preassigned information */

          if ((loc1 = (long)RDB[loc0 + SEARCH_MESH_CELL_CONTENT]) < VALID_PTR)
            {
              /* Check mode */

              if (mode == STL_SEARCH_MODE_FAST)
                {
                  /* Compare pointer */
                  
                  if (sld == -loc1)
                    return YES;
                  else if (loc1 == 0)
                    return NO;
                  else
                    {
                      /* Reset pointer to prevent going to next loop */

                      loc0 = -1;
                    }
                }
              else if (mode == STL_SEARCH_MODE_SAFE)
                {
                  /* Compare pointer */
                  
                  if (sld == -loc1)
                    {
                      /* Check number of crossings */
                      
                      if (cross % 2)
                        return NO;
                      else
                        return YES;
                    }
                  else if (loc1 == 0)
                    {
                      /* Check number of crossings */
                      
                      if (cross % 2)
                        return YES;
                      else
                        return NO;
                    }
                  else
                    {
                      /* Reset pointer to prevent going to next loop */

                      loc0 = -1;
                    }
                }
              else
                Die(FUNCTION_NAME, "Invalid mode");
            }
          
          /* Loop over content list */

          while (loc0 > VALID_PTR)
            {
              /* Pointer to content */
              
              loc1 = (long)RDB[loc0 + SEARCH_MESH_CELL_CONTENT];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
              
              /* Get pointer to solid */

              ptr = (long)RDB[loc1 + STL_FACET_PTR_SOLID];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              
              /* Check pointer */

              if (ptr == sld)
                {
                  /* Get surface normal */
                  
                  u0 = RDB[loc1 + STL_FACET_NORM_U];
                  v0 = RDB[loc1 + STL_FACET_NORM_V];
                  w0 = RDB[loc1 + STL_FACET_NORM_W];

                  /* Check values */

                  CheckValue(FUNCTION_NAME, "u0", "", u0, -1.0, 1.0);
                  CheckValue(FUNCTION_NAME, "v0", "", v0, -1.0, 1.0);
                  CheckValue(FUNCTION_NAME, "w0", "", w0, -1.0, 1.0);
                  
                  /* Calculate cosine between direction and normal */
                  
                  dir = u*u0 + v*v0 + w*w0;
                  
                  /* Check failure */
                  
                  if (fabs(dir) < 0.001)
                    {
                      /* Ray is too parallel to facet */

                      return STL_RAY_TEST_FAIL_PARA;
                    }
                  
                  /* Calculate distance to facet */
                  
                  if ((l = STLFacetDistance(loc1, x, y, z, u, v, w, YES, id)) 
                      < 0.0)
                    {
                      /* Intersection point is too close to edge */
                      
                      return STL_RAY_TEST_FAIL_EDGE;
                    }
                  else if (fabs(l - b) < 1E-6)
                    {
                      /* Facet is too close to search mesh cell boundary */

                      return STL_RAY_TEST_FAIL_MESH;
                    }
                  else if (mode == STL_SEARCH_MODE_SAFE)
                    {
                      /* Add to number of crossings */
                      
                      if (l < b)
                        cross++;
                    }
                  else if ((d != INFTY) && (l == d))
                    {
                      /* Two facets overlap */
                      
                      return STL_FACET_OVERLAP;
                    }
                  else if (l < d)
                    {
                      /* Remember distance and direction */
                      
                      d = l;
                      dir0 = dir;
                    }
                }
              
              /* Next */
              
              loc0 = NextItem(loc0);
            }
        }

      /* Check variables */

      CheckValue(FUNCTION_NAME, "d", "", d, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "b", "", b, -ZERO, INFTY);

      /* Check coordinates */
      
      CheckValue(FUNCTION_NAME, "x", "", x, xmin, xmax);
      CheckValue(FUNCTION_NAME, "y", "", y, ymin, ymax);
      CheckValue(FUNCTION_NAME, "z", "", z, zmin, zmax);

      /* Compare distances */

      if (d < b)
        {
          /* Check direction */

          CheckValue(FUNCTION_NAME, "dir0", "", dir0, -1.0, 1.0);

          /* Ray intersects a facet, check direction */

          if (dir0 < 0.0)
            return NO;
          else
            return YES;
        }
      else
        {
          /* Move over mesh cell boundary (toi 1E-9 toimii paremmin */
          /* kuin EXTRAP_L = 1E-6) */

          x = x + u*(b + 1E-9);
          y = y + v*(b + 1E-9);
          z = z + w*(b + 1E-9);
        }
    }

  /* Check infinite loop */

  if (i == 10000)
    {
      /* Test failed */

      return STL_RAY_TEST_FAIL_LOOP;
    }

  /* Check mode */

  if (mode == STL_SEARCH_MODE_SAFE)
    {
      /* Check number of crossings */

      if (cross % 2)
        return YES;
      else
        return NO;
    }
  else
    {
      /* Must be outside */

      return NO;
    }  
}

/*****************************************************************************/
