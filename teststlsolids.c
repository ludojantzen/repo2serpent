/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : teststlgeometry.c                              */
/*                                                                           */
/* Created:       2014/12/10 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Tests STL solids separately by random sampling               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TestSTLSolids:"

/*****************************************************************************/

void TestSTLSolids(long stl, long np, long nd)
{
  long sld, msh, ok0, ok, mode, n, m, nf, id;
  unsigned long seed;
  double xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, u, v, w;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(stl)", DATA_ARRAY, stl);

  fprintf(outp, "\nTesting STL solids...\n\n");

  /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */

  ExpandPrivateArrays();

  /* Set test mode */

  WDB[DATA_STL_GEOM_TEST_MODE] = (double)YES;

  /* Set mode */

  mode = STL_SEARCH_MODE_FAST;

  /* Get pointer to search mesh */
  
  msh = (long)RDB[stl + STL_PTR_FACET_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
  
  /* Loop over solids */
  
  sld = (long)RDB[stl + STL_PTR_SOLIDS];
  while (sld > VALID_PTR)
    {
      /* Get boundaries */
      
      xmin = RDB[sld + STL_SOLID_XMIN];
      xmax = RDB[sld + STL_SOLID_XMAX];
      ymin = RDB[sld + STL_SOLID_YMIN];
      ymax = RDB[sld + STL_SOLID_YMAX];
      zmin = RDB[sld + STL_SOLID_ZMIN];
      zmax = RDB[sld + STL_SOLID_ZMAX];

      /* Reset fail counter */
      
      nf = 0;
      
      /* Loop over points */

      for (n = 0; n < np; n++)
        {
          /* Get OpenMP thread num */

          id = 0;
            
          /* Init random number sequence */
          
          seed = ReInitRNG(n + 1);
          SEED[id*RNG_SZ] = seed;
          
          /* Sample point */
          
          x = RandF(id)*(xmax - xmin) + xmin;
          y = RandF(id)*(ymax - ymin) + ymin;
          z = RandF(id)*(zmax - zmin) + zmin;

          /* Reset check */

          ok0 = -1;
          
          /* Loop over directions */
          
          for (m = 0; m < nd; m++)
            {
              /* Sample direction */
              
              IsotropicDirection(&u, &v, &w, id);
              
              /* Perform ray test */
              
              ok = STLRayTest(sld, msh, x, y, z, u, v, w, mode, id);

              /* Check with previous */
        
              if (ok < 0)
                continue;
              if (ok0 == -1)
                ok0 = ok;
              else if (ok != ok0)
                {
                  /* Add to fail counter */
                  
                  nf++;
                  
                  /* Break loop */

                  break;
                }
            }
        }
      

        printf("FAIL: %s %ld\n", GetText(sld + STL_SOLID_PTR_FNAME), nf);
      
      /* Next solid */
      
      sld = NextItem(sld);
    }
}

/*****************************************************************************/
