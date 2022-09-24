/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : teststlgeometry.c                              */
/*                                                                           */
/* Created:       2014/03/07 (JLe)                                           */
/* Last modified: 2020/06/23 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Tests STL based geometries by random sampling                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TestSTLGeometry:"

/*****************************************************************************/

void TestSTLGeometry()
{
  long stl, loc0, loc1, loc2, ptr, n, m, np, nd, nf, id;
  unsigned long seed;
  double xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, u, v, w;

  /* Check pointer */

  if ((stl = (long)RDB[DATA_PTR_STL0]) < VALID_PTR)
    return;

  /* Get number of points */

  if ((np = (long)RDB[DATA_STL_TEST_N_PTS]) < 1)
    return;

  fprintf(outp, "Testing STL geometries...\n\n");

  /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */

  ExpandPrivateArrays();

  /* Get number of directions */

  nd = (long)RDB[DATA_STL_TEST_N_DIR];

  /* Loop over geometries */

  while (stl > VALID_PTR)
    {
      /* Get boundaries */

      xmin = RDB[stl + STL_XMIN];
      xmax = RDB[stl + STL_XMAX];
      ymin = RDB[stl + STL_YMIN];
      ymax = RDB[stl + STL_YMAX];
      zmin = RDB[stl + STL_ZMIN];
      zmax = RDB[stl + STL_ZMAX];

      /* Reset fail counter */

      nf = 0;

      fprintf(outp, "Testing universe %s by sampling %ld directions in %ld points:\n\n", GetText(stl + STL_PTR_NAME), nd, np);

#ifdef OPEN_MP
#pragma omp parallel private (n, m, id, seed, x, y, z, u, v, w, ptr, loc0, loc1, loc2)
#endif
      {

#ifdef OPEN_MP
#pragma omp for
#endif
        /* Loop over points */

        for (n = 0; n < np; n++)
          {
            /* Get OpenMP thread num */

            id = OMP_THREAD_NUM;

            /* Init random number sequence */

            seed = ReInitRNG(n + 1);
            SEED[id*RNG_SZ] = seed;

            /* Sample point */

            x = RandF(id)*(xmax - xmin) + xmin;
            y = RandF(id)*(ymax - ymin) + ymin;
            z = RandF(id)*(zmax - zmin) + zmin;

            /* Avoid compiler warning */

            loc0 = -1;

            /* Loop over directions */

            for (m = 0; m < nd; m++)
              {
                /* Sample direction */

                IsotropicDirection(&u, &v, &w, id);

                /* Find solid */

                ptr = FindSTLSolid(stl, x, y, z, u, v, w, YES, id);

                /* Check match */

                if (m == 0)
                  loc0 = ptr;
                else if (loc0 != ptr)
                  {
                    /* Check error type */
#ifdef OPEN_MP
#pragma omp critical
#endif
                    {
                      /*
                      if ((loc0 < VALID_PTR) && (ptr > VALID_PTR))
                        fprintf(outp, "- Point [%7.2f, %7.2f, %7.2f] is outside and inside solid \"%s\" (n = %ld)\n", x, y, z, GetText(ptr + STL_SOLID_PTR_FNAME), m);
                      else if ((loc0 > VALID_PTR) && (ptr < VALID_PTR))
                        fprintf(outp, "- Point [%7.2f, %7.2f, %7.2f] is outside and inside solid \"%s\" (n = %ld)\n", x, y, z, GetText(loc0 + STL_SOLID_PTR_FNAME), m);
                      else
                        fprintf(outp, "- Point [%7.2f, %7.2f, %7.2f] is inside solids \"%s\" and \"%s\" (n = %ld)\n", x, y, z, GetText(ptr + STL_SOLID_PTR_FNAME),
                                GetText(loc0 + STL_SOLID_PTR_FNAME), m);
                      */

                      /* Get pointer to STL bodies */

                      if (ptr > VALID_PTR)
                        loc1 = (long)RDB[ptr + STL_SOLID_PTR_BODY];
                      else
                        loc1 = NULLPTR;

                      if (loc0 > VALID_PTR)
                        loc2 = (long)RDB[loc0 + STL_SOLID_PTR_BODY];
                      else
                        loc2 = NULLPTR;

                      /* Print */

                      if ((loc2 < VALID_PTR) && (loc1 > VALID_PTR))
                        fprintf(outp, "- Point [%7.2f, %7.2f, %7.2f] is outside and inside solid \"%s\" (n = %ld)\n", x, y, z, GetText(loc1 + STL_BODY_PTR_BNAME), m);
                      else if ((loc2 > VALID_PTR) && (loc1 < VALID_PTR))
                        fprintf(outp, "- Point [%7.2f, %7.2f, %7.2f] is outside and inside solid \"%s\" (n = %ld)\n", x, y, z, GetText(loc2 + STL_BODY_PTR_BNAME), m);
                      else
                        fprintf(outp, "- Point [%7.2f, %7.2f, %7.2f] is inside solids \"%s\" and \"%s\" (n = %ld)\n", x, y, z, GetText(loc1 + STL_BODY_PTR_BNAME),
                                GetText(loc2 + STL_BODY_PTR_BNAME), m);

                      /* Add to counter */

                      nf++;
                    }

                    /* Break loop */

                    break;
                  }
              }
          }
      }

      if (nf == 0)
        fprintf(outp, "Consistency test passed in all random points.\n\n");
      else
        fprintf(outp, "\nConsistency test passed in %1.5f%% of random points.\n\n",
                100.0*(1.0 - (double)nf/((double)np)));


      /* Pointer to next geometry */

      stl = NextItem(stl);
    }

  /* Stop after testing */

  WDB[DATA_VOLUME_CALCULATION_MODE] = (double)YES;
}

/*****************************************************************************/
