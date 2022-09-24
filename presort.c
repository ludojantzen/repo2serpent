/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : presort.c                                      */
/*                                                                           */
/* Created:       2018/06/08 (JLe)                                           */
/* Last modified: 2018/06/08 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Samples points in geometry to pre-sort cell lists            */
/*                                                                           */
/* Comments: - Tää ei oikein toimi käytännössä, sillä useimmat search        */
/*             listit on cell-kohtaisia, eikä toi geometriahaku päivitä      */
/*             tässä niiden laskureita.                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PreSort:"

/*****************************************************************************/

void PreSort()
{
  long cell, uni, loc0, loc1, ptr, nb, np, idx, n, m, id;
  unsigned long seed;
  double xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, u, v, w;

  /* Get number of points */

  if ((np = (long)RDB[DATA_PRESORT_NP]) < 1)
    return;

  fprintf(outp, "Pre-sorting cell search lists...\n");

  /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */

  ExpandPrivateArrays();

  /* Get geometry boundaries */

  xmin = RDB[DATA_GEOM_MINX];
  xmax = RDB[DATA_GEOM_MAXX];
  ymin = RDB[DATA_GEOM_MINY];
  ymax = RDB[DATA_GEOM_MAXY];
  zmin = RDB[DATA_GEOM_MINZ];
  zmax = RDB[DATA_GEOM_MAXZ];

  /* Get number of batches */

  nb = (long)RDB[DATA_PRESORT_NB];

  /***************************************************************************/

  /***** Sample random points in geometry ************************************/

  /* Loop over batches */

  for (n = 0; n < nb; n++)
    {
#ifdef OPEN_MP
#pragma omp parallel private (m, idx, seed, uni, ptr, x, y, z, u, v, w, cell, id)
#endif
      {
        /* Loop over points */
        
#ifdef OPEN_MP
#pragma omp for      
#endif
        for (m = 0; m < np; m++)
          {
            /* Get OpenMP thread num */
            
            id = OMP_THREAD_NUM;
            
            /* Calculate index */
            
            idx = nb*np + m + 1;
            
            /* Init random number sequence */
            
            seed = ReInitRNG(idx);
            SEED[id*RNG_SZ] = seed;
            
            /* Sample point */
            
            x = RandF(id)*(xmax - xmin) + xmin;
            y = RandF(id)*(ymax - ymin) + ymin;
            
            if ((long)RDB[DATA_GEOM_DIM] == 3)
              z = RandF(id)*(zmax - zmin) + zmin;
            else
              z = 0.0;
            
            /* Sample direction (this is necessary for STL geometries) */
                
            IsotropicDirection(&u, &v, &w, id);

            /* Loop over universes */

            uni = (long)RDB[DATA_PTR_U0];
            while (uni > VALID_PTR)
              {
                /* Reset previous region */

                if ((ptr = (long)RDB[uni + UNIVERSE_PTR_PREV_REG]) > VALID_PTR)
                  PutPrivateData(ptr, -1.0, id);

                /* Next universe */

                uni = NextItem(uni);
              }
            
            /* Set cell search list option */

            ptr = (long)RDB[DATA_CELL_SEARCH_LIST];
            CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
            PutPrivateData(ptr, (double)CELL_SEARCH_LIST_UNI, id);

            /* Find position */
                
            if ((cell = WhereAmI(x, y, z, u, v, w, id)) < 0)
              Error(0, "Geometry error at %E %E %E", x, y, z);
          }
      }

      /* Sort lists */
      
      SortAll();
    }

  /***************************************************************************/

  /***** Reset counters ******************************************************/

  /* Loop over universes */

  uni = (long)RDB[DATA_PTR_U0];
  while (uni > VALID_PTR)
    {
      /* Pointer to cell list */
        
      if ((loc0 = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST]) < VALID_PTR)
        {
          /* Not a universe cell */

          uni = NextItem(uni);

          /* Cycle loop */

          continue;
        }

      /* Loop over list */

      while (loc0 > VALID_PTR)
        {
          /* Pointer to cell search list count */

          ptr = (long)RDB[loc0 + CELL_LIST_PTR_COUNT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          /* Reset */

          for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
            PutPrivateData(ptr, 0.0, id);

          /* Pointer to cell */

          cell = (long)RDB[loc0 + CELL_LIST_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Loop over search list */

          loc1 = (long)RDB[cell + CELL_PTR_SEARCH_LIST];
          while (loc1 > VALID_PTR)
            {
              /* Pointer to cell search list count */

              ptr = (long)RDB[loc1 + CELL_LIST_PTR_COUNT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              
              /* Reset */

              for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
                PutPrivateData(ptr, 0.0, id);

              /* Next */

              loc1 = NextItem(loc1);
            }          

          /* Next cell in list */

          loc0 = NextItem(loc0);
        }

      /* Next universe */

      uni = NextItem(uni);
    }          

  /***************************************************************************/

  /* Exit subroutine */

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
