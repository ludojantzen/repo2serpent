/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : universeboundaries.c                           */
/*                                                                           */
/* Created:       2010/11/02 (JLe)                                           */
/* Last modified: 2018/06/17 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Finds the outer boundaries for each universe                 */
/*                                                                           */
/* Comments: -Transformaatiot ei vaikuta mittoihin (processtransformations.c */
/*            printtaa varoituksen jos niitä käytetään alimmalle             */
/*            universumille)                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UniverseBoundaries:"

/*****************************************************************************/

void UniverseBoundaries()
{
  long uni, reg, cell, surf, ptr, loc0, nst, n;
  double minx, maxx, miny, maxy, minz, maxz;

  /* Loop over universes */

  uni = (long)RDB[DATA_PTR_U0];
  while (uni > 0)
    {
      /* Reset boundaries */

      minx =  INFTY;
      maxx = -INFTY;
      miny =  INFTY;
      maxy = -INFTY;
      minz =  INFTY;
      maxz = -INFTY;

      /* Check type */

      switch ((long)RDB[uni + UNIVERSE_TYPE])
        {
        case UNIVERSE_TYPE_NEST:
          {
            /***** Nest universe *********************************************/

            /* Pointer to nest */

            nst = (long)RDB[uni + UNIVERSE_PTR_NEST];
            CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);
            
            /* Pointer to outermost region */

            reg = (long)RDB[nst + NEST_PTR_REGIONS];
            CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

            reg = LastItem(reg);
            CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

            /* Pointer to surface (if any) */

            if ((surf = (long)RDB[reg + NEST_REG_PTR_SURF_OUT]) > 0)
              {
                /* Check pointer */
                
                CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);
                
                /* Get dimensions */

                MaxSurfDimensions(surf, &minx, &maxx, &miny, &maxy, 
                                  &minz, &maxz);
              }

            /* Break case */

            break;
            
            /*****************************************************************/
          }      

        case UNIVERSE_TYPE_CELL:
          {
            /***** Cell universe *********************************************/
        
            /* Pointer to cell list */
        
            ptr = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  

            /* Loop over cells */
  
            n = 0;

            while ((cell = ListPtr(ptr, n++)) > 0)
              {
                /* Pointer to cell */

                cell = (long)RDB[cell + CELL_LIST_PTR_CELL];
                CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

                /* Pointer to surface list */

                loc0 = (long)RDB[cell + CELL_PTR_SURF_LIST];
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                /* Loop over list */

                while ((surf = (long)RDB[loc0++]) > 0)
                  {
                    /* Check for infinite type */

                    if (((long)RDB[surf + SURFACE_TYPE] == SURF_INF) &&
                        (uni == (long)RDB[DATA_PTR_ROOT_UNIVERSE]))
                      Error(surf, "Infinite surface not allowed in universe %s",
                            GetText(uni + UNIVERSE_PTR_NAME));

                    /* Get dimensions */

                    MaxSurfDimensions(surf, &minx, &maxx, &miny, &maxy, 
                                      &minz, &maxz);
                  }
              }

            /* Break case */
            
            break;
            
            /*****************************************************************/
          }      
        }

      /* Reset boundaries if infintie */

      if ((minx == INFTY) && (maxx == -INFTY))
        {
          minx = 0.0;
          maxx = 0.0;
        }

      if ((miny == INFTY) && (maxy == -INFTY))
        {
          miny = 0.0;
          maxy = 0.0;
        }

      if ((minz == INFTY) && (maxz == -INFTY))
        {
          minz = 0.0;
          maxz = 0.0;
        }

      /* Put values */

      WDB[uni + UNIVERSE_MINX] = minx;
      WDB[uni + UNIVERSE_MAXX] = maxx;
      WDB[uni + UNIVERSE_MINY] = miny;
      WDB[uni + UNIVERSE_MAXY] = maxy;
      WDB[uni + UNIVERSE_MINZ] = minz;
      WDB[uni + UNIVERSE_MAXZ] = maxz;

      /* Put dimension */

      if ((minz == 0.0) && (maxz == 0.0))
        WDB[uni + UNIVERSE_DIM] = 2.0;
      else
        WDB[uni + UNIVERSE_DIM] = 3.0;
        
      /* Next universe */

      uni = NextItem(uni);
    }

  /* Get pointer to root universe */

  uni =  (long)RDB[DATA_PTR_ROOT_UNIVERSE];
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);
  
  /* Put geometry boundaries */

  WDB[DATA_GEOM_MINX] = RDB[uni + UNIVERSE_MINX];
  WDB[DATA_GEOM_MAXX] = RDB[uni + UNIVERSE_MAXX];
  WDB[DATA_GEOM_MINY] = RDB[uni + UNIVERSE_MINY];
  WDB[DATA_GEOM_MAXY] = RDB[uni + UNIVERSE_MAXY];
  WDB[DATA_GEOM_MINZ] = RDB[uni + UNIVERSE_MINZ];
  WDB[DATA_GEOM_MAXZ] = RDB[uni + UNIVERSE_MAXZ];

  /* Check */

  if ((RDB[DATA_GEOM_MINX] > RDB[DATA_GEOM_MAXX]) ||
      (RDB[DATA_GEOM_MINY] > RDB[DATA_GEOM_MAXY]) ||
      ((RDB[uni + UNIVERSE_DIM] == 3.0) &&
       (RDB[DATA_GEOM_MINZ] > RDB[DATA_GEOM_MAXZ])))
    Error(0, "Error in geometry boundaries (check definition of universe %s)",
          GetText(uni + UNIVERSE_PTR_NAME));

  /* Put geometry dimension */
  
  WDB[DATA_GEOM_DIM] = RDB[uni + UNIVERSE_DIM];
}

/*****************************************************************************/
