/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : insupercell.c                                  */
/*                                                                           */
/* Created:       2011/03/02 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Tests if point is in a super-imposed cell or universe        */
/*                                                                           */
/* Comments: - Karsittu versio whereami.c:stä                                */
/*                                                                           */
/*           - TODO: lisää fysikaalisten cellien testaus ja muuta nimi       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InSuperCell:"

/*****************************************************************************/

long InSuperCell(long uni, long c0, double x, double y, double z, long id)
{
  long cell, ptr, dummy;
  double u, v, w;

  /* Dummy direction cosines */

  u = 0.0;
  v = 0.0; 
  w = 1.0;

  /* Search loop */

  while (1 != 2)
    {
      /* Check universe pointer */
  
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Do coordinate transformation */
      
      if ((ptr = (long)RDB[uni + UNIVERSE_PTR_TRANS]) > VALID_PTR)
        CoordTrans(ptr, &x, &y, &z, &u, &v, &w, id);

      /* Universe symmetry */

      if ((ptr = (long)RDB[uni + UNIVERSE_PTR_SYM]) > VALID_PTR)
        UniSym(ptr, &x, &y, &z, &u, &v, &w);

      /* Check type */

      if ((long)RDB[uni + UNIVERSE_TYPE] != UNIVERSE_TYPE_SUPER)
        Die(FUNCTION_NAME, "Invalid type");
      
      /* Find cell */

      if ((cell = FindUniverseCell(uni, x, y, z, &dummy, id)) < VALID_PTR)
        return NO;

      /* Check fill pointer and call recursively */
            
      if ((ptr = (long)RDB[cell + CELL_PTR_FILL]) > VALID_PTR)
        {
          /* Filled cell, update universe pointer */
                
          uni = ptr;
        }
      else
        {
          /* Compare to cell pointer if given */

          if ((c0 > VALID_PTR) && (c0 == cell))
            return YES;

          /* All universe, check cell type */

          if ((c0 < VALID_PTR) && 
              ((long)RDB[cell + CELL_TYPE] != CELL_TYPE_OUTSIDE))
            return YES;

          /* Break loop */

          break;
        }
    }
  
  /* Not in */

  return NO;

  /***************************************************************************/
}

/*****************************************************************************/
