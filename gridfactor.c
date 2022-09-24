/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : gridfactor.c                                   */
/*                                                                           */
/* Created:       2010/12/15 (JLe)                                           */
/* Last modified: 2013/10/15 (JLe)                                           */
/* Version:       2.1.16                                                      */
/*                                                                           */
/* Description: Calculates interpolation factor for energy grid              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GridFactor:"

/*****************************************************************************/

double GridFactor(long erg, double E, long id)
{
  long i, ptr;
  double f, E0, E1;

  /* Check Pointer */

  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Test existing data */

  if ((f = TestValuePair(erg + ENERGY_GRID_PTR_PREV_VAL, E, id)) > -INFTY)
    return f;

  /* Find grid index */

  if ((i = GridSearch(erg, E)) < 0)
    f = -1.0;
  else
    {
      /* Pointer to data */

      ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Boundary values */
      
      E0 = RDB[ptr + i];
      E1 = RDB[ptr + i + 1];
      
      /* Interpolate (nimittäjän nollakohdan tarkistus 15.10.2013 / 2.1.16) */
      
      if (E0 != E1)
        f = (E - E0)/(E1 - E0);
      else
        f = 0.0;

      /* Check value */
      
      CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);
      
      /* Add index */
      
      f = f + (double)i;
    }

  /* Remember value */

  StoreValuePair(erg + ENERGY_GRID_PTR_PREV_VAL, E, f, id);

  /* Return factor */

  return f;
}

/*****************************************************************************/
