/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : removevoidcells.c                              */
/*                                                                           */
/* Created:       2013/10/15 (JLe)                                           */
/* Last modified: 2016/09/08 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Removes void cells                                           */
/*                                                                           */
/* Comments: - Used as a special trick to accelerate traking in complex      */
/*             geometries with a lot of void.                                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RemoveVoidCells:"

/*****************************************************************************/

void RemoveVoidCells()
{
  long cell, ptr;

  /* Check option */

  if ((long)RDB[DATA_REMOVE_VOID_CELLS] == NO)
    return;

  /* Loop over cells */

  cell = (long)RDB[DATA_PTR_C0];
  while (cell > VALID_PTR)
    {
      /* Copy pointer */

      ptr = cell;

      /* Next item */

      cell = NextItem(cell);

      /* Check for void and remove */
      
      if ((long)RDB[ptr + CELL_PTR_MAT] > VALID_PTR)
        if (!strcmp(GetText(ptr + CELL_PTR_MAT), "void"))
          RemoveItem(ptr);
    }

  /* Set mode and maximum number of lost particles */

  WDB[DATA_IGNORE_UNDEFINED_CELLS] = (double)YES;
  WDB[DATA_MAX_UNDEF_POS] = -1.0;
}

/*****************************************************************************/
