/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkunused.c                                  */
/*                                                                           */
/* Created:       2013/02/19 (JLe)                                           */
/* Last modified: 2019/11/05 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Checks and removes unused stuff                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckUnused:"

/*****************************************************************************/

void CheckUnused()
{
  long ptr;

  /* Cells */

  ptr = (long)RDB[DATA_PTR_C0];
  while (ptr > VALID_PTR)
    {
      /* Check */

      if (!((long)RDB[ptr + CELL_OPTIONS] & OPT_USED))
        Note(ptr, "Cell %s is not used in geometry",
             GetText(ptr + CELL_PTR_NAME));

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Remove */

  ptr = (long)RDB[DATA_PTR_C0];
  RemoveFlaggedItems(ptr, CELL_OPTIONS, OPT_USED, NO);

  /* Nests */

  ptr = (long)RDB[DATA_PTR_NST0];
  while (ptr > VALID_PTR)
    {
      /* Check */

      if (!((long)RDB[ptr + NEST_OPTIONS] & OPT_USED))
        {
          /* Check if nest contains burnable material */

          if (((long)RDB[ptr + NEST_OPTIONS] & OPT_BURN_MAT) &&
              ((long)RDB[DATA_PTR_MVOL0] > VALID_PTR))
            Note(ptr,
                 "Nest %s is not used in geometry, depletion zone\n                indexing used in the \"set mvol\" card may be off", GetText(ptr + NEST_PTR_NAME));
          else
            Note(ptr, "Nest %s is not used in geometry",
                 GetText(ptr + NEST_PTR_NAME));
        }

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Remove */

  ptr = (long)RDB[DATA_PTR_NST0];
  RemoveFlaggedItems(ptr, NEST_OPTIONS, OPT_USED, NO);

  /* Remove unused lattices */

  ptr = (long)RDB[DATA_PTR_L0];
  while (ptr > VALID_PTR)
    {
      /* Check */

      if (!((long)RDB[ptr + LAT_OPTIONS] & OPT_USED))
        Note(ptr, "Lattice %s is not used in geometry",
             GetText(ptr + LAT_PTR_NAME));

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Remove */

  ptr = (long)RDB[DATA_PTR_L0];
  RemoveFlaggedItems(ptr, LAT_OPTIONS, OPT_USED, NO);

  /* Remove unused transformations */

  ptr = (long)RDB[DATA_PTR_TR0];
  while (ptr > VALID_PTR)
    {
      /* Check */

      if (!((long)RDB[ptr + TRANS_OPTIONS] & OPT_USED))
        Note(ptr, "Transformation %s is not used in geometry",
             GetText(ptr + TRANS_PTR_NAME));

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Remove */

  ptr = (long)RDB[DATA_PTR_TR0];
  RemoveFlaggedItems(ptr, TRANS_OPTIONS, OPT_USED, NO);
}

/*****************************************************************************/
