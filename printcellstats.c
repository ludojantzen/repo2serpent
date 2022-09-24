/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printcellstats.c                               */
/*                                                                           */
/* Created:       2018/06/06 (JLe)                                           */
/* Last modified: 2018/06/06 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Prints cell search statistics etc.                           */
/*                                                                           */
/* Comments: - Used for debugging and optimization.                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintCellStats:"

/*****************************************************************************/

void PrintCellStats()
{
  long uni, loc0, cell, ptr, sum, n;
  return;
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

      /* Print name */

      fprintf(outp, "\nUniverse %s:\n\n", GetText(uni + UNIVERSE_PTR_NAME));

      /* Reset sum */

      sum = 0.0;

      /* Loop over list */

      while (loc0 > VALID_PTR)
        {
          /* Pointer to cell */

          cell = (long)RDB[loc0 + CELL_LIST_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Print name */

          fprintf(outp, "%10s ", GetText(cell + CELL_PTR_NAME));

          /* Pointer to surface list */

          ptr = (long)RDB[cell + CELL_PTR_SURF_LIST];

          /* Count number of surfaces */

          n = 0;
          while ((long)RDB[ptr + n] > VALID_PTR)
            n++;

          /* Print number of surfaces */

          fprintf(outp, "%4ld ", n);

          /* Pointer to cell search list count */

          ptr = (long)RDB[loc0 + CELL_LIST_PTR_COUNT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Add to sum */

          sum = sum + (long)SumPrivateData(ptr);

          /* Print */

          fprintf(outp, "%ld\n", (long)SumPrivateData(ptr));

          /* Next cell in list */

          loc0 = NextItem(loc0);
        }

      /* Print total */

      fprintf(outp, "Total:          %ld\n", sum);
            
      /* Next universe */

      uni = NextItem(uni);
    }          
}

/*****************************************************************************/
