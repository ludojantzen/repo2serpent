/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processnests.c                                 */
/*                                                                           */
/* Created:       2010/10/07 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Creates nets cells, etc.                                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessNests:"

/*****************************************************************************/

void ProcessNests()
{
  long nst, ptr, n, reg, cell, prev;
  char str[MAX_STR];

  /* Loop over nests */

  nst = (long)RDB[DATA_PTR_NST0];
  while (nst > VALID_PTR)
    {
      /* Get pointer to regions */
              
      ptr = (long)RDB[nst + NEST_PTR_REGIONS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Allocate memory for collision region */

      AllocValuePair(nst + NEST_PTR_COL_REG);

      /* Loop over regions */ 
              
      n = 0;
      while ((reg = ListPtr(ptr, n++)) > 0)
        {
          /* Create new cell */

          cell = NewItem(DATA_PTR_C0, CELL_BLOCK_SIZE);

          /* Header data */
          
          WDB[cell + PARAM_PTR_NAME] = RDB[nst + PARAM_PTR_NAME];
          WDB[cell + PARAM_PTR_FNAME] = RDB[nst + PARAM_PTR_FNAME];
          WDB[cell + PARAM_LINE] = RDB[nst + PARAM_LINE];

          /* Put name */

          sprintf(str, "nst%sc%ld", GetText(nst + NEST_PTR_NAME), n);
          WDB[cell + CELL_PTR_NAME] = (double)PutText(str);

          /* Put universe */

          WDB[cell + CELL_PTR_UNI] = RDB[nst + NEST_PTR_UNI];
          
          /* Copy direct material pointer (from divideburnmat.c) */

          WDB[cell + CELL_PTR_REG_MAT] = RDB[reg + NEST_REG_TMP_PTR];
          
          /* Material name and fill pointers (material name is temporarily */
          /* stored in cell pointer) */

          WDB[cell + CELL_PTR_MAT] = RDB[reg + NEST_REG_PTR_CELL];
          WDB[cell + CELL_PTR_FILL] = RDB[reg + NEST_REG_PTR_FILL];

          /* Allocate memory for collision counter */

          AllocValuePair(cell + CELL_COL_COUNT);

          /* Put cell pointer */

          WDB[reg + NEST_REG_PTR_CELL] = (double)cell;

          /* Pointer to out surface */

          if ((prev = PrevItem(reg)) > 0)
            WDB[reg + NEST_REG_PTR_SURF_OUT] =
              RDB[prev + NEST_REG_PTR_SURF_IN];
        }

      /* Next nest */

      nst = NextItem(nst);
    }
}

/*****************************************************************************/
