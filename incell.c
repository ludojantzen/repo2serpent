/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : incell.c                                       */
/*                                                                           */
/* Created:       2010/10/12 (JLe)                                           */
/* Last modified: 2014/02/09 (JLe)                                           */
/* Version:       2.1.17                                                     */
/*                                                                           */
/* Description: Checks if point is inside cell                               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InCell:"

/*****************************************************************************/

long InCell(long cell, double x, double y, double z, long on, long id)
{
  long n, ptr, loc0, surf, side, a, b, stack[10000];

  /* Check cell pointer */

  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

  /* Check type */

  if ((ptr = (long)RDB[cell + CELL_PTR_SURF_INSC]) > VALID_PTR)
    {
      /***********************************************************************/

      /***** Intersection list given *****************************************/

      /* Loop over list */

      n = 0;
      while ((loc0 = ListPtr(ptr, n++)) > VALID_PTR)
        {
          /* Pointer to surface */

          surf = (long)RDB[loc0 + CELL_INSC_PTR_SURF];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Surface side */

          side = (long)RDB[loc0 + CELL_INSC_SIDE];

          /* Test surface */

          if (TestSurface(surf, x, y, z, on, id) == YES)
            side = -side;

          /* Check result */

          if (side < 0)
            {
              /* Add count */

              ptr = (long)RDB[loc0 + CELL_INSC_PTR_OUT_COUNT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddPrivateData(ptr, 1.0, id);
              
              /* Point is out */

              return NO;
            }
        }
      
      /* Point is in */

      return YES;

      /***********************************************************************/
    }
  else if ((ptr = (long)RDB[cell + CELL_PTR_SURF_COMP]) > VALID_PTR)
    {
      /***********************************************************************/
      
      /***** Composition list given ******************************************/

      /* Reset counter */

      n = 0;

      /* Loop over composition */
      
      while ((long)RDB[ptr] != 0)
        {
          /* Check parameter type */
          
          if ((long)RDB[ptr] == SURF_OP_OR)
            {
              /* Check counter */
              
              CheckValue(FUNCTION_NAME, "n", "(1)", n, 2, 10000);
              
              /* Pop last two values */
              
              a = stack[--n];
              b = stack[--n];
             
              /* Push union to stack */
              
              stack[n++] = a | b;
            }
          else if ((long)RDB[ptr] == SURF_OP_NOT)
            {
              /* Check counter */
              
              CheckValue(FUNCTION_NAME, "n", "(2)", n, 1, 10000);
              
              /* Pop last value */
              
              a = stack[--n];
              
              /* Take complement */
              
              a ^= 1;
              
              /* Push value to stack */
              
              stack[n++] = a;
            }
          else if ((long)RDB[ptr] == SURF_OP_AND)
            {
              /* Check counter */
              
              CheckValue(FUNCTION_NAME, "n", "(3)", n, 2, 10000);
              
              /* Pop last two values */
              
              a = stack[--n];
              b = stack[--n];
              
              /* Push intersection to stack */
              
              stack[n++] = a & b;
            }
          else
            {
              /* Check counter */
              
              CheckValue(FUNCTION_NAME, "n", "(4)", n, 0, 10000);
              
              /* Test and push result to stack */
              
              if (TestSurface((long)RDB[ptr], x, y, z, on, id) == YES)
                stack[n++] = 0;
              else
                stack[n++] = 1;
            }
          
          /* Next */
          
          ptr++;
        }
      
      /* Check stack status */

      if (n != 1)
        Die(FUNCTION_NAME, "Horrible error in cell %s (n = %ld)", 
            GetText(cell + CELL_PTR_NAME), n);
      
      /* Return final result */
      
      if (stack[0] == 1)
        return YES;
      else
        return NO;
  
      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "No lists");

  /* Avoid compiler warning */

  return 0;
}

/*****************************************************************************/
