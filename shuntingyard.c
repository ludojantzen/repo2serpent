/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : shuntingyard.c                                 */
/*                                                                           */
/* Created:       2010/10/18 (JLe)                                           */
/* Last modified: 2019/03/26 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Converts infix notation to postfix using the shunting yard   */
/*              algorithm.                                                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ShuntingYard:"

/*****************************************************************************/

void ShuntingYard(long cell, long *infix, long ni)
{
  long np,  ptr, n;
  static long stack[10000];

  /* Check cell pointer */

  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

  /* Allocate memory for composition list */
              
  ptr = ReallocMem(DATA_ARRAY, ni + 1);
  WDB[cell + CELL_PTR_SURF_COMP] = (double)ptr;
              
  /* Reset count */
  
  np = 0;

  /* Loop over values */
  
  for (n = 0; n < ni; n++)
    {
      /* Check error */
      
      if ((n > 0) && (infix[n] == SURF_OP_OR) &&
          (infix[n - 1] < 0) && (infix[n - 1] != SURF_OP_RIGHT))
        Error(cell, "Error in surface list");

      /* Check overflow */

      if (np > 9999)
        Die(FUNCTION_NAME, "Array overflow");
      
      /* Handle parameters */
      
      if ((infix[n] < 0) && (infix[n] != SURF_OP_RIGHT) && 
          (infix[n] != SURF_OP_LEFT))
        {
          /* Pop operators with higher precedence */
          
          while ((np > 0) && (stack[np - 1] > infix[n]) )
            WDB[ptr++] = (double)stack[--np];

          /* Push operator */
          
          stack[np++] = infix[n];
        }
      else if (infix[n] == SURF_OP_LEFT)
        {
          /* Push left parenthesis */
          
          stack[np++] = infix[n];
        }
      else if (infix[n] == SURF_OP_RIGHT)
        {
          /* Check parenthesis error */
          
          if (--np < 0)
            Error(cell, "Missing left parenthesis");
          
          /* Pull stack */
          
          while (stack[np] != SURF_OP_LEFT)
            {
              WDB[ptr++] = (double)stack[np];
              
              /* Check error */

              if (--np < 0)
                Error(cell, "Missing left parenthesis");
            }
        }
      else
        WDB[ptr++] = (double)infix[n];
    }
  
  /* Check parenthesis error */
  
  for (n = 0; n < np; n++)
    if (stack[n] == SURF_OP_LEFT)
      Error(cell, "Missing right parenthesis");
  
  /* Pull stack */
  
  while (--np > -1)
    WDB[ptr++] = (double)stack[np];
  
  /* Put null terminator */
  
  WDB[ptr] = 0.0;
}


/*****************************************************************************/
