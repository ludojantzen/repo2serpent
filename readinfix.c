/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readinfix.c                                    */
/*                                                                           */
/* Created:       2010/10/18 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Converts surface list from text string to infix notation.    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadInfix:"

/*****************************************************************************/

long ReadInfix(long cell, long *infix, long *nspec)
{
  long i0, i, nparam;
  char str[10000], param[MAX_STR];
  
  /* Check cell pointer */

  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

  /* Check surface list pointer */

  if ((long)RDB[cell + CELL_PTR_SURF_LIST] < VALID_PTR)
    Error(cell, "Cell %s has no surface definitions", 
          GetText(cell + CELL_PTR_NAME));
  
  /* Copy string to a variable */

  sprintf(str, "%s", GetText(cell + CELL_PTR_SURF_LIST));

  /* Reset number of special operations */

  *nspec = 0;

  /* Read parameters in infix form and add intersections */

  i0 = 0;
  nparam = 0;

  while ((i = NextWord(&str[i0], param)) > 1)
    {
      if (param[0] == ':')
        {
          /* Union operator */
          
          infix[nparam++] = SURF_OP_OR;

          /* Add special operation count */
          
          (*nspec)++;
        }
      else if (param[0] == '-')
        {
          /* Complement operator, check preceding term */
          
          if ((nparam > 0) && ((infix[nparam - 1] == SURF_OP_RIGHT) ||
                            (infix[nparam - 1] > 0)))
            {
              /* Add intersection */
              
              infix[nparam++] = SURF_OP_AND;
            }
          
          /* Add complement */
          
          infix[nparam++] = SURF_OP_NOT;
        }
      else if (param[0] == '(')
        {
          /* Left parenthesis, check preceding term */
          
          if ((nparam > 0) && ((infix[nparam - 1] == SURF_OP_RIGHT) ||
                            (infix[nparam - 1] > 0)))
            {
              /* Add intersection */
              
              infix[nparam++] = SURF_OP_AND;
            }
          
          /* Add left parenthesis */
          
          infix[nparam++] = SURF_OP_LEFT;
          
          /* Add special operation count */
          
          (*nspec)++;
        }
      else if (param[0] == ')')
        {
          /* Right parenthesis */
          
          infix[nparam++] = SURF_OP_RIGHT;
          
          /* Add special operation count */
          
          (*nspec)++;
        }
      else if ((nparam > 0) && ((infix[nparam - 1] == SURF_OP_RIGHT) ||
                             (infix[nparam - 1] > 0)))
        {
          /* Add intesection */
          
          infix[nparam++] = SURF_OP_AND;
          
          /* Add surface name */
          
          infix[nparam++] = PutText(param);
        }
      else
        {
          /* Add surface name */
          
          infix[nparam++] = PutText(param);
        }
      
      /* Next word */
      
      i0 = i0 + i;
    }

  /* Return number of parameters */

  return nparam;
}

/*****************************************************************************/
