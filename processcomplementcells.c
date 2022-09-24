/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processcomplementcells.c                       */
/*                                                                           */
/* Created:       2015/07/17 (JLe)                                           */
/* Last modified: 2017/05/21 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Processes MCNP style cell complements in surface lists       */
/*                                                                           */
/* Comments: - Tota new-arrayn kokoa ei testata, antaa segfaultin jos        */
/*             se ylittyy.                                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessComplementCells:"

/*****************************************************************************/

void ProcessComplementCells()
{
  long cell, ptr, i, i0, n, con, nested;
  char *lst, new[20000], word[MAX_STR];

  /* Outer loop */

  do
    {
      /* Reset nested */

      nested = NO;

      /* Loop over cells */
      
      cell = (long)RDB[DATA_PTR_C0];
      while (cell > VALID_PTR)
        {
          /* Pointer to surface list */
          
          if ((long)RDB[cell + CELL_PTR_SURF_LIST] > VALID_PTR)
            {
              /* Get string */
              
              lst = GetText(cell + CELL_PTR_SURF_LIST);
              
              /* Reset new and conversion flags */
              
              new[0] = '\0';
              con = NO;
              
              /* Loop over string */
              
              i0 = 0;
              while ((i = NextWord(&lst[i0], word)) > 0)
                {
                  /* Update index */
                  
                  i0 = i0 + i;

                  /* Check hashtag */
                  
                  if (word[0] != '#')
                    sprintf(&new[strlen(new)], " %s", word);
                  else
                    {
                      /* Loop over cells in parenthesis */
                      
                      n = 0;
                      
                      do
                        {
                          /* Get next word in list */
                          
                          if ((i = NextWord(&lst[i0], word)) < 1)
                            Error(cell, "Error in surface list");
                          else
                            {
                              /* Update index */
                              
                              i0 = i0 + i;
                              
                              /* Check */
                              
                              if (word[0] == '(')
                                {
                                  if (n == 0)
                                    sprintf(&new[strlen(new)], " - (");
                                  else
                                    sprintf(&new[strlen(new)], " (");
                                  n++;
                                }
                              else if (word[0] == ')')
                                {
                                  sprintf(&new[strlen(new)], " )");
                                  n--;
                                }
                              else if (n > 0)
                                sprintf(&new[strlen(new)], " %s", word);
                              else
                                {
                                  /* Find cell */
                                  
                                  ptr = (long)RDB[DATA_PTR_C0];
                                  while (ptr > VALID_PTR)
                                    {
                                      /* Compare names */
                                      
                                      if (!strcmp(GetText(ptr + CELL_PTR_NAME), 
                                                  word))
                                        break;
                                      
                                      /* Next cell */
                                      
                                      ptr = NextItem(ptr);
                                    }
                                  
                                  /* Check */
                                  
                                  if (ptr < VALID_PTR)
                                    Error(cell, 
                                          "Cell %s in complement not defined",
                                          word);
                                  
                                  /* Print cell list in new */
                                  
                                  sprintf(&new[strlen(new)], " - ( %s )", 
                                          GetText(ptr + CELL_PTR_SURF_LIST));
                                }
                            }
                          
                          /* Set conversion flag */
                          
                          con = YES;
                        }
                      while (n > 0);
                    }
                }
              
              /* Check if conversion was made */

              if (con == YES)
                {
                  /* Check that no hastags remain */
                  
                  i = 0;
                  while (new[i] != '\0')
                    if (new[i++] == '#')
                      nested = YES;
                  
                  /* Put new cell list */
                  
                  WDB[cell + CELL_PTR_SURF_LIST] = (double)PutText(new);
                }
            }
          
          /* Next cell */
          
          cell = NextItem(cell);
        }
    }
  while (nested == YES);
}

/*****************************************************************************/
