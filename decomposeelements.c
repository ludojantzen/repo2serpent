/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : decomposeelements.c                            */
/*                                                                           */
/* Created:       2017/12/12 (JLe)                                           */
/* Last modified: 2017/12/12 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Decomposes elemental data into isotopic in material          */
/*              compositions.                                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "element_data.h"

#define FUNCTION_NAME "DecomposeElements:"

/*****************************************************************************/

void DecomposeElements()
{
  long mat, iso, ptr, n;
  char name[MAX_STR], id[MAX_STR], dens[MAX_STR];

  /* Check decomposition */

  if (((long)RDB[DATA_ELEM_DECOMP] == NO) && 
      ((long)RDB[DATA_ELEM_DECOMP_PTR_LIST] < VALID_PTR))
    return;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Loop over composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Get name */

          sprintf(name, "%s", GetText(iso + COMPOSITION_PTR_NUCLIDE));
          
          /* Check length */

          if (strlen(name) < 7)
            {
              /* Next nuclide */

              iso = NextItem(iso);

              /* Cycle loop */

              continue;
            }
          
          /* Separate name and id */

          sprintf(id, "%s", &name[strlen(name) - 3]);
          name[strlen(name) - 4] = '\0';

          /* Check for neutron data */

          if (id[2] != 'c')
            {
              /* Next nuclide */

              iso = NextItem(iso);

              /* Cycle loop */

              continue;
            }

          /* Check for natural element */
          
          if (!strcmp(&name[strlen(name) - 3], "000"))
            name[strlen(name) - 3] = '\0';
          else if (!strcmp(&name[strlen(name) - 3], "nat"))
            name[strlen(name) - 4] = '\0';
          else
            {
              /* Next nuclide */

              iso = NextItem(iso);

              /* Cycle loop */

              continue;
            }

          /* Check if list is given */

          if ((ptr = (long)RDB[DATA_ELEM_DECOMP_PTR_LIST]) > VALID_PTR)
            {
              /* Loop over list and find match */

              while ((long)RDB[ptr] > VALID_PTR)
                {
                  /* Compare */

                  if (!strcmp(GetText(ptr), name))
                    break;

                  /* Loop over element symbols to find match */

                  for (n = 0; n < NUMBER_OF_ELEMENTS; n++)
                    if (!strcasecmp(element_symbols[n], GetText(ptr)))
                      break;
                  
                  /* Check if found */
                  
                  if (n == atol(name))
                    break;                        

                  /* Loop over element names to find match */

                  for (n = 0; n < NUMBER_OF_ELEMENTS; n++)
                    if (!strcasecmp(element_names[n], GetText(ptr)))
                      break;
                      
                  /* Check if found */
                  
                  if (n == atol(name))
                    break;                        

                  /* Update pointer */

                  ptr++;
                }
              
              /* Check if found */
              
              if ((((long)RDB[DATA_ELEM_DECOMP] == YES) &&
                   ((long)RDB[ptr] < VALID_PTR)) ||
                  (((long)RDB[DATA_ELEM_DECOMP] == NO) &&
                   ((long)RDB[ptr] > VALID_PTR)))
                {
                  /* Next nuclide */

                  iso = NextItem(iso);
                  
                  /* Cycle loop */
                  
                  continue;
                }
            }

          /* Get density */

          sprintf(dens, "%1.10E", RDB[iso + COMPOSITION_ADENS]);
          
          /* Decompose */

          fprintf(outp, "Decomposing %s in into isotopes...\n",
                  GetText(iso + COMPOSITION_PTR_NUCLIDE));

          Element(iso, name, dens, id);

          /* Next */

          iso = NextItem(iso);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  fprintf(outp, "\n");
}

/*****************************************************************************/
