/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : Addchains.c                                    */
/*                                                                           */
/* Created:       2010/12/28 (JLe)                                           */
/* Last modified: 2013/02/12 (JLe)                                           */
/* Version:       2.1.13                                                     */
/*                                                                           */
/* Description: Adds nuclides in decay and transmutation chains to material  */
/*              compositions.                                                */
/*                                                                           */
/* Comments:  TODO: tähän pitää lisätä se transmutaatiolistan generointi     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddChains:"

/*****************************************************************************/

void AddChains(long mat, long nuc, long lvl)
{
  long rea, new, yld, ptr;

  /* Check nuclide used-flag */

  if ((long)RDB[nuc + NUCLIDE_OPTIONS] & OPT_USED)
    return;
  else
    SetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

  /* Check for lost nuclide */

  if (nuc == (long)RDB[DATA_PTR_NUCLIDE_LOST])
    return;

  /* Check exists-flag */
  
  if (!((long)RDB[nuc + NUCLIDE_OPTIONS] & OPT_EXISTS))
    {
      /* Add nuclide to composition */

      ptr = NewItem(mat + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);

      /* Put pointer */
      
      WDB[ptr + COMPOSITION_PTR_NUCLIDE] = (double)nuc;

      /* Set exists-flag */

      SetOption(nuc + NUCLIDE_OPTIONS, OPT_EXISTS);
      
      /* Reset density */
      
      WDB[ptr + COMPOSITION_ADENS] = 0.0;
    }

  /* Loop over reactions (NOTE: tätä voisi nopeuttaa lukemalla reaktiot */
  /* etukäteen omana listaansa) */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];

  while (rea > VALID_PTR)
    {
      /* Get pointer to target */
      
      if ((new = (long)RDB[rea + REACTION_PTR_TGT]) > VALID_PTR)
        {
          /* Call recursively */
          
          AddChains(mat, new, lvl + 1);
        }
      else if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
        {
          /* Get pointer to distribution */
          
          yld = (long)RDB[yld + FISSION_YIELD_PTR_DISTR];
          CheckPointer(FUNCTION_NAME, "yld", DATA_ARRAY, yld);
          
          /* Loop over distribution */ 

          while (yld > VALID_PTR)
            {
              /* Get pointer to target */
              
              if ((new = (long)RDB[yld + FY_PTR_TGT]) > VALID_PTR)
                {
                  /* Call recursively */
                  
                  AddChains(mat, new, lvl + 1);
                }
              
              /* Next */

              yld = NextItem(yld);
            }
        }
          
      /* Next reaction */
      
      rea = NextItem(rea);
    }  
}

/*****************************************************************************/
