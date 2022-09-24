/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : replacephotondata.c                            */
/*                                                                           */
/* Created:       2015/05/15 (JLe)                                           */
/* Last modified: 2018/03/13 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Replaces isotopic compositions with atomic compositions in   */
/*              photon transport simulation.                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReplacePhotonData:"

/*****************************************************************************/

void ReplacePhotonData()
{
  long mat, iso, nuc, n, nucl[110];
  double dens[110];

  /* Check mode */

  if (((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES) ||
      ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == NO))
    return;

  /* Check domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    Die(FUNCTION_NAME, "Domain decomposition in use");

  /* Loop over materials */
 
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Reset density and nuclide arrays */

      for (n = 0; n < 110; n++)
        {
          nucl[n] = 0;
          dens[n] = 0.0;
        }

      /* loop over composition and read data */
            
      if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP]) < VALID_PTR)
        iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */
          
          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Get index */

          n = (long)RDB[nuc + NUCLIDE_Z];
          CheckValue(FUNCTION_NAME, "n", "", n, 1, 110);

          /* Add to density */

          dens[n] = dens[n] + RDB[iso + COMPOSITION_ADENS];

          /* Put pointer */
          
          if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
            nucl[n] = nuc;
          else
            nucl[n] = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_DATA];

          /* Next */

          iso = NextItem(iso);
        }

      /* Check that at least one element was found */

      for (n = 0; n < 110; n++)
        if (nucl[n] > 0)
          break;

      /* Check */

      if ((n == 110) && ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] < VALID_PTR))
        Error(mat, "Material %s has no elements with photon data",
              GetText(mat + MATERIAL_PTR_NAME));

      /* Check if original composition exists */

      if ((long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP] > VALID_PTR)
        {
          /* Pointer to composition */          

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          
          /* Loop over elements*/
          
          for (n = 0; n < 110; n++)
            if (nucl[n] > VALID_PTR)
              {
                /* Check pointer */

                if (nucl[n] != (long)RDB[iso + COMPOSITION_PTR_NUCLIDE])
                  Die(FUNCTION_NAME, "Pointer error");
                
                /* Put pointer density */

                WDB[iso + COMPOSITION_ADENS] = dens[n];

                /* Pointer to next */

                iso = NextItem(iso);
              }
        }
      else
        {
          /* Copy pointer and reset */

          WDB[mat + MATERIAL_PTR_ORIG_NUC_COMP] = RDB[mat + MATERIAL_PTR_COMP]; 
          WDB[mat + MATERIAL_PTR_COMP] = NULLPTR;
          
          /* Put data */
          
          for (n = 0; n < 110; n++)
            if (nucl[n] > VALID_PTR)
              {
                /* Create structure */
                
                iso = NewItem(mat + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);
                
                /* Put pointer and density */
                
                WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)nucl[n];
                WDB[iso + COMPOSITION_ADENS] = dens[n];
              }
        }
      
      /* Next material */

      mat = NextItem(mat);
    }
}

/*****************************************************************************/
