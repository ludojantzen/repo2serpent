/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : clearreltransmuxs.c                            */
/*                                                                           */
/* Created:       2014/07/23 (VVa)                                           */
/* Last modified: 2015/06/24 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Clears relaxed one-group cross sections from material dep-   */
/*              lists between burnup iterations                              */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ClearRelTransmuXS:"

/*****************************************************************************/

void ClearRelTransmuXS()
{
  long dep, mat;

  /* Check burnup mode */

  if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO) 
    return;

  /* Loop over materials */
    
  mat = (long)RDB[DATA_PTR_M0];

  while (mat > VALID_PTR)
    {
      /* Check burn flag and test parallel id's */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        {
          /**** Clear flux ****/

          WDB[mat + MATERIAL_BURN_FLUX_REL] = 0.0;

          /**** Clear transmutation list ****/

          /* Pointer to transmutation list (pointteri voi olla null, jos */
          /* on vain decay-askelia, silloin ei tarvitse nollata) */
  
          dep = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST];

          /* Loop over list */

          while (dep > VALID_PTR)
            {

              /* Reset relaxed xs */

              WDB[dep + DEP_TRA_REL] = 0.0;

              /* Next reaction */

              dep = NextItem(dep);
            }

          /**** Clear fission list ****/

          /* Pointer to fission list (pointteri voi olla null jos materiaali */
          /* on ei-fissiili, silloin ei tarvitse nollata) */

          dep = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST];
  
          /* Loop over list */

          while (dep > VALID_PTR)
            {

              /* Reset relaxed xs */

              WDB[dep + DEP_TRA_REL] = 0.0;

              /* Next reaction */

              dep = NextItem(dep);
            }

        }

      /* Next material */

      mat = NextItem(mat);

    }

}
