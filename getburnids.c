/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getburnids.c                                   */
/*                                                                           */
/* Created:       2013/01/24 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Creates a list of library id's and temperatures of nuclides  */
/*              in burnable materials.                                       */
/*                                                                           */
/* Comments: - Muuta ton listan nimi ACTINIDE_ID listiks, tms.               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetburnIDs:"

/*****************************************************************************/

void GetBurnIDs()
{
  long mat, iso, nuc, ptr, TMS;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check burn-flag and division */
          
      if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) &&
          ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] < VALID_PTR))
        {
          /* Loop over composition */
              
          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
              
              /* Check fissile flag (5.1.2014 / 2.1.17) */

              if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE))
                {
                  /* Pointer to next */

                  iso = NextItem(iso);

                  /* Cycle loop */

                  continue;
                }

              /* Get TMS flag */

              TMS = ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS);

              /* Loop over existing */

              ptr = (long)RDB[DATA_PTR_FP_LIB_ID_LIST];
              while (ptr > VALID_PTR)
                {
                  /* Compare library ID's and temperatures */
              
                  if ((RDB[ptr + FP_IDENT_TEMP] == RDB[nuc + NUCLIDE_TEMP])
                      && CompareStr(ptr + FP_IDENT_PTR_ID, nuc + 
                                    NUCLIDE_PTR_LIB_ID) &&
                      ((long)RDB[ptr + FP_IDENT_TMS] == TMS))
                    break;
              
                  /* Next */
              
                  ptr = NextItem(ptr);
                }

              /* Add new if not found */

              if (ptr < VALID_PTR)
                {
                  /* Add item */
                  
                  ptr = NewItem(DATA_PTR_FP_LIB_ID_LIST, FP_IDENT_BLOCK_SIZE);
              
                  /* Put id and temperature */
                  
                  WDB[ptr + FP_IDENT_PTR_ID] = RDB[nuc + NUCLIDE_PTR_LIB_ID];
                  WDB[ptr + FP_IDENT_TEMP] = RDB[nuc + NUCLIDE_TEMP];
                  WDB[ptr + FP_IDENT_TMS] = (double)TMS;
                }

              /* Next nuclide in composition */
      
              iso = NextItem(iso);
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Close list */

  if ((ptr = (long)RDB[DATA_PTR_FP_LIB_ID_LIST]) > VALID_PTR)
    CloseList(ptr);
}

/*****************************************************************************/
