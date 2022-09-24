/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : linkinterfacematerials.c                       */
/*                                                                           */
/* Created:       2018/11/07 (VVa)                                           */
/* Last modified: 2018/11/07 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Links material pointers to interface                         */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "LinkInterfaceMaterials:"

/*****************************************************************************/

void LinkInterfaceMaterials(long ifc)
{
  long marr, mat, nmat, n;

  /***********************************************************************/

  /* Check interface pointer */

  CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

  /* Get material array pointer */

  marr = (long)RDB[ifc + IFC_PTR_MAT_ARR];
  CheckPointer(FUNCTION_NAME, "(marr)", DATA_ARRAY, marr);

  /* Get number of interface materials */

  nmat = (long)RDB[ifc + IFC_N_MAT];

  /* Link materials and set TMS limits */

  for (n = 0; n < nmat; n++)
    {
      /* Loop over materials and find match */

      mat = (long)RDB[DATA_PTR_M0];

      while (mat > VALID_PTR)
        {
          /* Compare name */

          if (CompareStr(mat + MATERIAL_PTR_NAME, marr + n))
            break;

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check if found */

      if (mat < VALID_PTR)
        Error(ifc, "Material %s linked to interface not defined",
              GetText(marr + n));

      /* Material found*/

      /* Link interface and put flag */

      WDB[mat + MATERIAL_PTR_IFC] = (double)ifc;
      WDB[mat + MATERIAL_USE_IFC] = (double)YES;

      WDB[marr + n] = (double)mat;

      /* Next interface material */
    }

  /* Link output material (needed in findinterfaceregion) */

  if ((long)RDB[ifc + IFC_PTR_OUTPUT_MAT] > VALID_PTR)
    {
      /* Is pointer to material name specifically given */

      /* Loop over materials and find match */

      mat = (long)RDB[DATA_PTR_M0];

      while (mat > VALID_PTR)
        {
          /* Compare name */

          if (CompareStr(mat + MATERIAL_PTR_NAME, ifc + IFC_PTR_OUTPUT_MAT))
            break;

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check if found */

      if (mat < VALID_PTR)
        Error(ifc, "Material %s linked to interface not defined",
              GetText(ifc + IFC_PTR_OUTPUT_MAT));

      /* Material found*/

      /* Link material */

      WDB[ifc + IFC_PTR_OUTPUT_MAT] = (double)mat;
    }
  else
    {
      /* No output material separately given, use the first of the interface */
      /* materials. */

      WDB[ifc + IFC_PTR_OUTPUT_MAT] = RDB[marr + 0];
    }

  return;
}
