/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : updateifcdensmax.c                             */
/*                                                                           */
/* Created:       2019/02/12 (RTu)                                           */
/* Last modified: 2019/02/12 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Updates maximum density of the interface                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UpdateIFCDensMax:"

/*****************************************************************************/

void UpdateIFCDensMax(long ifc)
{
  long type, dflist, ndens, i;
  double d, dmax;

  /* Check interface pointer */

  CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

  /* Reset maximum density */

  dmax = 0.0;

  /* Get interface type */

  type = (long)RDB[ifc + IFC_TYPE];

  switch(type)
    {
    case IFC_TYPE_REG_MESH:
    case IFC_TYPE_REG_MESH_MULTIMAT:
    case IFC_TYPE_REG_MESH_MULTILVL:
      {
        /* Get number of densities */

        ndens = (long)RDB[ifc + IFC_NP];

        /* Get pointer to density list */

        dflist = (long)RDB[ifc + IFC_PTR_DF_LIST];
        CheckPointer(FUNCTION_NAME, "(dflist)", DATA_ARRAY, dflist);

        /* Loop over densities */

        for (i = 0; i < ndens; i++)
          {
            /* Get the density */

            d = RDB[dflist + i];

            /* Compare to maximum */

            if (fabs(d) > fabs(dmax))
              dmax = d;
          }

        break;
      }
    default:
      Die(FUNCTION_NAME,
          "Interface type %ld is not supported (interface file: %s)",
          type, GetText(ifc + IFC_PTR_INPUT_FNAME));
    }

  /* Update maximum density */

  WDB[ifc + IFC_MAX_DENSITY] = dmax;
}
