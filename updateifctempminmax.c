/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : updateifctempminmax.c                           */
/*                                                                           */
/* Created:       2018/12/05 (RTu)                                           */
/* Last modified: 2019/02/13 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Updates minimum and maximum temperatures of the interface    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UpdateIFCTempMinMax:"

/*****************************************************************************/

void UpdateIFCTempMinMax(long ifc)
{
  long type, tmplist, ntmp, i;
  double Tmin, Tmax, T;

  /* Check interface pointer */

  CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

  /* Reset limiting values */

  Tmax = -INFTY;
  Tmin = INFTY;

  /* Get interface type */

  type = (long)RDB[ifc + IFC_TYPE];

  switch(type)
    {
    case IFC_TYPE_REG_MESH:
    case IFC_TYPE_REG_MESH_MULTIMAT:
    case IFC_TYPE_REG_MESH_MULTILVL:
      {
        /* Get number of temperatures */

        ntmp = (long)RDB[ifc + IFC_NP];

        /* Get pointer to temperature list */

        tmplist = (long)RDB[ifc + IFC_PTR_TMP_LIST];
        CheckPointer(FUNCTION_NAME, "(tmplist)", DATA_ARRAY, tmplist);

        /* Loop over temperatures */

        for (i = 0; i < ntmp; i++)
          {
            /* Get the temperature */

            T = RDB[tmplist + i];

            /* Compare to minimum and maximum */

            if (T < Tmin)
              Tmin = T;

            if (T > Tmax)
              Tmax = T;
          }

        break;
      }
    default:
      Die(FUNCTION_NAME,
          "Interface type %ld is not supported (interface file: %s)",
          type, GetText(ifc + IFC_PTR_INPUT_FNAME));
    }

  /* Update minimum and maximum */

  WDB[ifc + IFC_MIN_TEMP] = Tmin;
  WDB[ifc + IFC_MAX_TEMP] = Tmax;
}
