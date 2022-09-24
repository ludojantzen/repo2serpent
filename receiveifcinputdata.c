/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : receiveifcinputdata.c                          */
/*                                                                           */
/* Created:       2019/02/11 (VVa)                                           */
/* Last modified: 2019/03/29 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Receives interface T/rho data from Cerberus.                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReceiveIFCInputData:"

/*****************************************************************************/

void ReceiveIFCInputData()
{
  long ifc, nval, i, arr, temperature, density;
  char tmpstr[MAX_STR];

  /* Get name of field to receive */

  SocketReceiveString(tmpstr);

  /* Initialize some values */

  temperature = 0;
  density = 0;

  /* Find the correct interface and temperature/density */

  ifc = (long)RDB[DATA_PTR_IFC0];

  while (ifc > VALID_PTR)
    {
      /* Compare the field name to interface temperature and density fields */

      if ((long)RDB[ifc + IFC_PTR_T_FIELD_NAME] > VALID_PTR)
        if (!strcmp(tmpstr, GetText(ifc + IFC_PTR_T_FIELD_NAME)))
          {
            temperature = 1;
            density = 0;
            break;
          }

      if ((long)RDB[ifc + IFC_PTR_RHO_FIELD_NAME] > VALID_PTR)
        if (!strcmp(tmpstr, GetText(ifc + IFC_PTR_RHO_FIELD_NAME)))
          {
            temperature = 0;
            density = 1;
            break;
          }

      ifc = NextItem(ifc);
    }

  /* Check if found */

  if (ifc < VALID_PTR)
    Die(FUNCTION_NAME, "Could not find interface for input field %s",
        tmpstr);

  /* Get field size */

  SocketReceiveLong(&nval);

  /* Check field size */

  if (nval != (long)RDB[ifc + IFC_NP])
    Die(FUNCTION_NAME, "Was supposed to receive %ld values for field %s, but IFC "
        " requires %ld values",
        nval, tmpstr, (long)RDB[ifc + IFC_NP]);

  /* Get pointer to correct data array */

  if (density)
    arr = (long)RDB[ifc + IFC_PTR_DF_LIST];
  else if (temperature)
    arr = (long)RDB[ifc + IFC_PTR_TMP_LIST];
  else
    {
      arr = NULLPTR;
      Die(FUNCTION_NAME, "Something strange happened: %ld %ld",
          density, temperature);
    }

  CheckPointer(FUNCTION_NAME, "(arr)", DATA_ARRAY, arr);

  fprintf(outp, "Receiving input data for field %s.\n\n", tmpstr);

  /* Receive data */

  for (i = 0; i < nval; i++)
    SocketReceiveDouble(&WDB[arr + i]);

  /* If we received density data, we'll still need to convert it */

  /* Call process interface for this interface */

  ProcessSingleInterface(ifc, temperature, density);
}


/*****************************************************************************/
