/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sendifcoutputdata.c                            */
/*                                                                           */
/* Created:       2019/02/11 (VVa)                                           */
/* Last modified: 2019/02/13 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sends interface output data to Cerberus.                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SendIFCOutputData:"

/*****************************************************************************/

void SendIFCOutputData()
{
  long ifc, i, arr, nval;
  char tmpstr[MAX_STR];

  /* Get name of output field to send */

  SocketReceiveString(tmpstr);

  /* Find the correct interface */

  ifc = (long)RDB[DATA_PTR_IFC0];

  while (ifc > VALID_PTR)
    {
      /* Check if the interface has an output file associated with itself */
      /* If it does, compare the field name */

      if ((long)RDB[ifc + IFC_PTR_POWER_FIELD_NAME] > VALID_PTR)
        if (!strcmp(tmpstr, GetText(ifc + IFC_PTR_POWER_FIELD_NAME)))
            break;

      ifc = NextItem(ifc);
    }

  /* Check if found */

  if (ifc < VALID_PTR)
    Die(FUNCTION_NAME, "Could not find interface providing output field %s",
        tmpstr);

  /* Get output size */

  nval = (long)RDB[ifc + IFC_STAT_NREG];

  /* Send output size */

  SocketSendLong(nval);

  /* Get pointer to relaxed power data */

  arr = (long)RDB[ifc + IFC_PTR_STAT_REL];
  CheckPointer(FUNCTION_NAME, "(arr)", DATA_ARRAY, arr);

  fprintf(outp, "Sending output data for field %s.\n\n", tmpstr);

  /* Send data */

  for (i = 0; i < nval; i++)
    SocketSendDouble(RDB[arr + i]);
 }

/*****************************************************************************/
