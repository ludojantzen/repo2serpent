/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sendifcoutputtemplates.c                       */
/*                                                                           */
/* Created:       2019/02/06 (VVa)                                           */
/* Last modified: 2019/03/29 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sends interface templates to Cerberus.                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "krakenmeshtypes.h"

#define FUNCTION_NAME "SendIFCOutputTemplates:"

void SendIFCOutputMesh(long ifc)
{
  long type, msh;

  /* Get interface type */

  type = (long)RDB[ifc + IFC_TYPE];

  switch(type)
    {
    case IFC_TYPE_REG_MESH:
    case IFC_TYPE_REG_MESH_MULTIMAT:

      if (RDB[ifc + IFC_OUTPUT_TYPE] == IFC_OUTPUT_SAME_MESH)
        {
          /* Get pointer to mesh */

          msh = (long)RDB[ifc + IFC_PTR_SEARCH_MESH_LIST];

          SendIFCMesh(msh);
        }
      else
        {
          Die(FUNCTION_NAME,
              "IFC output mesh sender not implemented for non-same mesh reg mesh");

        }
      break;

      /*********************************************************************/

    case IFC_TYPE_REG_MESH_MULTILVL:

      if (RDB[ifc + IFC_OUTPUT_TYPE] != IFC_OUTPUT_SAME_MESH)
        Die(FUNCTION_NAME, "Don't know how to send non same-mesh output field "
            "for multi-level-meshes.");

      /* Send mesh type (nested mesh) */

      SocketSendLong(KRAKEN_MESH_TYPE_NESTED);

      /* Send number of nested meshes */

      SocketSendLong((long)RDB[ifc + IFC_N_MESH_LVL]);

      /* Get pointer to first mesh */

      msh = (long)RDB[ifc + IFC_PTR_SEARCH_MESH_LIST];

      /* Loop over meshes and send each in its turn */

      while (msh > VALID_PTR)
        {
          SendIFCMesh(msh);

          msh = NextItem(msh);
        }

      break;

      /*********************************************************************/

    default:
      Die(FUNCTION_NAME, "IFC input mesh sender not implemented for type %ld",
          type);

      break;

    }

  /* Send integer to indicate no remapping */

  SocketSendLong(0);
}
/*****************************************************************************/

void SendIFCOutputTemplates()
{
  long ifc, nifc;
  double mul;
  char tmpstr[MAX_STR];

  /* In the case of no interface definitions, simply send zero */

  if ((ifc = (long)RDB[DATA_PTR_IFC0]) < VALID_PTR)
    {
      SocketSendLong(0);

      return;
    }

  /* Count the number of interfaces */

  nifc = 0;

  while (ifc > VALID_PTR)
    {
      /* Check if the interface has an output file associated with itself */

      if ((long)RDB[ifc + IFC_PTR_POWER_FIELD_NAME] > VALID_PTR)
        nifc++;

      ifc = NextItem(ifc);
    }

  /* Send the number of interfaces providing output */

  SocketSendLong(nifc);

  /* Loop over the interfaces to send templates */

  ifc = (long)RDB[DATA_PTR_IFC0];

  while (ifc > VALID_PTR)
    {

      if ((long)RDB[ifc + IFC_PTR_POWER_FIELD_NAME] > VALID_PTR)
        {
          /* Create power field name */

          sprintf(tmpstr, "%s", GetText(ifc + IFC_PTR_POWER_FIELD_NAME));

          /* Send power field name */

          SocketSendString(tmpstr);

          /* Send power field units (Watts integral) */

          SocketSendLong(1);
          SocketSendLong(2);
          SocketSendLong(-3);
          SocketSendLong(0);
          SocketSendLong(0);

          /* Send multiplier */

          mul = 1;
          SocketSendDouble(mul);

          /* Send mesh data*/

          SendIFCOutputMesh(ifc);
        }

      /* Next interface */

      ifc = NextItem(ifc);
    }

 }

/*****************************************************************************/
