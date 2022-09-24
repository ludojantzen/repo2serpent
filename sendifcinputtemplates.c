/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sendifcinputtemplates.c                        */
/*                                                                           */
/* Created:       2019/02/05 (VVa)                                           */
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

#define FUNCTION_NAME "SendIFCInputTemplates:"

void SendIFCInputMesh(long ifc)
{
  long type, msh;

  /* Get interface type */

  type = (long)RDB[ifc + IFC_TYPE];

  switch(type)
    {
    case IFC_TYPE_REG_MESH:
    case IFC_TYPE_REG_MESH_MULTIMAT:

      /* Get pointer to mesh */

      msh = (long)RDB[ifc + IFC_PTR_SEARCH_MESH_LIST];

      /* Send mesh */

      SendIFCMesh(msh);

      break;

      /*********************************************************************/

    case IFC_TYPE_REG_MESH_MULTILVL:

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

void SendIFCInputTemplates()
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
      /* Need to figure out if something is simply temperature or density? */

      nifc++;

      ifc = NextItem(ifc);
    }

  /* Send two times the number of interfaces */

  SocketSendLong(2*nifc);

  /* Loop over the interfaces to send templates */

  ifc = (long)RDB[DATA_PTR_IFC0];

  while (ifc > VALID_PTR)
    {

      /* Create temperature field name */

      sprintf(tmpstr, "%s", GetText(ifc + IFC_PTR_T_FIELD_NAME));

      /* Send temperature field name */

      SocketSendString(tmpstr);

      /* Send temperature field units */

      SocketSendLong(0);
      SocketSendLong(0);
      SocketSendLong(0);
      SocketSendLong(1);
      SocketSendLong(0);

      /* Send multiplier */

      mul = 1;
      SocketSendDouble(mul);

      /* Send mesh data*/

      SendIFCInputMesh(ifc);

      /* Create density field name */

      sprintf(tmpstr, "%s", GetText(ifc + IFC_PTR_RHO_FIELD_NAME));

      /* Send density field name */

      SocketSendString(tmpstr);

      /* Send density field units */

      SocketSendLong(1);
      SocketSendLong(-3);
      SocketSendLong(0);
      SocketSendLong(0);
      SocketSendLong(0);

      /* Send multiplier */

      mul = -0.001;
      SocketSendDouble(mul);

      /* Send mesh data*/

      SendIFCInputMesh(ifc);

      /* Next interface */

      ifc = NextItem(ifc);
    }

 }

/*****************************************************************************/
