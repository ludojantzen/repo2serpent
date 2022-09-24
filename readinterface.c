/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readinterface.c                                */
/*                                                                           */
/* Created:       2012/02/14 (JLe)                                           */
/* Last modified: 2018/11/05 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads multi-physics interfaces                               */
/*                                                                           */
/* Comments:  - Reading of different interface types is done in their        */
/*              separate subroutines                                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadInterface:"

/*****************************************************************************/

void ReadInterface(long loc0, long update)
{
  long type, ptr, prev;
  FILE *fp;

  /* If reading interfaces from input file */

  if (loc0 < 0)
    fprintf(outp, "Reading multi-physics interface:\n\n");

  /***********************************************************************/

  /***** Common data *****************************************************/

  /* Test file format */

  TestDOSFile(GetText(loc0 + IFC_PTR_INPUT_FNAME));

  /* Open file for reading */

  if ((fp = fopen(GetText(loc0 + IFC_PTR_INPUT_FNAME), "r")) == NULL)
    Error(loc0, "Multi-physics interface file \"%s\" does not exist",
          GetText(loc0 + IFC_PTR_INPUT_FNAME));
  else
    fprintf(outp, "Reading multi-physics interface \"%s\"...\n",
            GetText(loc0 + IFC_PTR_INPUT_FNAME));

  /* Read interface type */

  if (fscanf(fp, "%ld", &type) != 1)
    Die(FUNCTION_NAME, "Could not read interface type from file \"%s\".",
        GetText(loc0 + IFC_PTR_INPUT_FNAME));

  /* Check interface type */

  CheckValue(FUNCTION_NAME, "type", "", type, IFC_TYPE_MIN, IFC_TYPE_MAX);

  /* Close file for now */

  fclose(fp);

  /* Put interface type */

  WDB[loc0 + IFC_TYPE] = (double)type;

  /* Put interface index */

  if(!update)
    {
      if((prev = PrevItem(loc0)) > VALID_PTR)
        {
          WDB[loc0 + IFC_IDX] = RDB[prev + IFC_IDX] + 1.0;
        }
      else
        WDB[loc0 + IFC_IDX] = (double)1;
    }

  switch(type)
    {
    case IFC_TYPE_PT_AVG:

      /* Read interface file with dedicated subroutine */

      ReadIFCPtAvg(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_REG_MESH:
    case IFC_TYPE_REG_MESH_MULTIMAT:
    case IFC_TYPE_REG_MESH_MULTILVL:
    case IFC_TYPE_REG_MESH_GENERAL:

      /* Read interface file with dedicated subroutine */

      ReadIFCRegMesh(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_FUNC:

      /* Read interface file with dedicated subroutine */

      ReadIFCFunc(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_FET_DENSITY:
    case IFC_TYPE_FET_TEMP:

      /* Read interface file with dedicated subroutine */

      ReadIFCFE(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_TET_MESH:

      /* Read interface file with dedicated subroutine */
      /* Not currently up to date (maybe should be removed) */
      /* since there is the OpenFOAM interface              */

      Die(FUNCTION_NAME, "Interface type %ld not currently supported", IFC_TYPE_TET_MESH);

      ReadIFCTetMesh(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_FUEP:

      /* Read interface file with dedicated subroutine */

      ReadIFCFB(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_FPIP:

      /* Read interface file with dedicated subroutine */

      ReadIFCFB(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_OPENFOAM:

      /* Read interface file with dedicated subroutine */

      ReadIFCOFMesh(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_OF_MAT:

      /* Read interface file with dedicated subroutine */

      ReadIFCOFMesh(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_OF_SOLID:

      /* Read interface file with dedicated subroutine */

      ReadIFCOFMesh(loc0, update);

      break;
      /*************************************************/

    default:
      Die(FUNCTION_NAME,
          "Unknown interface type %ld in interface file: %s\n",
          type, GetText(loc0 + IFC_PTR_INPUT_FNAME));

    }

  /* Allocate memory for storing of previous interaction data */

  if (!update)
    {

      /* Allocate memory for previous cell pointer */

      ptr = AllocPrivateData(1, PRIVA_ARRAY);
      WDB[loc0 + IFC_PTR_PREV_CELL] = (double)ptr;

      /* Allocate memory for previous collison */

      AllocValuePair(loc0 + IFC_PTR_PREV_COL_CELL);

    }
}

/*****************************************************************************/
