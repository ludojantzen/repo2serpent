/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processsingleinterface.c                       */
/*                                                                           */
/* Created:       2018/11/29 (VVa)                                           */
/* Last modified: 2019/02/13 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes a single multi-physic interface                    */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessSingleInterface:"

/*****************************************************************************/

void ProcessSingleInterface(long loc0, long updateT, long updateRho)
{
  long type, update;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

  /* TODO: Change all interface processing routines to accept two update */
  /*       parameters */

  if (updateT+updateRho > 0)
    update = 1;

  /* Get interface type */

  type = (long)RDB[loc0 + IFC_TYPE];

  /***** Processing for each interface type ************/

  switch(type)
    {
    case IFC_TYPE_PT_AVG:

      /***** Average of points *************************/

      ProcessIFCPtAvg(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_REG_MESH:
    case IFC_TYPE_REG_MESH_MULTIMAT:
    case IFC_TYPE_REG_MESH_MULTILVL:

      /***** Regular mesh based distribution ***********/

      ProcessIFCRegMesh(loc0, updateT, updateRho);

      break;
      /*************************************************/

    case IFC_TYPE_FUNC:

      /***** Function based interface ******************/

      ProcessIFCFunc(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_FET_DENSITY:
    case IFC_TYPE_FET_TEMP:

      /***** FET-based interface ***********************/


      /* The mat corresponding to IFC_FE_EFFECT_MAT_PTR we can use the provided */
      /* FE coefficients to calculate the density or temperature fields for     */
      /* multiphysics feedback.                                                 */

      ProcessIFCFETMaterial(loc0, update);

      /* Process material specification for outgoing fission energy FET scoring */

      if ((long)RDB[loc0 + IFC_CALC_OUTPUT] == YES)
        {
          /* The mat corresponding to IFC_PTR_MAT is used for configuring the FET */
          /* scoring interface. Otherwise FindInterfaceRegions() will not find    */
          /* any associated region, causing ScoreInterfacePower() to not call any */
          /* Score...() functions---including ScoreFET().                         */

          ProcessIFCFETMaterial(loc0, update);
        }

      break;
      /*************************************************/

    case IFC_TYPE_TET_MESH:

      /***** Unstructured tetrahedral mesh *************/
      /* Can be OpenFOAM based */

      ProcessIFCTetMesh(loc0, update);

      break;
      /*************************************************/

    case IFC_TYPE_FUEP:
    case IFC_TYPE_FPIP:

      /***** Fuel behavior interface *******************/

      ProcessIFCFB(loc0, update);

      break;
      /*************************************************/

    default:
      Die(FUNCTION_NAME,
          "Unknown interface type %ld in interface file: %s\n",
          type, GetText(loc0 + IFC_PTR_INPUT_FNAME));

    }

  /***************************************************************************/

  /***************************************************************************/

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
