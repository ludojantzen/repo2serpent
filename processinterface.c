/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processinterface.c                             */
/*                                                                           */
/* Created:       2012/02/14 (JLe)                                           */
/* Last modified: 2019/02/06 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes multi-physics interfaces                           */
/*                                                                           */
/* Comments: - Stats are allocated in allocinterfacestat.c                   */
/*                                                                           */
/*           - Polttoaineinterfacen aksiaalijako lisätty 3.4.2013            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessInterface:"

/*****************************************************************************/

void ProcessInterface(long update)
{
  long ifc, type, mat, mat0, beg, end, i;
  double Tdop;
  char tmpstr[MAX_STR], fieldname[MAX_STR];

  /* Check pointer */

  if ((ifc = (long)RDB[DATA_PTR_IFC0]) < VALID_PTR)
    return;

  fprintf(outp, "Processing multi-physics interfaces...\n");

  /* Loop over interfaces */

  while (ifc > VALID_PTR)
    {

      /* Get interface type */

      type = (long)RDB[ifc + IFC_TYPE];

      /***** Processing for each interface type ************/

      switch(type)
        {
        case IFC_TYPE_PT_AVG:

          /***** Average of points *************************/

          ProcessIFCPtAvg(ifc, update);

          break;
          /*************************************************/

        case IFC_TYPE_REG_MESH:
        case IFC_TYPE_REG_MESH_MULTIMAT:
        case IFC_TYPE_REG_MESH_MULTILVL:

          /***** Regular mesh based distribution ***********/

          ProcessIFCRegMesh(ifc, update, update);

          break;
          /*************************************************/

        case IFC_TYPE_FUNC:

          /***** Function based interface ******************/

          ProcessIFCFunc(ifc, update);

          break;
          /*************************************************/

        case IFC_TYPE_FET_DENSITY:
        case IFC_TYPE_FET_TEMP:

          /***** FET-based interface ***********************/

          ProcessIFCFETMaterial(ifc, update);

          break;
          /*************************************************/

        case IFC_TYPE_TET_MESH:

          /***** Unstructured tetrahedral mesh *************/
          /* Can be OpenFOAM based */

          ProcessIFCTetMesh(ifc, update);

          break;
          /*************************************************/

        case IFC_TYPE_FUEP:
        case IFC_TYPE_FPIP:

          /***** Fuel behavior interface *******************/

          ProcessIFCFB(ifc, update);

          break;
          /*************************************************/

        default:
          Die(FUNCTION_NAME,
              "Unknown interface type %ld in interface file: %s\n",
              type, GetText(ifc + IFC_PTR_INPUT_FNAME));

        }

      /* Next interface */

      ifc = NextItem(ifc);
    }

  /***************************************************************************/

  /***** Common stuff (JLe 1.10.2015 / 2.1.25) *******************************/

  if (!update)
    {

      /***********************/
      /* Set interface names */
      /***********************/

      /* Loop over interfaces */

      ifc = (long)RDB[DATA_PTR_IFC0];

      while (ifc > VALID_PTR)
        {
          /* Use the file name base as the name for the interface */

          sprintf(tmpstr, "%s", GetText(ifc + IFC_PTR_INPUT_FNAME));

          /* Find the last slash '/' from the file name */

          beg = 0;

          for (i = 0; i < strlen(tmpstr); i++)
            if (tmpstr[i] == '/')
              beg = i+1;

          /* Find the first dot '.' after beg from the file name */

          end = strlen(tmpstr)-1;

          for (i = strlen(tmpstr)-1; i > beg; i--)
            if (tmpstr[i] == '.')
              end = i;

          /* Copy the name */

          memcpy(tmpstr, GetText(ifc + IFC_PTR_INPUT_FNAME)+beg, (end - beg)*sizeof(char));

          /* Null terminate the string */

          tmpstr[end-beg] = '\0';

          /* Store names */

          WDB[ifc + IFC_PTR_NAME] = PutText(tmpstr);

          sprintf(fieldname, "%s_temperature", tmpstr);
          WDB[ifc + IFC_PTR_T_FIELD_NAME] = PutText(fieldname);

          sprintf(fieldname, "%s_density", tmpstr);
          WDB[ifc + IFC_PTR_RHO_FIELD_NAME] = PutText(fieldname);

          /******************************************************************/
          /* Use the output file name base as the name for the output field */
          /******************************************************************/

          if ((long)RDB[ifc + IFC_PTR_OUTPUT_FNAME] > VALID_PTR)
            {
              sprintf(tmpstr, "%s", GetText(ifc + IFC_PTR_OUTPUT_FNAME));

              /* Find the last slash '/' from the file name */

              beg = 0;

              for (i = 0; i < strlen(tmpstr); i++)
                if (tmpstr[i] == '/')
                  beg = i+1;

              /* Find the first dot '.' after beg from the file name */

              end = strlen(tmpstr)-1;

              for (i = strlen(tmpstr)-1; i > beg; i--)
                if (tmpstr[i] == '.')
                  end = i;

              /* Copy the name */

              memcpy(tmpstr, GetText(ifc + IFC_PTR_OUTPUT_FNAME)+beg, (end - beg)*sizeof(char));

              /* Null terminate the string */

              tmpstr[end-beg] = '\0';

              /* Store field name */

              WDB[ifc + IFC_PTR_POWER_FIELD_NAME] = PutText(tmpstr);
            }

          /* Next interface */

          ifc = NextItem(ifc);
        }

      /********************************************/
      /* Copy interface data to divided materials */
      /********************************************/

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {

          /* Check if doppler preprocessor should be set on based on interface */
          /* distribution */

          if ((Tdop = RDB[mat + MATERIAL_DOPPLER_TEMP]) < -1.0)
            {

              /* Remove minus sign */

              Tdop = -Tdop;

              /* Check Doppler temperature */

              CheckValue(FUNCTION_NAME, "Tdop", "", Tdop, 0.0, 100000.0);

              /* Store Doppler temperature */

              WDB[mat + MATERIAL_DOPPLER_TEMP] = Tdop;
            }

          /* Copy pointer from parent to children if it exists */

          if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
            if ((long)RDB[mat0 + MATERIAL_PTR_IFC] > VALID_PTR)
              {
                fprintf(outp, "Copying interface information from %s to %s\n",
                        GetText(mat0 + MATERIAL_PTR_NAME),
                        GetText(mat + MATERIAL_PTR_NAME));

                /* Copy pointer */

                WDB[mat + MATERIAL_PTR_IFC] = RDB[mat0 + MATERIAL_PTR_IFC];

                /* Copy mode and limits */

                WDB[mat + MATERIAL_TMS_TMIN] = RDB[mat0 + MATERIAL_TMS_TMIN];
                WDB[mat + MATERIAL_TMS_TMAX] = RDB[mat0 + MATERIAL_TMS_TMAX];
                WDB[mat + MATERIAL_TMS_MODE] = RDB[mat0 + MATERIAL_TMS_MODE];

                /* Copy material density (may be increased based on interface) */

                WDB[mat + MATERIAL_ADENS] = WDB[mat0 + MATERIAL_ADENS];

                /* Copy Doppler-preprocessor temperature */

                WDB[mat + MATERIAL_DOPPLER_TEMP] = RDB[mat0 + MATERIAL_DOPPLER_TEMP];

                /* Copy knowledge of interface */

                WDB[mat + MATERIAL_USE_IFC] = RDB[mat0 + MATERIAL_USE_IFC];

              }

          /* Check if material has tms-limits defined but no interface data */

          if (RDB[mat + MATERIAL_TMS_TMIN] < RDB[mat + MATERIAL_TMS_TMAX])
            if (RDB[mat + MATERIAL_USE_IFC] != (double)YES)
              Warn(FUNCTION_NAME, "TMS-temperature limits set for material %s, but no interface linked", GetText(mat + MATERIAL_PTR_NAME));

          /* Next material */

          mat = NextItem(mat);
        }
    }
  /***************************************************************************/

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
