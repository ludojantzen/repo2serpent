/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processifcregmesh.c                            */
/*                                                                           */
/* Created:       2015/02/02 (VVa)                                           */
/* Last modified: 2019/02/13 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes regular mesh based multi-physics interfaces        */
/*                                                                           */
/* Comments: - Stats are allocated in allocinterfacestat.c                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessIFCRegMesh:"

#define MAX_REMAP_CELLS 20

/*****************************************************************************/

void ProcessIFCRegMesh(long ifc, long updateT, long updateRho)
{
  long mat, np, n, nmat, marr;
  long tmplist, tmplist0, dflist, dflist0;
  double dmax, matdens;

  /***********************************************************************/

  /***** Link material ***************************************************/

  mat = -1;

  /* Avoid compiler warning */

  matdens = 0.0;

  /* Link materials */

  if (updateT + updateRho == 0)
    LinkInterfaceMaterials(ifc);

  /* Update maximum and minimum temperatures of the interface */
  /* This actually only necessary for internal ifcs */
  /* For external ifcs max and min are updated when reading the ifc file */

  if (updateT == YES)
    UpdateIFCTempMinMax(ifc);

  /* Set or check TMS limits for linked materials */

  if ((updateT + updateRho == 0) || (updateT == YES))
    SetIFCTMSLimits(ifc, updateT);

  /**********************************/
  /* Process interface density data */
  /**********************************/

  /* Update interface maximum density */
  /* This actually only necessary for internal ifcs */
  /* For external ifcs max is updated when reading the ifc file */

  if (updateRho == YES)
    UpdateIFCDensMax(ifc);

  /* Get pointer to density list */

  dflist = (long)RDB[ifc + IFC_PTR_DF_LIST];

  /* Get size */

  np = (long)RDB[ifc + IFC_NP];
  CheckValue(FUNCTION_NAME, "np", "", np, 1, 100000000);

  /* Get interface maximum density */

  dmax = RDB[ifc + IFC_MAX_DENSITY];

  /* Get number of materials in interface */

  nmat = (long)RDB[ifc + IFC_N_MAT];

  /* Check if density was updated or this is initialization */

  if ((updateRho == YES) || (updateRho + updateT == 0))
    {
      if (nmat == 1)
        {
          /* Get material array pointer */

          marr = (long)RDB[ifc + IFC_PTR_MAT_ARR];
          CheckPointer(FUNCTION_NAME, "(marr)", DATA_ARRAY, marr);

          /* Get pointer to the only material */

          mat = (long)RDB[marr];

          /* Process depending whether this is an update or not */

          if (updateRho == YES)
            {
              /* When the interface is updated, the MATERIAL_ADENS */
              /* already contains the atomic density */

              if (dmax < 0.0)
                {
                  /* Mass densities given */

                  matdens = -RDB[mat + MATERIAL_MDENS];

                  /* If this function was called from ProcessInternalIFCs() after */
                  /* InitInternal() in main MATERIAL_MDENS has not been set yet */

                  if (matdens == 0.0)
                    matdens = RDB[mat + MATERIAL_ADENS];
                }
              else if (dmax > 0.0)
                {
                  /* Atomic densities given */

                  matdens = RDB[mat + MATERIAL_ADENS];
                }
              else
                {
                  matdens = 0.0;
                  Error(ifc, "No non-zero densities in distribution");
                }

              /* Check maximum density against majorant */

              if (fabs(matdens) < fabs(dmax))
                Die(FUNCTION_NAME, "Material density exceeds material majorant for %s: %f vs %f",
                    GetText(mat + MATERIAL_PTR_NAME), dmax, matdens);
            }
          else
            {
              /* If this is the first time this interface is processed we will use */
              /* MATERIAL_ADENS for the density regardless of IFC units */
              /* They will be converted to atomic density later         */

              matdens = RDB[mat + MATERIAL_ADENS];

              /* Only increase the material atomic density before calculating majorants */
              /* (not on updates) */

              /* Only change it upwards, this way the normal density card can be used as */
              /* a maximum density, when setting up a coupled calculation */

              /* If base material density is given in different units than the */
              /* IFC-density, we cannot use base density as a majorant */
              /* We will just override the base density with the interface density */

              if ((fabs(matdens) < fabs(dmax)) || (matdens*dmax < 0.0))
                {
                  /* IFC density greater than material base density */

                  WDB[mat + MATERIAL_ADENS] = dmax;
                  matdens = dmax;
                }
            }

          /* Convert densities to density factors */

          if (matdens == 0.0)
            Die(FUNCTION_NAME, "Material majorant density is zero for %s",
                GetText(mat + MATERIAL_PTR_NAME));

          for (n = 0; n < np; n++)
            {
              if (RDB[dflist + n]*matdens < 0.0)
                Error(ifc, "Inconsistent densities given in distribution %f %f", RDB[dflist + n], matdens);
              else
                WDB[dflist + n] = RDB[dflist + n]/matdens;
            }
        }
      else
        {
          /* Check that DF_LIST constains density factors */

          for (n = 0; n < np; n++)
            {
              if (RDB[dflist + n] < 0.0 || RDB[dflist + n] > 1.0)
              {
                Error(ifc, "Density must be given as density factors for multimaterial interface. Found %f",
                      RDB[dflist + n]);
              }
            }

        }
    }

  /**********************************/
  /* Print information              */
  /**********************************/

  fprintf(outp, "\nRegular mesh based interface \"%s\":\n\n",
          GetText(ifc + IFC_PTR_INPUT_FNAME));

  /* Print outp the maximum density of the interface */

  fprintf(outp, " - Maximum density of interface: %E\n",
          dmax);

  /***** Print materials in the interface *****/

  /* Get pointer to array */

  marr = (long)RDB[ifc + IFC_PTR_MAT_ARR];
  CheckPointer(FUNCTION_NAME, "(marr)", DATA_ARRAY, marr);

  fprintf(outp, " - Interface is linked to the following materials:\n");

  for (n = 0; n < nmat; n++)
    {
      /* Get pointer to material */

      mat = (long)RDB[marr + n];
      CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

      fprintf(outp, "   * %s\n", GetText(mat + MATERIAL_PTR_NAME));
    }

  /****** Print limits for different materials ******/

  for (n = 0; n < nmat; n++)
    {
      /* Get pointer to material */

      mat = (long)RDB[marr + n];
      CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

      /* Get material density */

      if (updateT + updateRho != 0)
        {
          /* Get mass density or atomic density based on which one was */
          /* given */

          if (dmax < 0)
            {
              matdens = -RDB[mat + MATERIAL_MDENS];

              /* If this function was called from ProcessInternalIFCs() after */
              /* InitInternal() in main MATERIAL_MDENS has not been set yet */

              if (matdens == 0.0)
                matdens = RDB[mat + MATERIAL_ADENS];
            }
          else if(dmax > 0)
            matdens = RDB[mat + MATERIAL_ADENS];
        }
      else
        matdens = RDB[mat + MATERIAL_ADENS];

      /* Print outp the density used for calculating the majorant */

      fprintf(outp, " - Majorant density of material %s: %E\n",
              GetText(mat + MATERIAL_PTR_NAME), matdens);

      /* Print material temperature TMS limits */

      if ((long)RDB[mat + MATERIAL_TMS_MODE] == NO)
        fprintf(outp, " - No TMS treatment for material %s\n",
                GetText(mat + MATERIAL_PTR_NAME));
      else
        fprintf(outp, " - Material %s TMS limits: %.2f K and %.2f K\n",
                GetText(mat + MATERIAL_PTR_NAME),
                RDB[mat + MATERIAL_TMS_TMIN], RDB[mat + MATERIAL_TMS_TMAX]);
    }

  /* In time dependent coupled calculation, we'll create separate lists for BOI */
  /* To allow interpolation of temperature/density data inside the interval     */

  if ((RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN) && (updateRho + updateT == 0))
    {
      /* Get pointer to temperature list */

      tmplist = (long)RDB[ifc + IFC_PTR_TMP_LIST];

      /* Allocate memory for temperatures */

      tmplist0 = ReallocMem(DATA_ARRAY, np);

      /* Put temperature list pointer to memory */

      WDB[ifc + IFC_PTR_TMP_LIST_BOI] = (double)tmplist0;

      /* Copy temperature data from EOI to BOI */

      memcpy(&WDB[tmplist0], &RDB[tmplist], np*sizeof(double));

      /* Allocate memory for densities */

      dflist0 = ReallocMem(DATA_ARRAY, np);

      /* Put density list pointer to memory */

      WDB[ifc + IFC_PTR_DF_LIST_BOI] = (double)dflist0;

      /* Copy density data from EOI to BOI */

      memcpy(&WDB[dflist0], &RDB[dflist], np*sizeof(double));
    }

  return;

}
