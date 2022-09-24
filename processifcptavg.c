/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processifcptavg.c                              */
/*                                                                           */
/* Created:       2015/01/23 (VVa)                                           */
/* Last modified: 2018/11/07 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes point average multi-physics interfaces             */
/*                                                                           */
/* Comments: - Stats are allocated in allocinterfacestat.c                   */
/*           - Point average interface does not do interpolation wrt. time   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessIFCPtAvg:"

/*****************************************************************************/

void ProcessIFCPtAvg(long ifc, long update)
{
  long loc1, mat, msh, dim, nx, ny, nz, override, nmat, marr, n;
  double x, y, z, rad, xmin, xmax, ymin, ymax, zmin, zmax, dmax;
  double lims[6], matdens;

  /***********************************************************************/

  if (update == (long)NO)
    {
      /* Link materials */

      LinkInterfaceMaterials(ifc);
    }

  /* Set or check TMS limits for linked materials */

  SetIFCTMSLimits(ifc, update);

  /* Avoid compiler warning */

  matdens = 0.0;

  /**********************************/
  /* Process interface density data */
  /**********************************/

  /* Get number of materials in interface */

  nmat = (long)RDB[ifc + IFC_N_MAT];

  /* Process densities only for interfaces with one material */
  /* multi-material interfaces interpret density data as     */
  /* density factors from the beginning */

  if (nmat == 1)
    {
      /* Get material array pointer */

      marr = (long)RDB[ifc + IFC_PTR_MAT_ARR];
      CheckPointer(FUNCTION_NAME, "(marr)", DATA_ARRAY, marr);

      /* Get pointer to the only material */

      mat = (long)RDB[marr];

      /* Get interface maximum density */

      dmax = RDB[ifc + IFC_MAX_DENSITY];

      /* Get material density */

      override = 0;

      if (update)
        {
          /* When the interface is updated, the MATERIAL_ADENS */
          /* Already contains the atomic density */

          if (dmax < 0)
            /* Mass densities given*/
            matdens = -RDB[mat + MATERIAL_MDENS];
          else if(dmax > 0)
            /* Atomic densities given */
            matdens = RDB[mat + MATERIAL_ADENS];
          else
            {
              matdens = 0;
              Error(ifc, "No non-zero densities in distribution");
            }
          fprintf(outp, "Material %s, matdens %E, dmax %E\n",
                  GetText(mat + MATERIAL_PTR_NAME), matdens, dmax);

        }
      else
        {
          /* If this is the first time this interface is processed we will use */
          /* MATERIAL_ADENS for the density regardless of IFC units */
          /* They will be converted to atomic density later         */

          matdens = RDB[mat + MATERIAL_ADENS];

          /* If base material density is given in different units than the */
          /* IFC-density, we cannot use base density as a majorant */
          /* We will just override the base density with the interface density */

          if (matdens*dmax < 0)
            override = 1;
        }

      /* Only increase the material atomic density before calculating majorants */
      /* (not on updates) */

      /* Only change it upwards, this way the normal density card can be used as */
      /* a maximum density, when setting up a coupled calculation */

      if ((fabs(matdens) < fabs(dmax)) || (override == 1))
        {
          /* IFC density greater than material base density */

          if (!update)
            {
              WDB[mat + MATERIAL_ADENS] = dmax;
              matdens = RDB[mat + MATERIAL_ADENS];
            }
          else
            Die(FUNCTION_NAME, "Material density exceeds material majorant for %s",
                GetText(mat + MATERIAL_PTR_NAME));
        }
      else
        {
          /* Material base density greater than IFC density */

          dmax = matdens;

        }

      /* Loop over points to set density factors */

      loc1 = (long)RDB[ifc + IFC_PTR_POINTS];

      while (loc1 > VALID_PTR)
        {

          /* Convert densities to density factors */

          if (RDB[loc1 + IFC_PT_DF]*dmax < 0.0)
            Error(ifc, "Inconsistent densities given in distribution");
          else
            WDB[loc1 + IFC_PT_DF] = RDB[loc1 + IFC_PT_DF]/dmax;

          loc1 = NextItem(loc1);
        }
    }

  /***** Create the mesh that the points lie on *******************************/

  if (update == (long)NO)
    {
      /* Get limits */

      xmin = RDB[ifc + IFC_MESH_XMIN];
      xmax = RDB[ifc + IFC_MESH_XMAX];
      ymin = RDB[ifc + IFC_MESH_YMIN];
      ymax = RDB[ifc + IFC_MESH_YMAX];
      zmin = RDB[ifc + IFC_MESH_ZMIN];
      zmax = RDB[ifc + IFC_MESH_ZMAX];

      /* Get exclusion radius and dimension */

      rad = RDB[ifc + IFC_EXCL_RAD];
      dim = (long)RDB[ifc + IFC_DIM];

      /* Calculate search mesh size */

      if ((nx = (long)(xmax - xmin)/(2.0*rad)) < 1)
        nx = 1;

      if ((ny = (long)(ymax - ymin)/(2.0*rad)) < 1)
        ny = 1;

      if ((nz = (long)(zmax - zmin)/(2.0*rad)) < 1)
        nz = 1;

      /* Adjust boundaries */

      xmin = xmin - 1E-6;
      xmax = xmax + 1E-6;
      ymin = ymin - 1E-6;
      ymax = ymax + 1E-6;
      zmin = zmin - 1E-6;
      zmax = zmax + 1E-6;

      /* 1D and 2D distributions */

      if (dim == 1)
        {
          xmin = -INFTY;
          xmax = INFTY;
          nx = 1;

          ymin = -INFTY;
          ymax = INFTY;
          ny = 1;
        }
      else if (dim == 2)
        {
          zmin = -INFTY;
          zmax = INFTY;
          nz = 1;
        }

      /* Put mesh variables */

      lims[0] = xmin;
      lims[1] = xmax;
      lims[2] = ymin;
      lims[3] = ymax;
      lims[4] = zmin;
      lims[5] = zmax;

      /* Create mesh structure */

      msh = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_PTR, -1,
                       50, 50, 50, lims, -1);

      /* Put pointer */

      WDB[ifc + IFC_PTR_SEARCH_MESH_LIST] = (double)msh;

      /* Loop over points */

      loc1 = (long)RDB[ifc + IFC_PTR_POINTS];
      while (loc1 > VALID_PTR)
        {
          /* Get coordinates */

          x = RDB[loc1 + IFC_PT_X];
          y = RDB[loc1 + IFC_PT_Y];
          z = RDB[loc1 + IFC_PT_Z];

          /* Add point in search mesh */

          AddSearchMesh(msh, loc1, x - rad, x + rad, y - rad, y + rad,
                        z - rad, z + rad);

          /* Next point */

          loc1 = NextItem(loc1);
        }
    }

  /********************************************************************/

  /* Print out basic information about the interface */

  fprintf(outp, "\nPoint average based interface \"%s\":\n\n",
          GetText(ifc + IFC_PTR_INPUT_FNAME));

  /* Print outp the maximum density of the interface */

  fprintf(outp, " - Maximum density of interface: %E\n",
          RDB[ifc + IFC_MAX_DENSITY]);

  /***** Print materials in the interface *****/

  /* Get number of materials */

  nmat = (long)RDB[ifc + IFC_N_MAT];

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

  /* Get interface maximum density */

  dmax = RDB[ifc + IFC_MAX_DENSITY];

  for (n = 0; n < nmat; n++)
    {
      /* Get pointer to material */

      mat = (long)RDB[marr + n];
      CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

      /* Get material density */

      if (update == (long)YES)
        {
          /* Get mass density or atomic density based on which one was */
          /* given */

          if (dmax < 0)
            matdens = -RDB[mat + MATERIAL_MDENS];
          else if(dmax > 0)
            matdens = RDB[mat + MATERIAL_ADENS];
        }
      else
        matdens = RDB[mat + MATERIAL_ADENS];

      /* Print outp the density used for calculating the majorant */

      fprintf(outp, " - Majorant density of material %s: %E\n",
              GetText(mat + MATERIAL_PTR_NAME), matdens);

      /* Print material temperature TMS limits */

      if (RDB[mat + MATERIAL_TMS_MODE] == (double)NO)
        fprintf(outp, " - No TMS treatment for material %s\n",
                GetText(mat + MATERIAL_PTR_NAME));
      else
        fprintf(outp, " - Material %s TMS limits: %.2f K and %.2f K\n",
                GetText(mat + MATERIAL_PTR_NAME),
                RDB[mat + MATERIAL_TMS_TMIN], RDB[mat + MATERIAL_TMS_TMAX]);
    }


  /********************************************************************/
}

/*****************************************************************************/
