/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processifctetmesh.c                            */
/*                                                                           */
/* Created:       2015/01/15 (VVa)                                           */
/* Last modified: 2018/11/07 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes tetra mesh based multi-physics interfaces          */
/*                                                                           */
/* Comments: - Stats are allocated in allocinterfacestat.c                   */
/*           - Split from processinterface.c                                 */
/*           - TODO: Loop over tet cells to set materials and densities not  */
/*             working                                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessIFCTetMesh:"

/*****************************************************************************/

void ProcessIFCTetMesh(long ifc, long update)
{
  long loc1, mat, ptr, msh, np, idx;
  long cell, count, nc, ncc, tmplist, dflist;
  long tmplist0, dflist0, prev, tetlist, tet, ntet, i;
  double xmin, xmax, ymin, ymax, zmin, zmax, dmax;
  double last, frac, matdens, d, T, df, dfmin, dfmax;
  double lims[6], bb[6];

  /***********************************************************************/

  mat = -1;

  matdens = 0.0;

  /* Get pointer to temperature list */

  tmplist = (long)RDB[ifc + IFC_PTR_TMP_LIST];

  /* Get pointer to density list */

  dflist = (long)RDB[ifc + IFC_PTR_DF_LIST];

  /* Get pointer to tet list */

  tetlist = (long)RDB[ifc + IFC_PTR_TET_LIST];

  /* First do things that are only done the first time */

  if (!update)
    {
      /*************************************************************************/
      /***** Link materials, interface and geometry cells **********************/
      /*************************************************************************/

      /***** Link material to interface and interface to material **************/
      /***** for tet mesh with a single material                  **************/

      if ((long)RDB[ifc + IFC_PTR_OF_MFILE] < VALID_PTR)
        {
          /* Link interface to material and vice versa */

          LinkInterfaceMaterials(ifc);

          /* Get material array pointer */

          mat = (long)RDB[ifc + IFC_PTR_MAT_ARR];
          CheckPointer(FUNCTION_NAME, "(mat_arr)", DATA_ARRAY, mat);

          /* Get pointer to the material that was linked */

          mat = (long)RDB[mat];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Get pointer to linked geometry cell */

          cell = (long)RDB[ifc + IFC_PTR_GCELL_LIST];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Put material pointer there too */

          WDB[cell + CELL_PTR_MAT] = (double)mat;

        }
      else
        {
          /* An interface with multiple materials */
          /* Link materials to geometry cells by looping over tet-cells */

          loc1 = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];

          /* Reset previous material */

          prev = -1;

          while (loc1 > VALID_PTR)
            {
              /* Get geometry cell */

              cell = (long)RDB[loc1 + IFC_TET_PRNT_PTR_CELL];
              CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

              /* Negative pointer for CELL_PTR_MAT means it is a */
              /* pointer to material name (not yet linked) */

              if (RDB[cell + CELL_PTR_MAT] < 0)
                {
                  /* Not yet linked */

                  /* Remove minus sign */

                  WDB[cell + CELL_PTR_MAT] = -RDB[cell + CELL_PTR_MAT];

                  /* Reset found */

                  mat = -1;

                  /* Check previous */

                  if (prev > VALID_PTR)
                    {
                      /* Check direct match with previous */

                      if (CompareStr(cell + CELL_PTR_MAT, prev + MATERIAL_PTR_NAME))
                        mat = prev;

                    }

                  /* Check match */

                  if (mat < VALID_PTR)
                    {
                      /* Find this material from global materials */

                      mat = (long)RDB[DATA_PTR_M0];
                      if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, GetText(cell + CELL_PTR_MAT)))
                          < VALID_PTR)
                        Error(ifc, "Material %s is not defined", GetText(cell + CELL_PTR_MAT));
                    }

                  /* Check match */

                  if (mat < VALID_PTR)
                    Error(ifc, "Material %s is not defined", GetText(cell + CELL_PTR_MAT));

                  /* Set pointer */

                  WDB[cell + CELL_PTR_MAT] = (double)mat;

                  /* Keep pointer for quickly checking the previous material */

                  prev = mat;

                }

              /* Next tet-cell */

              loc1 = NextItem(loc1);
            }
        }

      /* Loop over tet cells to modify material densities and temperatures */
      /* e.g. majorant density and TMS-limits */

      loc1 = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];

      while (loc1 > VALID_PTR)
        {
          /* Get geometry cell */

          cell = (long)RDB[loc1 + IFC_TET_PRNT_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Get material from cell or link material to cell */

          mat = (long)RDB[cell + CELL_PTR_MAT];

          /* Get tet-cell index */

          idx = (long)RDB[loc1 + IFC_TET_PRNT_IDX];

          /******************************************************************/
          /* Check that density format corresponds with material definition */
          /* and update material maximum density if needed                  */
          /******************************************************************/

          /* If no density file is given, we don't have to check this */

          if (strcmp(GetText(ifc + IFC_PTR_OF_RFILE), "-1"))
            {
              /* Get tet-cell density */

              d = RDB[dflist + idx];

              /* Check that the densities are consistent (mass/atomic) */

              if (RDB[mat + MATERIAL_ADENS]*d < 0)
                Error(ifc, "Inconsistent densities (mass/atomic) given in distribution:\n"
                      "Density given in material card: %E\n"
                      "Density given for interface cell: %E");

              /**************************************************************/
              /* Increase the material density before calculating majorants */
              /* (not on updates)                                           */
              /**************************************************************/

              /* Only change it upwards, this way the normal density card can be used as */
              /* a maximum density, when setting up a coupled calculation */

              /* Get material density */

              /* As this is the first time this interface is processed we will use */
              /* MATERIAL_ADENS for the density regardless of IFC units */
              /* They will be converted to atomic density later         */

              matdens = RDB[mat + MATERIAL_ADENS];

              if (fabs(matdens) < fabs(d))
                {
                  /* IFC density greater than material base density */
                  /* Increase the material density */

                  WDB[mat + MATERIAL_ADENS] = d;
                }
            }

          /****************************************************************/
          /* Modify material TMS limits if interface contains temperature */
          /* information                                                  */
          /****************************************************************/

          if (RDB[ifc + IFC_MIN_TEMP] != RDB[ifc + IFC_MAX_TEMP])
            {
              /* Get tet-cell temperature */

              T = RDB[tmplist + idx];

              /* Compare to material maximum and minimum */

              if (T > RDB[mat + MATERIAL_TMS_TMAX])
                WDB[mat + MATERIAL_TMS_TMAX] = T;

              /* Put minimum temperature */

              if (T < RDB[mat + MATERIAL_TMS_TMIN])
                WDB[mat + MATERIAL_TMS_TMIN] = T;

            }

          loc1 = NextItem(loc1);
        }

      /******************************************************/
      /* Loop over GEOMETRY cells to set material TMS modes */
      /******************************************************/

      /* Get pointer to linked geometry cells */

      cell = (long)RDB[ifc + IFC_PTR_GCELL_LIST];
      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

      /* Loop over cell list */

      while (cell > VALID_PTR)
        {

          /* Get pointer to cell material */

          mat = (long)RDB[cell + CELL_PTR_MAT];

          /* Check that material doppler temperature is not already */
          /* set by tmp-card */

          if (RDB[mat + MATERIAL_DOPPLER_TEMP] >= 0)
            Error(mat, "Material temperature set by tmp-card but a temperature "
                  "distribution is also given by interface %s",
                  GetText(ifc + IFC_PTR_INPUT_FNAME));

          /* Set on-the-fly Doppler-broadening mode */
          /* or Doppler preprocessor temperature    */

          if ((RDB[mat + MATERIAL_TMS_TMIN] <
               RDB[mat + MATERIAL_TMS_TMAX]) ||
              ((long)RDB[mat + MATERIAL_PTR_MIX] > VALID_PTR))
            {
              /* Set TMS mode on */

              WDB[mat + MATERIAL_TMS_MODE] = (double)YES;

              /* Set Doppler-preprocessor off */

              WDB[mat + MATERIAL_DOPPLER_TEMP] = -1.0;

            }
          else
            {
              /* No tms-limits given in material card and constant temperature */
              /* distribution in interface -> Use preprocessor instead         */

              WDB[mat + MATERIAL_DOPPLER_TEMP] = RDB[ifc + IFC_MIN_TEMP];

              /* Set Doppler preprocessor on */

              WDB[DATA_USE_DOPPLER_PREPROCESSOR] = (double)YES;
            }

          /* Link interface and put flag */

          WDB[mat + MATERIAL_PTR_IFC] = (double)ifc;
          WDB[mat + MATERIAL_USE_IFC] = (double)YES;

          /* Next cell */

          cell = NextItem(cell);
        }
    }
  else
    {
      /* Check that the densities and temperatures stay below the majorants */
      /* and above the minorants */

      loc1 = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];

      while (loc1 > VALID_PTR)
        {
          /* Get geometry cell */

          cell = (long)RDB[loc1 + IFC_TET_PRNT_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Get material from cell or link material to cell */

          mat = (long)RDB[cell + CELL_PTR_MAT];

          /* Get tet-cell index */

          idx = (long)RDB[loc1 + IFC_TET_PRNT_IDX];

          /******************************************************************/
          /* Check that density format corresponds with material definition */
          /* and update material maximum density if needed                  */
          /******************************************************************/

          if (!update)
            {
              /* If no density file is given, we don't have to check this */

              if (strcmp(GetText(ifc + IFC_PTR_OF_RFILE), "-1"))
                {
                  /* Get tet-cell density */

                  d = RDB[dflist + idx];

                  /* Check that tet-cell density is below material density */

                  matdens = RDB[mat + MATERIAL_ADENS];

                  if (fabs(matdens) < fabs(d))
                    {
                      /* IFC density greater than material base density */
                      /* Increase the material density */

                      Warn(FUNCTION_NAME, "Interface density for material %s "
                           "larger than material majorant density (%E > %E)",
                           GetText(mat + MATERIAL_PTR_NAME), d, matdens);
                    }
                }
            }

          /*****************************/
          /* Check material TMS limits */
          /*****************************/

          if (RDB[ifc + IFC_MIN_TEMP] != RDB[ifc + IFC_MAX_TEMP])
            {
              /* Get tet-cell temperature */

              T = RDB[tmplist + idx];

              /* Compare to material maximum and minimum */

              if (T > RDB[mat + MATERIAL_TMS_TMAX])
                {
                  Warn(FUNCTION_NAME,
                       "Material temperature above TMS majorant for material %s (%E > %E)",
                       GetText(mat + MATERIAL_PTR_NAME), T, RDB[mat + MATERIAL_TMS_TMAX]);
                }

              /* Put minimum temperature */

              if (T < RDB[mat + MATERIAL_TMS_TMIN])
                {
                  Warn(FUNCTION_NAME,
                       "Material temperature below TMS minorant for material %s (%E < %E)",
                       GetText(mat + MATERIAL_PTR_NAME), T, RDB[mat + MATERIAL_TMS_TMIN]);
                }
            }

          loc1 = NextItem(loc1);
        }

    }

  /* Create search mesh */

  if (!update)
    {

      /* Get limits */

      xmin = RDB[ifc + IFC_MESH_XMIN];
      xmax = RDB[ifc + IFC_MESH_XMAX];
      ymin = RDB[ifc + IFC_MESH_YMIN];
      ymax = RDB[ifc + IFC_MESH_YMAX];
      zmin = RDB[ifc + IFC_MESH_ZMIN];
      zmax = RDB[ifc + IFC_MESH_ZMAX];

      /* Check boundaries */

      if ((xmin >= xmax) || (ymin >= ymax) || (zmin >= zmax))
        Error(ifc, "Structure is not 3D");

      /* Adjust boundaries */

      xmin = xmin - 1E-6;
      xmax = xmax + 1E-6;
      ymin = ymin - 1E-6;
      ymax = ymax + 1E-6;
      zmin = zmin - 1E-6;
      zmax = zmax + 1E-6;

      /* Put mesh variables */

      lims[0] = xmin;
      lims[1] = xmax;
      lims[2] = ymin;
      lims[3] = ymax;
      lims[4] = zmin;
      lims[5] = zmax;

      /* Read mesh split criterion */

      np = (long)RDB[ifc + IFC_SEARCH_MESH_ADA_SPLIT];
      CheckValue(FUNCTION_NAME, "np", "", np, 1, INFTY);

      /* Get pointer to size vector */

      ptr = (long)RDB[ifc + IFC_SEARCH_MESH_ADA_PTR_SZ];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Create mesh structure */

      msh = CreateMesh(MESH_TYPE_ADAPTIVE, MESH_CONTENT_PTR,
                       MESH_CONTENT_DATA_TET, np, 0, 0, lims, ptr);

      /* Put pointer */

      WDB[ifc + IFC_PTR_SEARCH_MESH_LIST] = (double)msh;

      /* Print progress */

      fprintf(outp, "\nCreating search mesh for interface %s:\n\n",
              GetText(ifc + IFC_PTR_INPUT_FNAME));

      last =  0.0;
      count = 0;

      /* Loop over tet cells to create search mesh  */
      /* This obviously has to be done for children */

      tetlist = (long)RDB[ifc + IFC_PTR_TET_LIST];

      ntet = (long)RDB[ifc + IFC_NC];

      for (i = 0; i < ntet; i++)
        {
          /* Preallocate memory for counters */

          if ((long)RDB[DATA_REAL_PRIVA_SIZE] -
              (long)RDB[DATA_ALLOC_PRIVA_SIZE] < 100)
            PreallocMem(100, PRIVA_ARRAY);

          /* Calculate fraction and print progress */

          frac = (double)(count++)/((double)ntet);

          if (frac - last > 0.10)
            {
              fprintf(outp, " %3.0f%% complete\n", 100.0*frac);
              last = frac;
            }

          /* Get pointer to tet */

          tet = (long)RDB[tetlist + i];

          /* Calculate limits based on points */

          CalculateTetBoundingBox(tet, bb);

          /* Add to search mesh */

          AddSearchMesh(msh, tet, bb[0], bb[1], bb[2], bb[3], bb[4], bb[5]);

        }
    }

  fprintf(outp, " %3.0f%% complete\n\n", 100.0);

  fprintf(outp, "Processing densities for interface %s...\n",
          GetText(ifc + IFC_PTR_INPUT_FNAME));

  /* Loop over tet cells to set density factors */
  /* We'll only have to handle parents since the densities of the children */
  /* are read from the same list */

  loc1 = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];

  /* Reset minimum and maximum DF (used for printing) */

  dfmin = 1.0;
  dfmax = 0.0;

  /* If density file was not given, no need to update material */
  /* densities or majorant densities or anything */

  if (!strcmp(GetText(ifc + IFC_PTR_OF_RFILE), "-1"))
    {
      /* No density file given, density factor is always 1.0 */

      while (loc1 > VALID_PTR)
        {
          /* Get tet-cell index */

          idx = (long)RDB[loc1 + IFC_TET_PRNT_IDX];

          /* Simply set the density factor to 1.0 */

          WDB[dflist + idx] = 1.0;

          /* Get next tet-cell */

          loc1 = NextItem(loc1);
        }

      /* Set minimum and maximum density factors */

      dfmin = 1.0;
      dfmax = 1.0;
    }
  else
    {

      /* Density file was given, need to calculate density factors */

      while (loc1 > VALID_PTR)
        {
          /* Get tet-cell index */

          idx = (long)RDB[loc1 + IFC_TET_PRNT_IDX];

          /* Get tet-cell density */

          d = RDB[dflist + idx];

          /* Get geometry cell */

          cell = (long)RDB[loc1 + IFC_TET_PRNT_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Get material */

          mat = (long)RDB[cell + CELL_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Get majorant density to use for calculating the DF */

          if (update)
            {
              /* Check whether atomic or mass densities are given in interface */

              if (RDB[ifc + IFC_MAX_DENSITY] > 0)
                {
                  /* Atomic densities given */

                  matdens = RDB[mat + MATERIAL_ADENS];
                }
              else
                {
                  /* Mass densities given */

                  matdens = -RDB[mat + MATERIAL_MDENS];
                }
            }
          else
            {
              /* If this is the first time this interface is processed we will use */
              /* MATERIAL_ADENS for the density regardless of IFC units */
              /* They will be converted to atomic density later         */

              matdens = RDB[mat + MATERIAL_ADENS];
            }

          /* Calculate density factor */

          df = d/matdens;

          /* Store density factor */

          WDB[dflist + idx] = df;

          /* Compare to minimum and maximum */

          if (df > dfmax)
            dfmax = df;
          else if (df < dfmin)
            dfmin = df;

          /* Check that the density factor is below unity*/

          if (RDB[dflist + idx] > 1.0)
            Die(FUNCTION_NAME, "Larger than unity density factor for interface %s: %E",
                GetText(ifc + IFC_PTR_INPUT_FNAME), RDB[dflist + idx]);

          if (WDB[dflist + idx] < 0)
            Die(FUNCTION_NAME, "Negative density factor for interface %s: %E",
                GetText(ifc + IFC_PTR_INPUT_FNAME), WDB[dflist + idx]);

          /* Next */

          loc1 = NextItem(loc1);

        }
    }

  /* In time dependent coupled calculation, we'll create separate lists for BOI */
  /* To allow interpolation of temperature/density data inside the interval     */

  if ((RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN) && (update == 0))
    {
      /* Get number of parent cells */

      nc = (long)RDB[ifc + IFC_NC_PARENTS];

      /* Allocate memory for temperatures */

      tmplist0 = ReallocMem(DATA_ARRAY, nc);

      /* Put temperature list pointer to memory */

      WDB[ifc + IFC_PTR_TMP_LIST_BOI] = (double)tmplist0;

      /* Copy temperature data from EOI to BOI */

      memcpy(&WDB[tmplist0], &RDB[tmplist], nc*sizeof(double));

      /* Allocate memory for densities */

      dflist0 = ReallocMem(DATA_ARRAY, nc);

      /* Put density list pointer to memory */

      WDB[ifc + IFC_PTR_DF_LIST_BOI] = (double)dflist0;

      /* Copy density data from EOI to BOI */

      memcpy(&WDB[dflist0], &RDB[dflist], nc*sizeof(double));
    }

  /***************************************************************************/

  /****************************************************************************/
  /***** Print out basic things from the interface ****************************/
  /****************************************************************************/

  /* Get limits */

  if (!update)
    {
      xmin = RDB[ifc + IFC_MESH_XMIN];
      xmax = RDB[ifc + IFC_MESH_XMAX];
      ymin = RDB[ifc + IFC_MESH_YMIN];
      ymax = RDB[ifc + IFC_MESH_YMAX];
      zmin = RDB[ifc + IFC_MESH_ZMIN];
      zmax = RDB[ifc + IFC_MESH_ZMAX];

      fprintf(outp, "\nOpenFOAM mesh based interface \"%s\":\n\n",
              GetText(ifc + IFC_PTR_INPUT_FNAME));

      fprintf(outp, " - Dimensions: x = [%1.1f, %1.1f]; ",
              xmin, xmax);
      fprintf(outp, "y = [%1.1f, %1.1f]; ",
              ymin, ymax);
      fprintf(outp, "z = [%1.1f, %1.1f]\n",
              zmin, zmax);

      /* Get number of cells */

      nc = (long)RDB[ifc + IFC_NC_PARENTS];

      /* Get number of child cells */

      ncc = (long)RDB[ifc + IFC_NC];

      /* Print number of cells */

      if (nc != 0)
        {
          fprintf(outp, " - Initial number of cells: %ld\n", nc);

          /* Print number of child cells */

          fprintf(outp, " - Divided into a number of child cells: %ld\n", ncc);

        }
      else
        fprintf(outp, " - Number of cells: %ld\n", ncc);

      /* Get maximum density */

      if ((dmax = RDB[ifc + IFC_MAX_DENSITY]) == 0.0)
        if (strcmp(GetText(ifc + IFC_PTR_OF_RFILE), "-1"))
          Die(FUNCTION_NAME, "Zero maximum density");

      /* Print IFC-density limits */

      if (strcmp(GetText(ifc + IFC_PTR_OF_RFILE), "-1"))
        {
          /* Density file given */

          fprintf(outp, " - Maximum density of interface: %E\n",
                  RDB[ifc + IFC_MAX_DENSITY]);
        }
      else
        {

          /* Density file not given */

          if (RDB[ifc + IFC_PTR_OF_MFILE] < VALID_PTR)
            {
              fprintf(outp, " - Using mat-card material density:\n");
            }
          else
            fprintf(outp, " - Using mat-card material densities:\n");

          /* Print material densities based on geometry cells */

          /* Get pointer to linked geometry cells */

          cell = (long)RDB[ifc + IFC_PTR_GCELL_LIST];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Loop over cell list */

          while (cell > VALID_PTR)
            {

              /* Get pointer to cell material */

              mat = (long)RDB[cell + CELL_PTR_MAT];

              /* Print density */

              if ((matdens = RDB[mat + MATERIAL_ADENS]) < 0)
                {
                  /* Mass densities given*/

                  fprintf(outp, "%12s --- %.4E g/cm3\n",
                          GetText(mat + MATERIAL_PTR_NAME), -matdens);
                }
              else
                {
                  /* Print atomic density */

                  fprintf(outp, "%12s --- %.4E *10^24 atoms/cm3\n",
                          GetText(mat + MATERIAL_PTR_NAME), matdens);
                }

              /* Next geometry cell */

              cell = NextItem(cell);
            }
        }

      /* Print minimum and maximum density factors */

      fprintf(outp, " - Density factors between %.4f and %.4f.\n", dfmin, dfmax);

      /* Print maximum density and TMS limits of linked material(s) */

      if (RDB[ifc + IFC_PTR_OF_MFILE] < VALID_PTR)
        {

          /* Get material array pointer */

          mat = (long)RDB[ifc + IFC_PTR_MAT_ARR];
          CheckPointer(FUNCTION_NAME, "(mat_arr)", DATA_ARRAY, mat);

          /* Get pointer to the only material */

          mat = (long)RDB[mat];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Print material majorant density */

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
      else
        {
          fprintf(outp, " - Possibly multiple materials in interface\n");
        }
    }

  /***************************************************************************/

}

/*****************************************************************************/
