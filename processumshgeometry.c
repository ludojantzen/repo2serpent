/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processumshgeometry.c                          */
/*                                                                           */
/* Created:       2013/11/23 (JLe)                                           */
/* Last modified: 2019/10/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Processes unstructured mesh based geometry                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessUMSHGeometry:"

/*****************************************************************************/

void ProcessUMSHGeometry()
{

  long umsh, loc0, loc1, ptr, count, msh, np;
  long mat0, mat, cell, prnt, tetlist, tet, i, ntet;
  char *mname;
  double xmin, xmax, ymin, ymax, zmin, zmax, lims[6], frac, last, bb[6];

  /* Check pointer */

  if ((umsh = (long)RDB[DATA_PTR_UMSH0]) < VALID_PTR)
    return;

  fprintf(outp, "Processing unstructured mesh based geometries...\n");

  /* Loop over definitions */

  while (umsh > VALID_PTR)
    {

      /***********************************************************************/

      /***** Create search mesh **********************************************/

      /* Pointer to interface structure */

      loc0 = (long)RDB[umsh + UMSH_PTR_IFC];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Loop over tet cells to link materials */
      /* Only do this to parents */

      loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH_PARENTS];
      CheckPointer(FUNCTION_NAME, "(loc1_)", DATA_ARRAY, loc1);

      /* Reset previous material */

      mat0 = -1;

      while (loc1 > VALID_PTR)
        {
          /* Get pointer to parent (if exists) or self */

          prnt = loc1;
          CheckPointer(FUNCTION_NAME, "(prnt)", DATA_ARRAY, prnt);

          /* Get geometry cell */

          cell = (long)RDB[prnt + IFC_TET_PRNT_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          if ((RDB[cell + CELL_PTR_MAT]) < 0)
            {
              /* Not yet linked */

              /* Remove minus sign */

              WDB[cell + CELL_PTR_MAT] = -RDB[cell + CELL_PTR_MAT];

              /* Get material name */

              mname = GetText(cell + CELL_PTR_MAT);

              /* Check previous */

              mat = -1;

              if (mat0 > VALID_PTR)
                if (!strcmp(mname, GetText(mat0 + MATERIAL_PTR_NAME)))
                  mat = mat0;

              /* Find match */

              if (mat < VALID_PTR)
                {
                  mat = (long)RDB[DATA_PTR_M0];
                  while (mat > VALID_PTR)
                    {
                      /* Compare name */

                      if (!strcmp(mname, GetText(mat + MATERIAL_PTR_NAME)))
                        break;

                      /* Check if material was divided for burnup calculation */

                      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT])
                          > VALID_PTR)
                        if (!strcmp(mname, GetText(ptr + MATERIAL_PTR_NAME)))
                          break;

                      /* Next material */

                      mat = NextItem(mat);
                    }
                }

              /* Check match */

              if (mat < VALID_PTR)
                Error(loc0, "Material %s is not defined", mname);

              /* Set pointer */

              WDB[cell + CELL_PTR_MAT] = (double)mat;
              mat0 = mat;

              /* Put interface to material, if this is an interface solid */

              if(RDB[loc0 + IFC_PTR_OF_TFILE] > VALID_PTR)
                {
                  WDB[mat + MATERIAL_PTR_IFC] = (double)loc0;
                }
            }
          else
            {
              /* Already linked (OF_SOLID) */

              mat = (long)RDB[cell + CELL_PTR_MAT];

            }

          loc1 = NextItem(loc1);

        }

      /* Copy number of original cells */

      if ((long)RDB[loc0 + IFC_NC_PARENTS] > 0)
        WDB[umsh + UMSH_N_ORIG_CELLS] = RDB[loc0 + IFC_NC_PARENTS];
      else
        WDB[umsh + UMSH_N_ORIG_CELLS] = RDB[loc0 + IFC_NC];

      /* Check if interface based umsh */

      if((long)RDB[umsh + UMSH_PTR_FNAME] > 0)
        {
          /* Other things for these will be processed in processtetmesh.c */

          umsh = NextItem(umsh);

          continue;
        }


      /* Get limits */

      xmin = RDB[loc0 + IFC_MESH_XMIN];
      xmax = RDB[loc0 + IFC_MESH_XMAX];
      ymin = RDB[loc0 + IFC_MESH_YMIN];
      ymax = RDB[loc0 + IFC_MESH_YMAX];
      zmin = RDB[loc0 + IFC_MESH_ZMIN];
      zmax = RDB[loc0 + IFC_MESH_ZMAX];

      /* Check boundaries */

      if ((xmin >= xmax) || (ymin >= ymax) || (zmin >= zmax))
        Error(umsh, "Structure is not 3D");

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

      np = (long)RDB[loc0 + IFC_SEARCH_MESH_ADA_SPLIT];
      CheckValue(FUNCTION_NAME, "np", "", np, 1, INFTY);

      /* Get pointer to size vector */

      ptr = (long)RDB[loc0 + IFC_SEARCH_MESH_ADA_PTR_SZ];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Create mesh structure */

      msh = CreateMesh(MESH_TYPE_ADAPTIVE, MESH_CONTENT_PTR,
                       MESH_CONTENT_DATA_TET, np, 0, 0, lims, ptr);

      /* Put pointer */

      WDB[loc0 + IFC_PTR_SEARCH_MESH_LIST] = (double)msh;

      /* Print progress */

      fprintf(outp, "\nCreating search mesh for universe %s:\n\n",
              GetText(umsh + UMSH_PTR_NAME));

      last =  0.0;
      count = 0;

      /* Loop over tet cells to create search mesh  */
      /* This obviously has to be done for children */

      tetlist = (long)RDB[loc0 + IFC_PTR_TET_LIST];

      ntet = (long)RDB[loc0 + IFC_NC];

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

      fprintf(outp, " %3.0f%% complete\n\n", 100.0);

      /***********************************************************************/


      /***********************************************************************/

      /* Next geometry */

      umsh = NextItem(umsh);
    }
}

/*****************************************************************************/
