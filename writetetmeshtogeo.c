/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writetetmeshtogeo.c                            */
/*                                                                           */
/* Created:       2013/12/27 (JLe)                                           */
/* Last modified: 2018/01/26 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Used to write tet mesh format multi-physics interface to     */
/*              normal surface-cell type geometry. Used for testing and      */
/*              debugging purposes only.                                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WriteTetMeshtoGeo:"

/*****************************************************************************/

void WriteTetMeshtoGeo()
{
#ifdef blohblah
  long loc0, loc1, loc2, loc3, cell, surf, ptr, ptr1, n, mat, nuc;
  long nd, np, side, i, j, idx, prnt;
  long sslist, tmplist, dflist;
  FILE *fp;

  if ((long)RDB[DATA_PTR_IFC0] < VALID_PTR)
    return;

  fprintf(outp, "Writing tet mesh interface into geometry file \"geo.inp\"...\n");

  /* Open file for writing */

  if ((fp = fopen("geo.inp", "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Loop over interafaces */

  loc0 = (long)RDB[DATA_PTR_IFC0];

  while (loc0 > VALID_PTR)
    {
      /* Check type */

      if ((long)RDB[loc0 + IFC_TYPE] != IFC_TYPE_TET_MESH)
        {
          /* Next interface */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /* Get pointer to surface list */

      sslist = (long)RDB[loc0 + IFC_PTR_SURF_LIST];
      CheckPointer(FUNCTION_NAME, "(sslist)", DATA_ARRAY, sslist);

      /* Get pointer to temperature list */

      tmplist = (long)RDB[loc0 + IFC_PTR_TMP_LIST];

      /* Get pointer to density list */

      dflist = (long)RDB[loc0 + IFC_PTR_DF_LIST];

      /***********************************************************************/

      /***** Print surface data **********************************************/

      /* Loop over tets */

      loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];

      while (loc1 > VALID_PTR)
        {
          /* Get pointer to parent cell (if exists) or self */

          if ((prnt = (long)RDB[loc1 + IFC_TET_MSH_PTR_PARENT]) < VALID_PTR)
            prnt = loc1;

          /* Get number of faces */

          nd = (long)RDB[loc1 + IFC_TET_MSH_NF];

          /* Get pointer to geometry cell */

          cell = (long)RDB[prnt + IFC_TET_PRNT_PTR_CELL];

          /* Get pointer to material */
          /* NO: Material names are replaced with material pointers only at processumshgeom... */
          /*  mat0 = (long)RDB[cell + CELL_PTR_MAT];*/

          mat = (long)RDB[cell + CELL_PTR_MAT];
          printf("mat = %ld\n", mat);

          printf("****\nCell %ld, material %s, %ld faces\n",
                 loc1,
                 GetText(mat + MATERIAL_PTR_NAME), nd);

          ptr = (long)RDB[loc1 + IFC_TET_MSH_PTR_FACES];
          ptr1 = (long)RDB[loc1 + IFC_TET_MSH_PTR_SIDES];

          /* Loop over faces */

          for (i = 0; i < nd ; i++)
            {

              printf("Face %ld, surf %ld, side %ld\n", i, (long)RDB[ptr + i], (long)RDB[ptr1 + i]);

              /* Get surface index */

              n = (long)RDB[ptr + i];
              surf = ListPtr(sslist, n);

              side = (long)RDB[ptr1 + i];

              np = (long)RDB[surf + SURFACE_N_PARAMS];

              /* Print surface definition if owner of this face */

              if (side == -1)
                {
                  fprintf(fp, "surf i%lds%ld plane ", (long)RDB[loc0 + IFC_IDX], n);

                  /* Get pointer to surface parameters (point list) */

                  loc2 = (long)RDB[surf + SURFACE_PTR_PARAMS];

                  /* Loop over points */

                  for (j = 0; j < np; j++)
                    {
                      /* Get point */
                      loc3 = (long)RDB[loc2 + j];

                      /* Print point */
                      fprintf(fp, "%E %E %E ", RDB[loc3], RDB[loc3 + 1], RDB[loc3 + 2]);
                    }

                  fprintf(fp,"\n");

                }

            }

          loc1 = NextItem(loc1);
        }

      fprintf(fp, "\n");

      /* Create geometry cells */

      loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];

      while (loc1 > VALID_PTR)
        {

          /* Get pointer to parent (if exists) or self */

          if ((prnt = (long)RDB[loc1 + IFC_TET_MSH_PTR_PARENT]) < VALID_PTR)
            prnt = loc1;

          fprintf(fp, "cell i%ldc%ld 1 ", (long)RDB[loc0 + IFC_IDX],
                  loc1);
          fprintf(fp, "i%ldm%ld ", (long)RDB[loc0 + IFC_IDX],
                  loc1);

          /* Get number of faces */

          nd = (long)RDB[loc1 + IFC_TET_MSH_NF];

          /* Get pointer to facelist */

          ptr = (long)RDB[loc1 + IFC_TET_MSH_PTR_FACES];

          /* Get pointer to sidelist */

          ptr1 = (long)RDB[loc1 + IFC_TET_MSH_PTR_SIDES];

          /* Loop over faces */

          for (i = 0; i < nd ; i++)
            {

              /* Get surface index */

              n = (long)RDB[ptr + i];
              surf = ListPtr(sslist, n);

              side = (long)RDB[ptr1 + i];

              np = (long)RDB[surf + SURFACE_N_PARAMS];

              /* Print surface definition if owner of this face */

              if (side  == -1)
                fprintf(fp, "-i%lds%ld ", (long)RDB[loc0 + IFC_IDX], n);
              else
                fprintf(fp, "i%lds%ld ", (long)RDB[loc0 + IFC_IDX], n);

            }

          fprintf(fp, "\n");

          loc1 = NextItem(loc1);
        }

      fprintf(fp, "\n");

      /***** Print material data *********************************************/

      /* Loop over cells */

      loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
      while (loc1 > VALID_PTR)
        {

          /* Get parent cell index */

          if ((prnt = (long)RDB[loc1 + IFC_TET_MSH_PTR_PARENT]) < VALID_PTR)
            prnt = loc1;

          idx = (long)RDB[prnt + IFC_TET_PRNT_IDX];

          /* Get pointer to geometry cell */

          cell = (long)RDB[prnt + IFC_TET_PRNT_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          mat  = (long)RDB[cell + CELL_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Print material data */

          fprintf(fp, "\nmat i%ldm%ld %1.5E tms %1.5E\n",
                  (long)RDB[loc0 + IFC_IDX], loc1,
                  -RDB[mat + MATERIAL_MDENS]*RDB[dflist + idx],
                  RDB[tmplist + idx]);

          /* Loop over composition */

          ptr = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (ptr > VALID_PTR)
            {
              /* Pointer to nuclide */

              nuc = (long)RDB[ptr + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Print */

              fprintf(fp, "%10s %1.5E\n", GetText(nuc + NUCLIDE_PTR_NAME),
                      RDB[ptr + COMPOSITION_ADENS]/RDB[mat + MATERIAL_ADENS]);

              /* Next */

              ptr = NextItem(ptr);
            }

          fprintf(fp, "\n");

          /* Next cell */

          loc1 = NextItem(loc1);
        }

      fprintf(fp, "\n");

      /***** Write options ***************************************************/

      /* Search mesh */

      fprintf(fp, "set csm 1 1 ");

      fprintf(fp, "%E %E ", RDB[loc0 + IFC_MESH_XMIN],
              RDB[loc0 + IFC_MESH_XMAX]);

      if (RDB[loc0 + IFC_SEARCH_MESH_NX] > 10.0)
        fprintf(fp, "%ld ", (long)(RDB[loc0 + IFC_SEARCH_MESH_NX]/10.0));
      else
        fprintf(fp, "1 ");

      fprintf(fp, "%E %E ", RDB[loc0 + IFC_MESH_YMIN],
              RDB[loc0 + IFC_MESH_YMAX]);

      if (RDB[loc0 + IFC_SEARCH_MESH_NY] > 10.0)
        fprintf(fp, "%ld ", (long)(RDB[loc0 + IFC_SEARCH_MESH_NY]/10.0));
      else
        fprintf(fp, "1 ");

      fprintf(fp, "%E %E ", RDB[loc0 + IFC_MESH_ZMIN],
              RDB[loc0 + IFC_MESH_ZMAX]);

      if (RDB[loc0 + IFC_SEARCH_MESH_NZ] > 10.0)
        fprintf(fp, "%ld ", (long)(RDB[loc0 + IFC_SEARCH_MESH_NZ]/10.0));
      else
        fprintf(fp, "1 ");

      fprintf(fp, "\n");

      /* Treat undefined regions as void */

      fprintf(fp, "set voidc 1\n\n");

      /***********************************************************************/

      /* Next interface */

      loc0 = NextItem(loc0);
    }

  /* Close file */

  fclose(fp);

  fprintf(outp, "OK.\n\n");

  /* Terminate run */

  exit(0);
#endif
}

/*****************************************************************************/
