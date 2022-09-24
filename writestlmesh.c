/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writestlmesh.c                                 */
/*                                                                           */
/* Created:       2019/08/30 (JLe)                                           */
/* Last modified: 2019/10/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Writes STL search mesh into a binary file.                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WriteSTLMesh:"

/* Recursive write function */

void WriteSTLMesh0(FILE *fp, long msh, long mode);

/*****************************************************************************/

void WriteSTLMesh(long stl)
{
  long msh, ptr, loc0, loc1, i, j, n, pos0;
  FILE *fp;
  char outfile[MAX_STR], *stp;

  /* Check option */

  if ((long)RDB[DATA_STL_MESH_FILES] == NO)
    return;

  /* Print note */

  fprintf(outp, "Writing STL search mesh in file...\n");

  /* Output file name */

  sprintf(outfile, "%s.smh", GetText(DATA_PTR_INPUT_FNAME));

  /* Open data file */

  if (PrevItem(stl) < VALID_PTR)
    fp = fopen(outfile, "w");
  else
    fp = fopen(outfile, "a");

  /* Get initial position */

  pos0 = ftell(fp);

  /* Pointer to universe */

  ptr = (long)RDB[stl + STL_PTR_UNI];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get universe name */

  stp = GetText(ptr + UNIVERSE_PTR_NAME);

  /* Write universe name */

  n = strlen(stp);
  fwrite(&n, sizeof(long), 1, fp);
  fwrite(stp, sizeof(char), n, fp);

  /* Write split limit */

  n = (long)RDB[stl + STL_SEARCH_MESH_ADA_SPLIT];
  fwrite(&n, sizeof(long), 1, fp);

  /* Loop over split criteria */

  ptr = (long)RDB[stl + STL_SEARCH_MESH_PTR_SZ];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Loop over items */

  while (1 != 2)
    {
      /* Write size */

      n = (long)RDB[ptr++];
      fwrite(&n, sizeof(long), 1, fp);

      /* Check terminator */

      if (n < 0)
        break;
    }

  /* Write merge radius */

  fwrite(&RDB[stl + STL_MERGE_RAD], sizeof(double), 1, fp);

  /* Reset solid count (zero is for bacground universe) */

  i = 1;

  /* Loop over stl solids */

  loc0 = (long)RDB[stl + STL_PTR_SOLIDS];
  while (loc0 > VALID_PTR)
    {
      /* Put index */

      if (i > INT_MAX - 1)
        Error(stl, "Maximum number of solids per universe exceeded");
      else
        WDB[loc0 + STL_SOLID_IDX] = (double)(i++);

      /* Check number of facets */

      if ((long)RDB[loc0 + STL_SOLID_N_FACETS] > INT_MAX - 1)
        Error(stl, "Maximum number of facets per solid exceeded");

      /* Pointer to facets */

      loc1 = (long)RDB[loc0 + STL_SOLID_PTR_FACETS];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop over facets */

      for (j = 0; j < (long)RDB[loc0 + STL_SOLID_N_FACETS]; j++)
        {
          /* Put index */

          WDB[loc1 + STL_FACET_IDX] = (double)j;

          /* Pointer to next */

          loc1 = loc1 + STL_FACET_BLOCK_SIZE;
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
        }

      /* Next solid */

      loc0 = NextItem(loc0);
    }

  /* Write solid mesh */

  msh = (long)RDB[stl + STL_PTR_SOLID_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
  WriteSTLMesh0(fp, msh, 1);

  /* Write facet mesh */

  msh = (long)RDB[stl + STL_PTR_FACET_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
  WriteSTLMesh0(fp, msh, 2);

  /* Write initial position */

  fwrite(&pos0, sizeof(long), 1, fp);

  /* Close file */

  fclose(fp);

  /* Exit OK */

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/

/***** Recursive write function **********************************************/

void WriteSTLMesh0(FILE *fp, long msh, long mode)
{
  long nx, ny, nz, i, j, k, ptr, lst, sld, fct;
  int idx;

  /* Check mesh pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get size */

  nx = (long)RDB[msh + MESH_N0];
  CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 10000000);

  ny = (long)RDB[msh + MESH_N1];
  CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 10000000);

  nz = (long)RDB[msh + MESH_N2];
  CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 10000000);

  /* Write size and boundaries */

  fwrite(&nx, sizeof(long), 1, fp);
  fwrite(&ny, sizeof(long), 1, fp);
  fwrite(&nz, sizeof(long), 1, fp);

  fwrite(&RDB[msh + MESH_MIN0], sizeof(double), 1, fp);
  fwrite(&RDB[msh + MESH_MAX0], sizeof(double), 1, fp);
  fwrite(&RDB[msh + MESH_MIN1], sizeof(double), 1, fp);
  fwrite(&RDB[msh + MESH_MAX1], sizeof(double), 1, fp);
  fwrite(&RDB[msh + MESH_MIN2], sizeof(double), 1, fp);
  fwrite(&RDB[msh + MESH_MAX2], sizeof(double), 1, fp);

  /* Loop over mesh */

  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
        {
          /* Pointer to data */

          ptr = ReadMeshPtr(msh, i, j, k);
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Check pointer */

          if ((lst = (long)RDB[ptr]) > VALID_PTR)
            {
              /* Write 1 to indicate that data follows */

              idx = 1;
              fwrite(&idx, sizeof(int), 1, fp);

              /* Check if preassigned or list of items */

              if ((sld = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT]) < VALID_PTR)
                {
                  /* Pre-assigned data */

                  if (-sld < VALID_PTR)
                    {
                      /* Background universe, write zero */

                      idx = 0;
                      fwrite(&idx, sizeof(int), 1, fp);
                    }
                  else
                    {
                      /* Write negative solid index */

                      idx = -(int)RDB[-sld + STL_SOLID_IDX];
                      fwrite(&idx, sizeof(int), 1, fp);
                    }
                }
              else
                {
                  /* List of items */

                  while (lst > VALID_PTR)
                    {
                      /* Check mode */

                      if (mode == 1)
                        {
                          /* Solid mesh */

                          sld = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
                          CheckPointer(FUNCTION_NAME, "(sld)",
                                       DATA_ARRAY, sld);

                          /* Write solid index */

                          idx = (int)RDB[sld + STL_SOLID_IDX];
                          fwrite(&idx, sizeof(int), 1, fp);
                        }
                      else
                        {
                          /* Facet mesh */

                          fct = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
                          CheckPointer(FUNCTION_NAME, "(fct)",
                                       DATA_ARRAY, fct);

                          /* Pointer to solid */

                          sld = (long)RDB[fct + STL_FACET_PTR_SOLID];
                          CheckPointer(FUNCTION_NAME, "(sld)",
                                       DATA_ARRAY, sld);

                          /* Write solid index */

                          idx = (int)RDB[sld + STL_SOLID_IDX];
                          fwrite(&idx, sizeof(int), 1, fp);

                          /* Write facet index */

                          idx = (int)RDB[fct + STL_FACET_IDX];
                          fwrite(&idx, sizeof(int), 1, fp);
                        }

                      /* Next */

                      lst = NextItem(lst);
                    }

                  /* Write -1 to terminate */

                  idx = -1;
                  fwrite(&idx, sizeof(int), 1, fp);
                }
            }
          else if (lst == 0)
            {
              /* Solid mesh may have empty cells */

              if (mode == 1)
                {
                  /* Write zero */

                  idx = 0;
                  fwrite(&idx, sizeof(int), 1, fp);
                }
              else
                Die(FUNCTION_NAME, "Zero in facet mesh");

            }
          else if (lst != NULLPTR)
            {
              /* Write -1 to indicate that another mesh follows */

              idx = -1;
              fwrite(&idx, sizeof(int), 1, fp);

              /* Call recursively */

              WriteSTLMesh0(fp, -lst, mode);
            }
          else
            Die(FUNCTION_NAME, "WTF?");
        }
}

/*****************************************************************************/
