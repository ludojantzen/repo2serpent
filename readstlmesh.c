/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writestlmesh.c                                 */
/*                                                                           */
/* Created:       2019/08/30 (JLe)                                           */
/* Last modified: 2020/05/23 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Writes STL search mesh into a binary file.                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadSTLMesh:"

/* Recursive write function */

long ReadSTLMesh0(FILE *fp, long stl, long mode);

/*****************************************************************************/

long ReadSTLMesh(long stl)
{
  long msh, loc0, loc1, ptr, i, j, sz, n, pos;
  double r;
  FILE *fp;
  char infile[MAX_STR], str[MAX_STR];

  /* Check option */

  if ((long)RDB[DATA_STL_MESH_FILES] == NO)
    return NO;

  /* Output file name */

  sprintf(infile, "%s.smh", GetText(DATA_PTR_INPUT_FNAME));

  /* Open data file */

  if ((fp = fopen(infile, "r")) == NULL)
    return NO;

  /* Loop to last value */

  if ((sz = fseek(fp, -sizeof(double), SEEK_END)) != 0)
    Error(stl, "Error in STL search mesh file");

  /* Loop until correct universe */

  for (i = 0; i < 1000000000; i++)
    {
      /* Read beginning position of data block */

      if ((sz = fread(&pos, sizeof(long), 1, fp)) == 0)
        Error(stl, "Error in STL search mesh file");

      /* Loop to beginning of data block */

      if ((sz = fseek(fp, pos, SEEK_SET)) != 0)
        Error(stl, "Error in STL search mesh file");

      /* Read universe data */

      if ((sz = fread(&n, sizeof(long), 1, fp)) == 0)
        Error(stl, "Error in STL search mesh file");

      if ((sz = fread(str, sizeof(char), n, fp)) == 0)
        Error(stl, "Error in STL search mesh file");

      /* Put null terminator */

      str[n] = '\0';

      /* Pointer to this universe */

      if ((ptr = (long)RDB[stl + STL_PTR_UNI]) < VALID_PTR)
        Error(stl, "Error in universe");

      /* Compare */

      if (!strcmp(str, GetText(ptr + UNIVERSE_PTR_NAME)))
        {
          /* Found match, break loop */

          break;
        }

      /* Position of previous last value */

      if ((pos = pos - sizeof(double)) < 0)
        {
          /* No more data, close file */

          fclose(fp);

          /* Exit */

          return NO;
        }

      /* Loop to data */

      if ((sz = fseek(fp, pos, SEEK_SET)) != 0)
        Error(stl, "Error in STL search mesh file");
    }

  /* Check for infinite loop */

  if (i == 1000000000)
    Die(FUNCTION_NAME, "Infinite loop");

  /* Read split limit */

  if ((sz = fread(&n, sizeof(long), 1, fp)) == 0)
    Error(stl, "Error in STL search mesh file");

  /* Compare */

  if (n != (long)RDB[stl + STL_SEARCH_MESH_ADA_SPLIT])
    Error(stl, "Change in mesh parameters: delete file %s and rerun", infile);

  /* Pointer to search mesh split criteria */

  ptr = (long)RDB[stl + STL_SEARCH_MESH_PTR_SZ];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Loop over items */

  while (1 != 2)
    {
      /* Read size */

      if ((sz = fread(&n, sizeof(long), 1, fp)) == 0)
        Error(stl, "Error in STL search mesh file");

      /* Compare */

      if (n != (long)RDB[ptr++])
        Error(stl, "Change in mesh parameters: delete file %s and rerun",
              infile);

      /* Check terminator */

      if (n < 0)
        break;
    }

  /* Read merge radius */

  if ((sz = fread(&r, sizeof(double), 1, fp)) == 0)
    Error(stl, "Error in STL search mesh file");

  /* Compare */

  if (r != RDB[stl + STL_MERGE_RAD])
    Error(stl, "Change in mesh parameters: delete file %s and rerun", infile);

  /* Reset solid count (zero is for bacground universe) */

  i = 1;

  /* Loop over stl solids */

  loc0 = (long)RDB[stl + STL_PTR_SOLIDS];
  while (loc0 > VALID_PTR)
    {
      /* Put index */

      WDB[loc0 + STL_SOLID_IDX] = (double)(i++);

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

  /* Read solid mesh */

  msh = ReadSTLMesh0(fp, stl, 1);
  WDB[stl + STL_PTR_SOLID_MESH] = (double)msh;

  /* Read facet mesh */

  msh = ReadSTLMesh0(fp, stl, 2);
  WDB[stl + STL_PTR_FACET_MESH] = (double)msh;

  /* Close file */

  fclose(fp);

  /* Exit OK */

  return YES;
}

/*****************************************************************************/

/***** Recursive read function ***********************************************/

long ReadSTLMesh0(FILE *fp, long stl, long mode)
{
  long sz, nx, ny, nz, i, j, k, ptr, loc0, msh, lst, sld, fct;
  int idx;
  double lims[6], V;

  /* read size */

  if ((sz = fread(&nx, sizeof(long), 1, fp)) == 0)
    Error(stl, "Error in STL search mesh file");

  if ((sz = fread(&ny, sizeof(long), 1, fp)) == 0)
    Error(stl, "Error in STL search mesh file");

  if ((sz = fread(&nz, sizeof(long), 1, fp)) == 0)
    Error(stl, "Error in STL search mesh file");

  /* Read boundaries */

  if ((sz = fread(&lims, sizeof(double), 6, fp)) == 0)
    Error(stl, "Error in STL search mesh file");

  /* Check sizes */

  CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 100000);
  CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 100000);
  CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 100000);

  /* Add to size */

  WDB[stl + STL_SEARCH_MESH_CELLS] =
    RDB[stl + STL_SEARCH_MESH_CELLS] + nx*ny*nz;

  /* Unit cell volume */

  V = (lims[1] - lims[0])*(lims[3] - lims[2])*(lims[5] - lims[4])/(nx*ny*nz);
  CheckValue(FUNCTION_NAME, "V", "", V, ZERO, INFTY);

  /* Create mesh */

  msh = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_PTR, -1, nx, ny, nz,
                   lims, -1);

  /* Loop over mesh */

  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
        {
          /* Pointer to data */

          loc0 = ReadMeshPtr(msh, i, j, k);
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          /* Read index */

          if ((sz = fread(&idx, sizeof(int), 1, fp)) == 0)
            Error(stl, "Error in STL search mesh data file");

          /* Check index */

          if (idx == 1)
            {
              /* Data follows */

              if ((sz = fread(&idx, sizeof(int), 1, fp)) == 0)
                Error(stl, "Error in STL search mesh data file");

              /* Check */

              if (idx < 0)
                {
                  /* Pre-assigned solid index, create item */

                  lst = NewItem(loc0, SEARCH_MESH_CELL_BLOCK_SIZE);

                  /* Find solid */

                  sld = (long)RDB[stl + STL_PTR_SOLIDS];
                  while (sld > VALID_PTR)
                    {
                      /* Compare */

                      if ((int)RDB[sld + STL_SOLID_IDX] == -idx)
                        break;

                      /* Pointer to next */

                      sld = NextItem(sld);
                    }

                  /* Check */

                  if (sld < VALID_PTR)
                    Error(stl, "Error in pre-generated file");
                  else
                    WDB[lst + SEARCH_MESH_CELL_CONTENT] = -(double)sld;

                  /* Add to volume */

                  if (mode == 2)
                    WDB[stl + STL_SEARCH_MESH_V] =
                      RDB[stl + STL_SEARCH_MESH_V] + V;
                }
              else if (idx == 0)
                {
                  /* Background universe, create item */

                  lst = NewItem(loc0, SEARCH_MESH_CELL_BLOCK_SIZE);

                  /* Put zero */

                  WDB[lst + SEARCH_MESH_CELL_CONTENT] = 0.0;

                  /* Add to volume */

                  if (mode == 2)
                    WDB[stl + STL_SEARCH_MESH_V] =
                      RDB[stl + STL_SEARCH_MESH_V] + V;
                }
              else if (mode == 1)
                {
                  /* List of solids */

                  do
                    {
                      /* Create new item */

                      lst = NewItem(loc0, SEARCH_MESH_CELL_BLOCK_SIZE);

                      /* Find solid */

                      sld = (long)RDB[stl + STL_PTR_SOLIDS];
                      while (sld > VALID_PTR)
                        {
                          /* Compare */

                          if ((int)RDB[sld + STL_SOLID_IDX] == idx)
                            break;

                          /* Pointer to next */

                          sld = NextItem(sld);
                        }

                      /* Check */

                      if (sld < VALID_PTR)
                        Error(stl, "Error in pre-generated file");
                      else
                        WDB[lst + SEARCH_MESH_CELL_CONTENT] = (double)sld;

                      /* Allocate memory for counter */

                      ptr = AllocPrivateData(1, PRIVA_ARRAY);
                      WDB[lst + SEARCH_MESH_PTR_CELL_COUNT] = (double)ptr;

                      /* Read next index */

                      if ((sz = fread(&idx, sizeof(int), 1, fp)) == 0)
                        Error(stl, "Error in STL search mesh data file");
                    }
                  while (idx > 0);
                }
              else
                {
                  /* List of facets solids */

                  do
                    {
                      /* Create new item */

                      lst = NewItem(loc0, SEARCH_MESH_CELL_BLOCK_SIZE);

                      /* Find solid */

                      sld = (long)RDB[stl + STL_PTR_SOLIDS];
                      while (sld > VALID_PTR)
                        {
                          /* Compare */

                          if ((int)RDB[sld + STL_SOLID_IDX] == idx)
                            break;

                          /* Pointer to next */

                          sld = NextItem(sld);
                        }

                      /* Check */

                      if (sld < VALID_PTR)
                        Error(stl, "Error in pre-generated file");

                      /* Read facet index */

                      if ((sz = fread(&idx, sizeof(int), 1, fp)) == 0)
                        Error(stl, "Error in STL search mesh data file");

                      /* Check */

                      if ((idx < 0) ||
                          (idx > (int)RDB[sld + STL_SOLID_N_FACETS] - 1))
                        Error(stl, "Error in STL search mesh data file");

                      /* Pointer to facets */

                      fct = (long)RDB[sld + STL_SOLID_PTR_FACETS];
                      CheckPointer(FUNCTION_NAME, "(fct)", DATA_ARRAY, fct);
                      fct = fct + idx*STL_FACET_BLOCK_SIZE;

                      /* Put facet pointer */

                      WDB[lst + SEARCH_MESH_CELL_CONTENT] = (double)fct;

                      /* Read next index */

                      if ((sz = fread(&idx, sizeof(int), 1, fp)) == 0)
                        Error(stl, "Error in STL search mesh data file");
                    }
                  while (idx > -1);
                }
            }
          else if (idx == -1)
            {
              /* Another mesh follows, call recursively */

              ptr = ReadSTLMesh0(fp, stl, mode);
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Put pointer */

              WDB[loc0] = -(double)ptr;
            }
          else if (idx == 0)
            {
              /* Solid mesh may have empty cells */

              if (mode != 1)
                Die(FUNCTION_NAME, "Zero in facet mesh");
            }
          else
            Die(FUNCTION_NAME, "WTF?");
        }

  /* Change type to adaptive */

  WDB[msh + MESH_TYPE] = (double)MESH_TYPE_ADAPTIVE;

  /* Set content type */

  if (mode == 1)
    WDB[msh + MESH_CONTENT_DATA_TYPE] = -1.0;
  else
    WDB[msh + MESH_CONTENT_DATA_TYPE] = (double)MESH_CONTENT_DATA_STL;

  /* Return mesh pointer */

  return msh;
}

/*****************************************************************************/
