/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setstlmeshpointers.c                           */
/*                                                                           */
/* Created:       2014/12/10 (JLe)                                           */
/* Last modified: 2019/09/01 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Reconfigures facet and solid pointer lists                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetSTLMeshPointers:"

/*****************************************************************************/

void SetSTLMeshPointers(long msh)
{
  long nx, ny, nz, i, j, k, loc0, loc1, loc2, ptr, sld;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Check type */

  if ((long)RDB[msh + MESH_TYPE] != MESH_TYPE_ADAPTIVE)
    Die(FUNCTION_NAME, "Invalid mesh type");

  /* Check content */

  if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_PTR)
    Die(FUNCTION_NAME, "Invalid content type");

  /* Check data type */

  if ((long)RDB[msh + MESH_CONTENT_DATA_TYPE] != MESH_CONTENT_DATA_STL)
    Die(FUNCTION_NAME, "Not an STL mesh");

  /* Get size */

  nx = (long)RDB[msh + MESH_N0];
  ny = (long)RDB[msh + MESH_N1];
  nz = (long)RDB[msh + MESH_N2];

   /* Loop over content */

  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
        {
          /* Get pointer to content  */

          loc0 = ReadMeshPtr(msh, i, j, k);
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          /* Get pointer to content */

          loc1 = (long)RDB[loc0];

          /* Check */

          if (-loc1 > VALID_PTR)
            {
              /* Pointer to a new mesh, call recursively */

              SetSTLMeshPointers(-loc1);
            }
          else if (loc1 > VALID_PTR)
            {
              /* Copy pointer */

              loc2 = loc1;

              /* Loop over contetnt */

              while (loc1 > VALID_PTR)
                {
                  /* Pointer to facet (negative if filled mesh is read  */
                  /* from file) */

                  if ((ptr = (long)RDB[loc1 + SEARCH_MESH_CELL_CONTENT])
                      < VALID_PTR)
                    break;

                  /* Pointer to solid */

                  sld = (long)RDB[ptr + STL_FACET_PTR_SOLID];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Loop over content and check if exists */

                  ptr = loc2;
                  while (ptr > VALID_PTR)
                    {
                      if (sld == (long)RDB[ptr + SEARCH_MESH_PTR_CELL_COUNT])
                        break;
                      else if ((long)RDB[ptr + SEARCH_MESH_PTR_CELL_COUNT]
                               < VALID_PTR)
                        {
                          /* Put pointer */

                          WDB[ptr + SEARCH_MESH_PTR_CELL_COUNT] = (double)sld;

                          /* Break loop */

                          break;
                        }

                      /* Next */

                      ptr = NextItem(ptr);
                    }

                  /* Next */

                  loc1 = NextItem(loc1);
                }
            }
        }
}

/*****************************************************************************/
