/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fillstlmesh.c                                  */
/*                                                                           */
/* Created:       2014/03/03 (JLe)                                           */
/* Last modified: 2019/08/31 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Fills STL mesh with cell pointers                            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FillSTLMesh:"

/*****************************************************************************/

void FillSTLMesh(long stl, long msh, double x0, double y0, double z0)
{
  long nx, ny, nz, i, j, k, new, ptr, sld;
  double xmin, xmax, ymin, ymax, zmin, zmax, px, py, pz, x, y, z, u, v, w, l;
  unsigned long seed;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
  CheckPointer(FUNCTION_NAME, "(stl)", DATA_ARRAY, stl);

  /* Check type */

  if ((long)RDB[msh + MESH_TYPE] != MESH_TYPE_ADAPTIVE)
    Die(FUNCTION_NAME, "Invalid mesh type");

  /* Check content */

  if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_PTR)
    Die(FUNCTION_NAME, "Invalid content type");

  /* Check data type */

  if ((long)RDB[msh + MESH_CONTENT_DATA_TYPE] != MESH_CONTENT_DATA_STL)
    Die(FUNCTION_NAME, "Not an STL mesh");

  /* Init random number sequence */

  seed = ReInitRNG(0);
  SEED[0] = seed;

  /* Get size */

  nx = (long)RDB[msh + MESH_N0];
  ny = (long)RDB[msh + MESH_N1];
  nz = (long)RDB[msh + MESH_N2];

  /* Get mesh boundaries */

  xmin = RDB[msh + MESH_MIN0];
  xmax = RDB[msh + MESH_MAX0];
  ymin = RDB[msh + MESH_MIN1];
  ymax = RDB[msh + MESH_MAX1];
  zmin = RDB[msh + MESH_MIN2];
  zmax = RDB[msh + MESH_MAX2];

  /* Calculate mesh pitch */

  px = (xmax - xmin)/((double)nx);
  py = (ymax - ymin)/((double)ny);
  pz = (zmax - zmin)/((double)nz);

  /* Add to size */

  WDB[stl + STL_SEARCH_MESH_CELLS] =
    RDB[stl + STL_SEARCH_MESH_CELLS] + (double)(nx*ny*nz);

   /* Loop over content */

  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
        {
          /* Get pointer to content  */

          ptr = ReadMeshPtr(msh, i, j, k);

          /* Check */

          if ((new = -(long)RDB[ptr]) > VALID_PTR)
            {
              /* Pointer to a new mesh, call recursively */

              FillSTLMesh(stl, new, x0, y0, z0);
            }
          else if (new == 0)
            {
              /* No surface intersections, calculate center coordinates */
              /* NOTE: Choosing the center point is not a good idea,    */
              /* because cells with even number of sub-regions are      */
              /* split from there. */

              x = xmin + ((double)i + RandF(0))*px;
              y = ymin + ((double)j + RandF(0))*py;
              z = zmin + ((double)k + RandF(0))*pz;

              /* Check point */

              CheckValue(FUNCTION_NAME, "x", "", x, xmin, xmax);
              CheckValue(FUNCTION_NAME, "y", "", y, ymin, ymax);
              CheckValue(FUNCTION_NAME, "z", "", z, zmin, zmax);

              /* Calculate direction cosines */

              u = (x - x0);
              v = (y - y0);
              w = (z - z0);

              if ((l = sqrt(u*u + v*v + w*w)) == 0.0)
                Die(FUNCTION_NAME, "Division by zero");

              u = u/l;
              v = v/l;
              w = w/l;

              /* Store previous point */

              x0 = x;
              y0 = y;
              z0 = z;

              /* Find solid */

              sld = FindSTLSolid(stl, x, y, z, u, v, w, NO, 0);

              /* Allocate memory for mesh cell */

              ptr = NewItem(ptr, SEARCH_MESH_CELL_BLOCK_SIZE);

              /* Store */

              if (sld > VALID_PTR)
                WDB[ptr + SEARCH_MESH_CELL_CONTENT] = -(double)sld;

              /* Add to volume */

              WDB[stl + STL_SEARCH_MESH_V] =
                RDB[stl + STL_SEARCH_MESH_V] + px*py*pz;
            }
        }
}

/*****************************************************************************/
