/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addsearchmesh.c                                */
/*                                                                           */
/* Created:       2012/02/18 (JLe)                                           */
/* Last modified: 2018/11/15 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Adds an item to search mesh                                  */
/*                                                                           */
/* Comments: - Handling of adaptivity depends on geometry type, currently    */
/*             works for unstructured tet mesh and STL.                      */
/*                                                                           */
/*           - NOTE: Tota kolmiofacettien lisäämistä tarkennettiin           */
/*                   15.11.2018 lisäämällä toi TestFacetOverlap(), joka      */
/*                   katsoo todelliset leikkaukset. Pitää vielä testata että */
/*                   toimii.                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddSearchMesh:"

/*****************************************************************************/

void AddSearchMesh(long msh, long loc0, double xmin, double xmax,
                   double ymin, double ymax, double zmin, double zmax)
{
  long nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax, i, j, k, loc1, ptr, new;
  long loc2, loc3, szp, dtype, ns, pts, n, m, side;
  double dx, dy, dz, lims[6], xmin0, xmax0, ymin0, ymax0, zmin0, zmax0;
  double xmin1, xmax1, ymin1, ymax1, zmin1, zmax1, x, y, z, x1, y1, z1, bb[6];
  double x2, y2, z2, x3, y3, z3, A, B, C, D;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Check mesh type and content */

  if (((long)RDB[msh + MESH_TYPE] != MESH_TYPE_CARTESIAN) &&
      ((long)RDB[msh + MESH_TYPE] != MESH_TYPE_ADAPTIVE))
    Die(FUNCTION_NAME, "Invalid mesh type");

  if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_PTR)
    Die(FUNCTION_NAME, "Invalid content type");

  /* Check order */

  if ((xmin > xmax) || (ymin > ymax) || (zmin > zmax))
    Die(FUNCTION_NAME, "Boundaries not in order %E %E %E %E %E %E", xmin, xmax, ymin, ymax, zmin, zmax);

  /* Get data type */

  dtype = (long)RDB[msh + MESH_CONTENT_DATA_TYPE];

  /* Get split criterion */

  ns = (long)RDB[msh + MESH_ADA_SPLIT];

  /* Get mesh size */

  nx = (long)RDB[msh + MESH_N0];
  ny = (long)RDB[msh + MESH_N1];
  nz = (long)RDB[msh + MESH_N2];

  /* Get limits */

  xmin0 = RDB[msh + MESH_MIN0];
  xmax0 = RDB[msh + MESH_MAX0];
  ymin0 = RDB[msh + MESH_MIN1];
  ymax0 = RDB[msh + MESH_MAX1];
  zmin0 = RDB[msh + MESH_MIN2];
  zmax0 = RDB[msh + MESH_MAX2];

  /* Calculate boundaries */

  imin = (long)(((double)nx)*(xmin - xmin0)/(xmax0 - xmin0));
  imax = (long)(((double)nx)*(xmax - xmin0)/(xmax0 - xmin0));
  jmin = (long)(((double)ny)*(ymin - ymin0)/(ymax0 - ymin0));
  jmax = (long)(((double)ny)*(ymax - ymin0)/(ymax0 - ymin0));
  kmin = (long)(((double)nz)*(zmin - zmin0)/(zmax0 - zmin0));
  kmax = (long)(((double)nz)*(zmax - zmin0)/(zmax0 - zmin0));

  /* Check with minimum and maximum */

  if (imax >= nx)
    imax = nx - 1;

  if (imin < 0)
    imin = 0;

  if (jmax >= ny)
    jmax = ny - 1;

  if (jmin < 0)
    jmin = 0;

  if (kmax >= nz)
    kmax = nz - 1;

  if (kmin < 0)
    kmin = 0;

  /* Calculate mesh cell size */

  dx = (xmax0 - xmin0)/((double)nx);
  dy = (ymax0 - ymin0)/((double)ny);
  dz = (zmax0 - zmin0)/((double)nz);

  /* Check content data type */

  if ((long)RDB[msh + MESH_CONTENT_DATA_TYPE] == MESH_CONTENT_DATA_STL)
    {
      /* Get surface coordinates */

      pts = (long)RDB[loc0 + STL_FACET_PTR_PT1];
      CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

      x1 = RDB[pts + STL_POINT_X];
      y1 = RDB[pts + STL_POINT_Y];
      z1 = RDB[pts + STL_POINT_Z];

      pts = (long)RDB[loc0 + STL_FACET_PTR_PT2];
      CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

      x2 = RDB[pts + STL_POINT_X];
      y2 = RDB[pts + STL_POINT_Y];
      z2 = RDB[pts + STL_POINT_Z];

      pts = (long)RDB[loc0 + STL_FACET_PTR_PT3];
      CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

      x3 = RDB[pts + STL_POINT_X];
      y3 = RDB[pts + STL_POINT_Y];
      z3 = RDB[pts + STL_POINT_Z];

      /* Calculate coefficients */

      A = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2);
      B = z2*x3 - z3*x2 - z1*(x3 - x2) + x1*(z3 - z2);
      C = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2);
      D = -x1*(y2*z3 - y3*z2) + y1*(x2*z3 - x3*z2) - z1*(x2*y3 - x3*y2);
    }
  else
    {
      /* Avoid compiler warning */

      A = 0.0;
      B = 0.0;
      C = 0.0;
      D = 0.0;
    }

  /* Loop over possible cells */

  for (i = imin; i < imax + 1; i++)
    for (j = jmin; j < jmax + 1; j++)
      for (k = kmin; k < kmax + 1; k++)
        {
          /*******************************************************************/
          
          /***** Exclusion by surface distance in STL ************************/

          /* Check content data type */

          if ((long)RDB[msh + MESH_CONTENT_DATA_TYPE] == MESH_CONTENT_DATA_STL)
            {
              /* Reset side */

              side = -1;

              /* Loop over corner points */

              for (n = 0; n < 8; n++)
                {
                  /* Set coordinates */

                  switch (n)
                    {
                    case 0:
                      {
                        x = xmin0;
                        y = ymin0;
                        z = zmin0;

                        break;
                      }

                    case 1:
                      {
                        x = xmax0;
                        y = ymin0;
                        z = zmin0;

                        break;
                      }

                    case 2:
                      {
                        x = xmin0;
                        y = ymax0;
                        z = zmin0;

                        break;
                      }

                    case 3:
                      {
                        x = xmin0;
                        y = ymin0;
                        z = zmax0;

                        break;
                      }

                    case 4:
                      {
                        x = xmax0;
                        y = ymax0;
                        z = zmin0;

                        break;
                      }

                    case 5:
                      {
                        x = xmax0;
                        y = ymin0;
                        z = zmax0;

                        break;
                      }

                    case 6:
                      {
                        x = xmin0;
                        y = ymax0;
                        z = zmax0;

                        break;
                      }

                    case 7:
                      {
                        x = xmax0;
                        y = ymax0;
                        z = zmax0;

                        break;
                      }

                    default:
                      {
                        /* Avoid compiler warning */

                        x = 0.0;
                        y = 0.0;
                        z = 0.0;

                        Die(FUNCTION_NAME, "WTF?");
                      }
                    }

                  /* Check side */

                  if (A*x + B*y + C*z + D > 0.0)
                    m = 1;
                  else
                    m = 0;

                  /* Compare */

                  if (n == 0)
                    side = m;
                  else if (side != m)
                    break;
                }

              /* Check count */

              if (n == 8)
                continue;
              else
                {
                  /* Calculate minimum coordinates */
                  
                  x = xmin0 + i*dx;
                  y = ymin0 + j*dy;
                  z = zmin0 + k*dz;
                  
                  /* Test overlap */
                  
                  if (TestFacetOverlap(loc0, x, x + dx, y, y + dy, z, z + dz)
                      == NO)
                    continue;
                }
            }
          
          /*******************************************************************/

          /***** Add content *************************************************/

          /* Get pointer */

          loc1 = ReadMeshPtr(msh, i, j, k);

          /* Check pointer */

          if ((new = -(long)RDB[loc1]) > VALID_PTR)
            {
              /* Mesh cell containts pointer to new mesh, call */
              /* recursively. */

              AddSearchMesh(new, loc0, xmin, xmax, ymin, ymax, zmin, zmax);

              /* Cycle loop */

              continue;
            }

          /* Check if buffer used for storing items is empty */

          if ((ptr = (long)RDB[DATA_ADA_MESH_BUF]) < VALID_PTR)
            NewItem(DATA_ADA_MESH_BUF, SEARCH_MESH_CELL_BLOCK_SIZE);

          /* Check list pointer */

          if ((long)RDB[loc1] < VALID_PTR)
            {
              /* No list created yet, must be created with NewItem() */

              loc2 = NewItem(loc1, SEARCH_MESH_CELL_BLOCK_SIZE);
            }
          else
            {
              /* Check if size is less than 2 (needs 1 item afer removing */
              /* to allow use of AddItem() later on) */

              if (ListSize(ptr) < 2)
                NewItem(DATA_ADA_MESH_BUF, SEARCH_MESH_CELL_BLOCK_SIZE);

              /* Get item from buffer */

              loc2 = (long)RDB[DATA_ADA_MESH_BUF];

              /* Remove and add to list */

              RemoveItem(loc2);
              AddItem(loc1, loc2);
            }

          /* Put pointer */

          WDB[loc2 + SEARCH_MESH_CELL_CONTENT] = (double)loc0;

          /* Allocate memory for counter (for STL meshes the pointer is */
          /* used for storing a list solids) */
          /* Only OMP id 0 will increment the counter to avoid excessive */
          /* memory size for search mesh */

          if ((long)RDB[msh + MESH_CONTENT_DATA_TYPE] == MESH_CONTENT_DATA_STL)
            WDB[loc2 + SEARCH_MESH_PTR_CELL_COUNT] = -1.0;
          else if (((long)RDB[loc2 + SEARCH_MESH_PTR_CELL_COUNT] < VALID_PTR)
                   && ((long)RDB[msh + MESH_CONTENT_DATA_TYPE] != 
                       MESH_CONTENT_DATA_TET))
            {
              ptr = AllocPrivateData(1, PRIVA_ARRAY);
              WDB[loc2 + SEARCH_MESH_PTR_CELL_COUNT] = (double)ptr;
            }

          /*******************************************************************/

          /***** Handle adaptivity *******************************************/

          /* Check type and split condition */

          if (((long)RDB[msh + MESH_TYPE] == MESH_TYPE_ADAPTIVE) &&
              (ListSize(loc2) >= (long)RDB[msh + MESH_ADA_SPLIT]))
            {
              /* Get pointer to size vector */

              szp = (long)RDB[msh + MESH_ADA_PTR_SZ];
              CheckPointer(FUNCTION_NAME, "(szp)", DATA_ARRAY, szp);

              /* Check size */

              if ((long)RDB[szp] < 0)
                continue;

              /* Calculate boundaries */

              lims[0] = xmin0 + (double)i*dx;
              lims[1] = lims[0] + dx;

              lims[2] = ymin0 + (double)j*dy;
              lims[3] = lims[2] + dy;

              lims[4] = zmin0 + (double)k*dz;
              lims[5] = lims[4] + dz;

              /* Create new mesh */

              new = CreateMesh(MESH_TYPE_ADAPTIVE, MESH_CONTENT_PTR, dtype,
                               ns, 0, 0, lims, szp);

              /* Remove items and add to new mesh */

              loc2 = (long)RDB[loc1];
              while (loc2 > VALID_PTR)
                {
                  /* Get content */

                  ptr = (long)RDB[loc2 + SEARCH_MESH_CELL_CONTENT];

                  /* Reset boundaries and avoid compiler warning */

                  xmin1 = INFTY;
                  xmax1 = -INFTY;
                  ymin1 = INFTY;
                  ymax1 = -INFTY;
                  zmin1 = INFTY;
                  zmax1 = -INFTY;

                  /* Check content data type */

                  if ((long)RDB[msh + MESH_CONTENT_DATA_TYPE] ==
                      MESH_CONTENT_DATA_TET)
                    {

                      /* Calculate limits based on points */

                      CalculateTetBoundingBox(ptr, bb);

                      /* Get pre-calculated bounding box */

                      xmin1 = bb[0];
                      xmax1 = bb[1];
                      ymin1 = bb[2];
                      ymax1 = bb[3];
                      zmin1 = bb[4];
                      zmax1 = bb[5];
                    }
                  else if ((long)RDB[msh + MESH_CONTENT_DATA_TYPE] ==
                      MESH_CONTENT_DATA_STL)
                    {
                      /* Calculate bounding box from points */

                      pts = (long)RDB[ptr + STL_FACET_PTR_PT1];
                      CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

                      if (RDB[pts + STL_POINT_X] < xmin1)
                        xmin1 = RDB[pts + STL_POINT_X];
                      if (RDB[pts + STL_POINT_X] > xmax1)
                        xmax1 = RDB[pts + STL_POINT_X];

                      if (RDB[pts + STL_POINT_Y] < ymin1)
                        ymin1 = RDB[pts + STL_POINT_Y];
                      if (RDB[pts + STL_POINT_Y] > ymax1)
                        ymax1 = RDB[pts + STL_POINT_Y];

                      if (RDB[pts + STL_POINT_Z] < zmin1)
                        zmin1 = RDB[pts + STL_POINT_Z];
                      if (RDB[pts + STL_POINT_Z] > zmax1)
                        zmax1 = RDB[pts + STL_POINT_Z];

                      pts = (long)RDB[ptr + STL_FACET_PTR_PT2];
                      CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

                      if (RDB[pts + STL_POINT_X] < xmin1)
                        xmin1 = RDB[pts + STL_POINT_X];
                      if (RDB[pts + STL_POINT_X] > xmax1)
                        xmax1 = RDB[pts + STL_POINT_X];

                      if (RDB[pts + STL_POINT_Y] < ymin1)
                        ymin1 = RDB[pts + STL_POINT_Y];
                      if (RDB[pts + STL_POINT_Y] > ymax1)
                        ymax1 = RDB[pts + STL_POINT_Y];

                      if (RDB[pts + STL_POINT_Z] < zmin1)
                        zmin1 = RDB[pts + STL_POINT_Z];
                      if (RDB[pts + STL_POINT_Z] > zmax1)
                        zmax1 = RDB[pts + STL_POINT_Z];

                      pts = (long)RDB[ptr + STL_FACET_PTR_PT3];
                      CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

                      if (RDB[pts + STL_POINT_X] < xmin1)
                        xmin1 = RDB[pts + STL_POINT_X];
                      if (RDB[pts + STL_POINT_X] > xmax1)
                        xmax1 = RDB[pts + STL_POINT_X];

                      if (RDB[pts + STL_POINT_Y] < ymin1)
                        ymin1 = RDB[pts + STL_POINT_Y];
                      if (RDB[pts + STL_POINT_Y] > ymax1)
                        ymax1 = RDB[pts + STL_POINT_Y];

                      if (RDB[pts + STL_POINT_Z] < zmin1)
                        zmin1 = RDB[pts + STL_POINT_Z];
                      if (RDB[pts + STL_POINT_Z] > zmax1)
                        zmax1 = RDB[pts + STL_POINT_Z];
                    }
                  else
                    Die(FUNCTION_NAME, "Invalid content data type");

                  /* Adjust boundaries */

                  xmin1 = xmin1 - 1E-6;
                  xmax1 = xmax1 + 1E-6;
                  ymin1 = ymin1 - 1E-6;
                  ymax1 = ymax1 + 1E-6;
                  zmin1 = zmin1 - 1E-6;
                  zmax1 = zmax1 + 1E-6;

                  /* Check */

                  if (xmax1 < xmin1)
                    Die(FUNCTION_NAME, "Error in bounding box");
                  if (ymax1 < ymin1)
                    Die(FUNCTION_NAME, "Error in bounding box");
                  if (zmax1 < zmin1)
                    Die(FUNCTION_NAME, "Error in bounding box");

                  /* Copy pointer and get pointer to next */

                  loc3 = loc2;
                  loc2 = NextItem(loc2);

                  /* Remove item and put to buffer */

                  RemoveItem(loc3);
                  AddItem(DATA_ADA_MESH_BUF, loc3);

                  /* Put to new mesh */

                  AddSearchMesh(new, ptr, xmin1, xmax1, ymin1, ymax1,
                                zmin1, zmax1);
                }

              /* Store pointer (negative value indicates that pointer is  */
              /* to another mesh. NOTE: Tässä on nyt sellainen vaara että */
              /* toi tallennettu arvo on sattumalta NULL) */

              if (new == -NULLPTR)
                Warn(FUNCTION_NAME, "null pointer here");

              WDB[loc1] = -(double)new;

              /* Check */

              ptr = ReadMeshPtr(msh, i, j, k);
              if ((long)RDB[ptr] != -new)
                Die(FUNCTION_NAME, "wtf?");
            }

          /*******************************************************************/
        }
}

/*****************************************************************************/
