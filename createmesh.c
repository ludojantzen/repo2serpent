/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : createmesh.c                                   */
/*                                                                           */
/* Created:       2011/05/14 (JLe)                                           */
/* Last modified: 2019/03/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Creates a new mesh structure                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CreateMesh:"

/*****************************************************************************/

long CreateMesh(long type, long cont, long dtype, long nx, long ny, long nz,
                const double *params, long szp)
{
  long msh, ptr, n, m;
  double xmin, xmax, ymin, ymax, zmin, zmax;

  /* Avoid compiler warning */

  xmin = INFTY;
  xmax = -INFTY;
  ymin = INFTY;
  ymax = -INFTY;
  zmin = INFTY;
  zmax = -INFTY;

  /* Allocate memory for structure */

  msh = ReallocMem(DATA_ARRAY, MESH_BLOCK_SIZE);

  /* NOTE: This is to avoid accidental null pointer in STL geometries */
  /* (huom! tää on todella typerä tapa hoitaa tää ongelma) */

  if (msh == -NULLPTR)
    msh = ReallocMem(DATA_ARRAY, MESH_BLOCK_SIZE);

  /* Check type */

  if ((type == MESH_TYPE_CARTESIAN) || (type == MESH_TYPE_HEXX) ||
      (type == MESH_TYPE_HEXY))
    {
      /***********************************************************************/

      /***** Cartesian and hex meshes ****************************************/

      /* Check dimensions */

      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 100000);
      CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 100000);
      CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 100000);

      /* Get parameters */

      xmin = params[0];
      xmax = params[1];
      ymin = params[2];
      ymax = params[3];
      zmin = params[4];
      zmax = params[5];

      /* Check */

      if (type == MESH_TYPE_CARTESIAN)
        {
          CheckValue(FUNCTION_NAME, "xmin", "", xmin, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "xmax", "", xmax, xmin, INFTY);
          CheckValue(FUNCTION_NAME, "ymin", "", ymin, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "ymax", "", ymax, ymin, INFTY);
          CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);
        }
      else
        {
          CheckValue(FUNCTION_NAME, "xmin", "", xmin, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "xmax", "", xmax, 0.0, INFTY);
          CheckValue(FUNCTION_NAME, "ymin", "", ymin, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "ymax", "", ymax, 0.0, 0.0);
          CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);
        }

      /***********************************************************************/
    }
   else if (type == MESH_TYPE_ADAPTIVE)
    {
      /***********************************************************************/

      /***** Adaptive mesh ***************************************************/

      /* Check dimensions (split criterion is in nx, ny and nz are zero) */

      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, INFTY);
      CheckValue(FUNCTION_NAME, "ny", "", ny, 0, 0);
      CheckValue(FUNCTION_NAME, "nz", "", nz, 0, 0);

      /* Get parameters */

      xmin = params[0];
      xmax = params[1];
      ymin = params[2];
      ymax = params[3];
      zmin = params[4];
      zmax = params[5];

      /* Check */

      CheckValue(FUNCTION_NAME, "xmin", "", xmin, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "xmax", "", xmax, xmin, INFTY);
      CheckValue(FUNCTION_NAME, "ymin", "", ymin, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "ymax", "", ymax, ymin, INFTY);
      CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);

      /* Check size poiner */

      CheckPointer(FUNCTION_NAME, "(szp)", DATA_ARRAY, szp);

      /* Put split criterion and pointer to next level size */

      WDB[msh + MESH_ADA_SPLIT] = (double)nx;
      WDB[msh + MESH_ADA_PTR_SZ] = (double)(szp + 1);

      /* Get size */

      nx = (long)RDB[szp];
      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 1000000);

      /* Put size */

      WDB[msh + MESH_N0] = (double)nx;
      WDB[msh + MESH_N1] = (double)nx;
      WDB[msh + MESH_N2] = (double)nx;

      /* Set ny and nz equal to nx */

      ny = nx;
      nz = nx;

      /***********************************************************************/
    }
  else if (type == MESH_TYPE_CYLINDRICAL)
    {
      /***********************************************************************/

      /***** Cylindrical mesh ************************************************/

      /* Check dimensions */

      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 100000);
      CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 100000);
      CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 100000);

      /* Get parameters */

      xmin = params[0];
      xmax = params[1];
      ymin = params[2];
      ymax = params[3];
      zmin = params[4];
      zmax = params[5];

      /* Check */

      CheckValue(FUNCTION_NAME, "xmin", "", xmin, 0.0, INFTY);
      CheckValue(FUNCTION_NAME, "xmax", "", xmax, xmin, INFTY);
      CheckValue(FUNCTION_NAME, "ymin", "", ymin, 0.0, PI);
      CheckValue(FUNCTION_NAME, "ymax", "", ymax, ymin, 2*PI);
      CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);

      /***********************************************************************/
    }
  else if (type == MESH_TYPE_ORTHOGONAL)
    {
      /***********************************************************************/

      /***** Unevenly-spaced orthogolan mesh *********************************/

      /* Check dimensions */

      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 100000);
      CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 100000);
      CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 100000);

      /* Reset index */

      n = 0;

      /* Allocate memory for x-limits */

      ptr = ReallocMem(DATA_ARRAY, nx + 1);
      WDB[msh + MESH_ORTHO_PTR_XLIM] = (double)ptr;

      /* Read values */

      for (m = 0; m < nx + 1; m++)
        {
          /* Put value */

          WDB[ptr++] = params[n + m];

          /* Check order */

          if ((nx > 1) && (m > 0))
            if (params[n + m - 1] >= params[n + m])
              Die(FUNCTION_NAME, "Error in x-limits");
        }

      /* Put minimum and maximum */

      xmin = params[n];
      xmax = params[n + nx];

      /* Update index */

      n = n + nx + 1;

      /* Allocate memory for y-limits */

      ptr = ReallocMem(DATA_ARRAY, ny + 1);
      WDB[msh + MESH_ORTHO_PTR_YLIM] = (double)ptr;

      /* Read values */

      for (m = 0; m < ny + 1; m++)
        {
          /* Put value */

          WDB[ptr++] = params[n + m];

          /* Check order */

          if ((ny > 1) && (m > 0))
            if (params[n + m - 1] >= params[n + m])
              Die(FUNCTION_NAME, "Error in y-limits");
        }

      /* Put minimum and maximum */

      ymin = params[n];
      ymax = params[n + ny];

      /* Update index */

      n = n + ny + 1;

      /* Allocate memory for z-limits */

      ptr = ReallocMem(DATA_ARRAY, nz + 1);
      WDB[msh + MESH_ORTHO_PTR_ZLIM] = (double)ptr;

      /* Read values */

      for (m = 0; m < nz + 1; m++)
        {
          /* Put value */

          WDB[ptr++] = params[n + m];

          /* Check order */

          if ((nz > 1) && (m > 0))
            if (params[n + m - 1] >= params[n + m])
              Die(FUNCTION_NAME, "Error in z-limits");
        }

      /* Put minimum and maximum */

      zmin = params[n];
      zmax = params[n + nz];

      /***********************************************************************/
    }
  else if (type == MESH_TYPE_ICYL)
    {
      /***********************************************************************/

      /***** Unevenly-spaced cylindrical mesh ********************************/

      /* Check dimensions */

      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 100000);
      CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 100000);
      CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 100000);

      /* Reset index */

      n = 0;

      /* Allocate memory for radial limits */

      ptr = ReallocMem(DATA_ARRAY, nx + 1);
      WDB[msh + MESH_ORTHO_PTR_XLIM] = (double)ptr;

      /* Read values */

      for (m = 0; m < nx + 1; m++)
        {
          /* Check value */

          CheckValue(FUNCTION_NAME, "r", "", params[n + m], 0.0, INFTY);

          /* Put value */

          WDB[ptr++] = params[n + m];

          /* Check order */

          if ((nx > 1) && (m > 0))
            if (params[n + m - 1] >= params[n + m])
              Die(FUNCTION_NAME, "Error in radial limits");
        }

      /* Put minimum and maximum */

      xmin = params[n];
      xmax = params[n + nx];

      /* Update index */

      n = n + nx + 1;

      /* Allocate memory for angular limits */

      ptr = ReallocMem(DATA_ARRAY, ny + 1);
      WDB[msh + MESH_ORTHO_PTR_YLIM] = (double)ptr;

      /* Read values */

      for (m = 0; m < ny + 1; m++)
        {
          /* Check value */

          CheckValue(FUNCTION_NAME, "phi", "", params[n + m], 0.0, 2.001*PI);

          /* Put value */

          if (params[n + m] < 2.0*PI)
            WDB[ptr++] = params[n + m];
          else
            WDB[ptr++] = 2.0*PI;

          /* Check order */

          if ((ny > 1) && (m > 0))
            if (params[n + m - 1] >= params[n + m])
              Die(FUNCTION_NAME, "Error in angular limits");
        }

      /* Put minimum and maximum */

      ymin = params[n];
      ymax = params[n + ny];

      /* Update index */

      n = n + ny + 1;

      /* Allocate memory for axial limits */

      ptr = ReallocMem(DATA_ARRAY, nz + 1);
      WDB[msh + MESH_ORTHO_PTR_ZLIM] = (double)ptr;

      /* Read values */

      for (m = 0; m < nz + 1; m++)
        {
          /* Check value */

          CheckValue(FUNCTION_NAME, "z", "", params[n + m], -INFTY, INFTY);

          /* Put value */

          WDB[ptr++] = params[n + m];

          /* Check order */

          if ((nz > 1) && (m > 0))
            if (params[n + m - 1] >= params[n + m])
              Die(FUNCTION_NAME, "Error in axial limits");
        }

      /* Put minimum and maximum */

      zmin = params[n];
      zmax = params[n + nz];

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid mesh type %ld", type);

  /* Put common parameters */

  WDB[msh + MESH_TYPE] = (double)type;
  WDB[msh + MESH_CONTENT] = (double)cont;
  WDB[msh + MESH_CONTENT_DATA_TYPE] = (double)dtype;
  WDB[msh + MESH_MIN0] = xmin;
  WDB[msh + MESH_MAX0] = xmax;
  WDB[msh + MESH_MIN1] = ymin;
  WDB[msh + MESH_MAX1] = ymax;
  WDB[msh + MESH_MIN2] = zmin;
  WDB[msh + MESH_MAX2] = zmax;
  WDB[msh + MESH_N0] = (double)nx;
  WDB[msh + MESH_N1] = (double)ny;
  WDB[msh + MESH_N2] = (double)nz;
  AllocValuePair(msh + MESH_PREV_COL_IDX);

  /* Check content */

  if (cont == MESH_CONTENT_RES)
    {
      /* Allocate memory for results */

      ptr = AllocPrivateData(nx*ny*nz + 1, RES2_ARRAY);

      /* Put pointer */

      WDB[msh + MESH_PTR_RES2] = (double)ptr;
    }
  else if (cont == MESH_CONTENT_DAT)
    {
      /* Allocate memory for data */

      ptr = ReallocMem(DATA_ARRAY, nx*ny*nz + 1);

      /* Put pointer */

      WDB[msh + MESH_PTR_DATA] = (double)ptr;
    }
  else if (cont == MESH_CONTENT_PTR)
    {
      /* Allocate memory for pointer */

      ptr = ReallocMem(DATA_ARRAY, nx*ny*nz);

      /* Put pointer */

      WDB[msh + MESH_PTR_PTR] = (double)ptr;
    }
  else if (cont == MESH_CONTENT_NONE)
    {
      /* No content */

      WDB[msh + MESH_PTR_PTR] = NULLPTR;
    }
  else
    Die(FUNCTION_NAME, "Invalid content type");

  /* Return mesh pointer */

  return msh;
}

/*****************************************************************************/
