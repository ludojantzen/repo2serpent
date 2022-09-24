/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writewwmesh.c                                  */
/*                                                                           */
/* Created:       2018/11/08 (JLe)                                           */
/* Last modified: 2019/09/03 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Writes weight window mesh into a binary file.                */
/*                                                                           */
/* Comments: - Format was changed 8.11.2018 / 2.1.31 from ASCII to binary.   */
/*             Reading routine readwwmesh.c was changed as well.             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WriteWWMesh:"

/* Recursive write function */

void WriteWWMesh0(FILE *fp, long msh);

/*****************************************************************************/

void WriteWWMesh(long rmx)
{
  long i, j, k, n, nx, ny, nz, ng, nmax, ptr, loc0, loc1, msh, type, idx;
  FILE *fp;
  char outfile[MAX_STR];

  /* Pointer to mesh */

  msh = (long)RDB[rmx + RMX_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  /* Output file name */

  if ((long)RDB[DATA_RUN_VR_ITER] == NO)
    sprintf(outfile, "%s.wwd", GetText(DATA_PTR_INPUT_FNAME));
  else
    sprintf(outfile, "%s.wwd%ld", GetText(DATA_PTR_INPUT_FNAME),
            (long)RDB[DATA_VR_ITER_IDX]);

  /* Open data file */

  fp = fopen(outfile, "w");

  /* Print particle type */

  if ((type = (long)RDB[rmx + RMX_PARTICLE_TYPE]) == 0)
    Die(FUNCTION_NAME, "Particle type not set");
  else
    fwrite(&type, sizeof(long), 1, fp);

  /***************************************************************************/

  /***** Write energy grid data **********************************************/

  /* Print energy group data */

  if ((ptr = (long)RDB[rmx + RMX_PTR_EGRID]) > VALID_PTR)
    {
      /* Print number of groups */

      fwrite(&ng, sizeof(long), 1, fp);

      /* Pointer to group data */

      ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Print data */

      fwrite(&RDB[ptr], sizeof(double), ng + 1, fp);
    }
  else
    {
      /* No group structure defined */

      n = -1;
      fwrite(&n, sizeof(long), 1, fp);
    }

  /***************************************************************************/

  /***** Write data **********************************************************/

  /* Get size of data array */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  n = ListSize(loc0);
  fwrite(&n, sizeof(long), 1, fp);

  /* Loop over data */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Write mesh index and cell index */

      i = (long)RDB[loc0 + RMX_CELL_IDX];
      fwrite(&i, sizeof(long), 1, fp);

      i = (long)RDB[loc0 + RMX_CELL_MESH_IDX];
      fwrite(&i, sizeof(long), 1, fp);

      /* Write volume */

      fwrite(&RDB[loc0 + RMX_CELL_RVOL], sizeof(double), 1, fp);

      /* Get maximum number of neighbours */

      nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);
      fwrite(&nmax, sizeof(long), 1, fp);

      /* Loop over bounds */

      loc1 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
      while (loc1 > VALID_PTR)
        {
          /* Get pointer to cell */

          ptr = (long)RDB[loc1 + RMX_CELL_BOUND_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

          /* Write cell index */

          i = (long)RDB[ptr + RMX_CELL_IDX];
          fwrite(&i, sizeof(long), 1, fp);

          /* Write forward and adjoint indexes */

          i = (long)RDB[loc1 + RMX_CELL_BOUND_FWD_IDX];
          fwrite(&i, sizeof(long), 1, fp);

          i = (long)RDB[loc1 + RMX_CELL_BOUND_ADJ_IDX];
          fwrite(&i, sizeof(long), 1, fp);

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Source importances */

      ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
      CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr);
      fwrite(&RDB[ptr], sizeof(double), ng, fp);

      ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC_KEEP];
      CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr);
      fwrite(&RDB[ptr], sizeof(double), ng, fp);

      /* Current importances */

      ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr);
      fwrite(&RDB[ptr], sizeof(double), ng, fp);

      ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR_KEEP];
      CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr);
      fwrite(&RDB[ptr], sizeof(double), ng, fp);

      /* Partial solutions */

      ptr = (long)RDB[loc0 + RMX_CELL_ADJ_SOL_IN_CURR];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

      /* Write currents */

      fwrite(&RDB[ptr], sizeof(double), nmax*ng, fp);

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Write mesh **********************************************************/

  /* Get type */

  type = (long)RDB[msh + MESH_TYPE];
  fwrite(&type, sizeof(long), 1, fp);

  /* Check type */

  if (type == MESH_TYPE_ADAPTIVE)
    {
      /* Recursive write routine */

      WriteWWMesh0(fp, msh);
    }
  else
    {
      /* Get size */

      nx = (long)RDB[msh + MESH_N0];
      CheckValue(FUNCTION_NAME, "ng", "", nx, 1, 100000);

      ny = (long)RDB[msh + MESH_N1];
      CheckValue(FUNCTION_NAME, "ng", "", ny, 1, 100000);

      nz = (long)RDB[msh + MESH_N2];
      CheckValue(FUNCTION_NAME, "ng", "", nz, 1, 100000);

      /* Write mesh data */

      if ((type == MESH_TYPE_CARTESIAN) || (type == MESH_TYPE_CYLINDRICAL) ||
          (type == MESH_TYPE_HEXX) || (type == MESH_TYPE_HEXY))
        {
          /* Cartesian, cylindrical or hex mesh */

          fwrite(&nx, sizeof(long), 1, fp);
          fwrite(&ny, sizeof(long), 1, fp);
          fwrite(&nz, sizeof(long), 1, fp);

          fwrite(&RDB[msh + MESH_MIN0], sizeof(double), 1, fp);
          fwrite(&RDB[msh + MESH_MAX0], sizeof(double), 1, fp);
          fwrite(&RDB[msh + MESH_MIN1], sizeof(double), 1, fp);
          fwrite(&RDB[msh + MESH_MAX1], sizeof(double), 1, fp);
          fwrite(&RDB[msh + MESH_MIN2], sizeof(double), 1, fp);
          fwrite(&RDB[msh + MESH_MAX2], sizeof(double), 1, fp);
        }
      else if ((type == MESH_TYPE_ORTHOGONAL) || (type == MESH_TYPE_ICYL))
        {
          /* Unevenly-spaced orthogonal or cylindrical mesh, print sizes */

          fwrite(&nx, sizeof(long), 1, fp);
          fwrite(&ny, sizeof(long), 1, fp);
          fwrite(&nz, sizeof(long), 1, fp);

          /* Print x-boundaries */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_XLIM];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          fwrite(&RDB[ptr], sizeof(double), nx + 1, fp);

          /* Print y-boundaries */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_YLIM];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          fwrite(&RDB[ptr], sizeof(double), ny + 1, fp);

          /* Print z-boundaries */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_ZLIM];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          fwrite(&RDB[ptr], sizeof(double), nz + 1, fp);
        }
      else
        Die(FUNCTION_NAME, "Invalid mesh type");

      /* Loop over mesh and write mesh parameters */

      for (k = 0; k < nz; k++)
        for (j = 0; j < ny; j++)
          for (i = 0; i < nx; i++)
            {
              /* Pointer to data */

              loc0 = ReadMeshPtr(msh, i, j, k);
              CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

              loc0 = (long)RDB[loc0];
              CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

              /* Write cell index */

              idx = (long)RDB[loc0 + RMX_CELL_IDX];
              CheckValue(FUNCTION_NAME, "idx", "", idx, 0, 1000000000);
              fwrite(&idx, sizeof(long), 1, fp);
            }

      /***********************************************************************/
    }

  /* Close file */

  fclose(fp);
}

/*****************************************************************************/

/***** Recursive write function **********************************************/

void WriteWWMesh0(FILE *fp, long msh)
{
  long nx, ny, nz, i, j, k, idx, ptr, loc0;

  /* Get size */

  nx = (long)RDB[msh + MESH_N0];
  CheckValue(FUNCTION_NAME, "ng", "", nx, 1, 100000);

  ny = (long)RDB[msh + MESH_N1];
  CheckValue(FUNCTION_NAME, "ng", "", ny, 1, 100000);

  nz = (long)RDB[msh + MESH_N2];
  CheckValue(FUNCTION_NAME, "ng", "", nz, 1, 100000);

  /* Write size and boundaries (adaptive mesh is now alwasy Cartesian) */

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

          /* Check content */

          if ((loc0 = (long)RDB[ptr]) > VALID_PTR)
            {
              /* Write cell index */

              idx = (long)RDB[loc0 + RMX_CELL_IDX];
              CheckValue(FUNCTION_NAME, "idx", "", idx, 0, 1000000000);
              fwrite(&idx, sizeof(long), 1, fp);
            }
          else if (loc0 != NULLPTR)
            {
              /* Write -1 to indicate that another mesh follows */

              idx = -1;
              fwrite(&idx, sizeof(long), 1, fp);

              /* Call recursively */

              WriteWWMesh0(fp, -loc0);
            }
          else
            Die(FUNCTION_NAME, "WTF?");
        }
}

/*****************************************************************************/
