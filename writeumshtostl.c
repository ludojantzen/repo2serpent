/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writeumshtostl.c                               */
/*                                                                           */
/* Created:       2014/03/12 (JLe)                                           */
/* Last modified: 2018/01/26 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Used to write tet-based unstructured mesh geometry into STL. */
/*              To be used for testing and debugging purposes only.          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WriteUMSHtoSTL:"

/*****************************************************************************/

void WriteUMSHtoSTL()
{
#ifdef blohblah
  long loc0, loc1, loc2, cell, surf, ptr, prnt;
  char mname[MAX_STR], fname[MAX_STR];
  double x0, y0, z0, x1, y1, z1, x2, y2, z2, u0, v0, w0, u1, v1, w1;
  double u2, v2, w2, l;
  FILE *fp;

  if ((long)RDB[DATA_PTR_UMSH0] < VALID_PTR)
    return;

  fprintf(outp, "Writing unstructured mesh geometries into STL file...\n");

  /* Loop over geometries */

  loc0 = (long)RDB[DATA_PTR_UMSH0];
  while (loc0 > VALID_PTR)
    {
      /* Get pointer to interface structure */

      loc1 = (long)RDB[loc0 + UMSH_PTR_IFC];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      loc1 = (long)RDB[loc1 + IFC_PTR_TET_MSH];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Avoid compiler warning */

      cell = -1;

      /* Loop over cells to find first material name */

      while (loc1 > VALID_PTR)
        {
          /* Get pointer to parent cell (if exists) or self */

          if ((prnt = (long)RDB[loc1 + IFC_TET_MSH_PTR_PARENT]) < VALID_PTR)
            prnt = loc1;

          /* Get pointer to cell structure */

          cell = (long)RDB[prnt + IFC_TET_PRNT_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Check material pointer */

          if ((cell = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
            break;

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Check pointer and get material name */

      if (loc1 < VALID_PTR)
        Error(loc0, "Geometry has no materials");
      else
        sprintf(mname, "%s", GetText(cell + MATERIAL_PTR_NAME));

      /* Get file name */

      sprintf(fname, "%s.stl", mname);

      /* Open file for writing */

      if ((fp = fopen(fname, "w")) == NULL)
        Die(FUNCTION_NAME, "Unable to open file for writing");

      /* Print title */

      fprintf(fp, " solid %s\n", mname);

      /* Get pointer to interface structure */

      loc1 = (long)RDB[loc0 + UMSH_PTR_IFC];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      loc1 = (long)RDB[loc1 + IFC_PTR_TET_MSH];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop over cells */

      while (loc1 > VALID_PTR)
        {
          /* Get pointer to parent cell (if exists) or self */

          if ((prnt = (long)RDB[loc1 + IFC_TET_MSH_PTR_PARENT]) < VALID_PTR)
            prnt = loc1;

          /* Get pointer to cell structure */

          cell = (long)RDB[prnt + IFC_TET_PRNT_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Pointer to intersection list */

          loc2 = (long)RDB[cell + CELL_PTR_SURF_INSC];
          while (loc2 > VALID_PTR)
            {
              /* Check if surface has neighbours */

              if ((long)RDB[loc2 + CELL_INSC_PTR_NEXT_TET_CELL] > VALID_PTR)
                {
                  /* Skip */

                  loc2 = NextItem(loc2);

                  /* Cycle loop */

                  continue;
                }

              /* Pointer to surface */

              surf = (long)RDB[loc2 + CELL_INSC_PTR_SURF];
              CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

              /* Get pointer to parameters */

              ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get coordinates */

              x0 = RDB[ptr++];
              y0 = RDB[ptr++];
              z0 = RDB[ptr++];
              x1 = RDB[ptr++];
              y1 = RDB[ptr++];
              z1 = RDB[ptr++];
              x2 = RDB[ptr++];
              y2 = RDB[ptr++];
              z2 = RDB[ptr++];

              /* Calculate surface normal */

              /* Calculate vectors */

              u1 = x1 - x0;
              v1 = y1 - y0;
              w1 = z1 - z0;

              u2 = x2 - x1;
              v2 = y2 - y1;
              w2 = z2 - z1;

              /* Cross product */

              u0 =  (v1*w2 - w1*v2);
              v0 = -(u1*w2 - w1*u2);
              w0 =  (u1*v2 - v1*u2);

              /* Normalize */

              if ((l = u0*u0 + u0*u0 + u0*u0) > 0.0)
                {
                  l = sqrt(l);
                  u0 = u0/l;
                  v0 = v0/l;
                  w0 = w0/l;
                }
              else
                Note(loc0, "Degenerate facet encountered");

              /* Print */

              fprintf(fp, "facet normal %13.6E %13.6E %13.6E\n", u0, v0, w0);
              fprintf(fp, " outer loop\n");
              fprintf(fp, "  vertex %13.6E %13.6E %13.6E\n", x0, y0, z0);
              fprintf(fp, "  vertex %13.6E %13.6E %13.6E\n", x1, y1, z1);
              fprintf(fp, "  vertex %13.6E %13.6E %13.6E\n", x2, y2, z2);
              fprintf(fp, " endloop\n");
              fprintf(fp, "endfacet\n");

              /* Next */

              loc2 = NextItem(loc2);
            }

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Print end */

      fprintf(fp, " endsolid %s\n", mname);

      /* Close file */

      fclose(fp);

      /* Next geometry */

      loc0 = NextItem(loc0);
    }

  fprintf(outp, "OK.\n\n");

  /* Terminate run */

  exit(0);
#endif
}

/*****************************************************************************/
