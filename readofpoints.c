/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readofpoints.c                                 */
/*                                                                           */
/* Created:       2018/01/26 (VVa)                                           */
/* Last modified: 2018/01/26 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Reads OpenFOAM multi-physics interface points file           */
/*                                                                           */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadOFPoints:"

/*****************************************************************************/

void ReadOFPoints(long ifc)
{
  double x, y, z, xmin, xmax, ymin, ymax, zmin, zmax, value;
  char *line, fname[MAX_STR];
  long n, nd, dim[7], pointlist, disttype;
  FILE *fp;

  /***** Read points *************************************************/

  /* Reset mesh boundaries */

  xmin =  INFTY;
  xmax = -INFTY;
  ymin =  INFTY;
  ymax = -INFTY;
  zmin =  INFTY;
  zmax = -INFTY;

  if ((long)RDB[ifc + IFC_PTR_OF_PFILE] < VALID_PTR)
    Error(ifc, "No points file defined");

  /* Get points file filename */

  sprintf(fname, "%s", GetText(ifc + IFC_PTR_OF_PFILE));

  /* Check file format */

  TestDOSFile(fname);

  /* Open points file for reading */

  if ((fp = fopen(fname, "r")) == NULL)
    Error(ifc, "Points file \"%s\" does not exist",
          fname);

  /* Read header data */

  ReadOFHeader(fp, &n, &nd, (long *)dim, &disttype, &value);

  /* Check number of points */

  CheckValue(FUNCTION_NAME, "nd", "", nd, 4, 100000000000);

  /* Store number of points */

  WDB[ifc + IFC_NP] = (double)nd;

  /* Allocate memory for points */

  pointlist = ReallocMem(DATA_ARRAY, nd*3);

  /* Store pointer to first point */

  WDB[ifc + IFC_PTR_POINT_LIST] = (double)pointlist;
  WDB[ifc + IFC_PTR_POINT_LIST_PARENTS] = (double)pointlist;

  /* Read points */

  for (n = 0; n < nd; n++)
    {
      /* Read coordinates */

      line = ReadOFData(fp, OF_FILE_POINTS);

      if (sscanf(line, "%lf %lf %lf", &x, &y, &z) == EOF)
        Error(ifc, "Not enough entries in points file");

      /* Convert to cm */

      x = x*100.0;
      y = y*100.0;
      z = z*100.0;

      /* Put data */

      WDB[pointlist++] = x;
      WDB[pointlist++] = y;
      WDB[pointlist++] = z;

      /* Compare to limits */

      if (x < xmin)
        xmin = x;
      if (x > xmax)
        xmax = x;

      if (y < ymin)
        ymin = y;
      if (y > ymax)
        ymax = y;

      if (z < zmin)
        zmin = z;
      if (z > zmax)
        zmax = z;
    }

  /* Close file */

  fclose(fp);

  /* Store boundaries */

  WDB[ifc + IFC_MESH_XMIN] = xmin;
  WDB[ifc + IFC_MESH_XMAX] = xmax;
  WDB[ifc + IFC_MESH_YMIN] = ymin;
  WDB[ifc + IFC_MESH_YMAX] = ymax;
  WDB[ifc + IFC_MESH_ZMIN] = zmin;
  WDB[ifc + IFC_MESH_ZMAX] = zmax;

}
