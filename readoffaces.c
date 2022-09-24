/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readoffaces.c                                  */
/*                                                                           */
/* Created:       2018/01/26 (VVa)                                           */
/* Last modified: 2018/01/26 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Reads OpenFOAM multi-physics interface faces  file           */
/*                                                                           */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadOFFaces:"

/*****************************************************************************/

void ReadOFFaces(long ifc)
{
  char *line, tmpstr[MAX_STR], fname[MAX_STR];
  long n, nd, dim[7], pointlist, solist, snlist, sslist, surf, ptr, i, j;
  long p, np, nf, disttype;
  double value;
  FILE *fp;

  /* Get faces file filename */

  sprintf(fname, "%s", GetText(ifc + IFC_PTR_OF_FFILE));

  /* Check file format */

  TestDOSFile(fname);

  /* Open faces file for reading */

  if ((fp = fopen(fname, "r")) == NULL)
    Error(ifc, "Faces file \"%s\" does not exist",
          fname);

  /* Read header data */

  ReadOFHeader(fp, &n, &nf, (long *)dim, &disttype, &value);

  /* Check number of faces */

  CheckValue(FUNCTION_NAME, "nf", "", nf, 4, 100000000000);

  /* Store number of faces */

  WDB[ifc + IFC_NF_PARENTS] = (double)nf;

  /* Preallocate memory for data */

  PreallocMem(nf*(UMSH_SURF_BLOCK_SIZE + 3) + 2*nf, DATA_ARRAY);

  /* Allocate memory for pointers */

  solist = ReallocMem(DATA_ARRAY, nf);
  WDB[ifc + IFC_PTR_OWNR_LIST_PARENTS] = (double)solist;

  snlist = ReallocMem(DATA_ARRAY, nf);
  WDB[ifc + IFC_PTR_NBR_LIST_PARENTS] = (double)snlist;

  sslist = ReallocMem(DATA_ARRAY, nf);
  WDB[ifc + IFC_PTR_SURF_LIST_PARENTS] = (double)sslist;

  /* Reset surface pointer */

  surf = -1;

  /* Get number of points */

  nd = (long)RDB[ifc + IFC_NP];

  /* Get pointer to first point */

  pointlist = (long)RDB[ifc + IFC_PTR_POINT_LIST];

  /* Loop over faces */

  for (n = 0; n < nf; n++)
    {

      /* Create surface and put pointer to its place */

      surf = ReallocMem(DATA_ARRAY, UMSH_SURF_BLOCK_SIZE);

      /* Set/reset pointers */

      WDB[sslist + n] = surf;
      WDB[solist + n] = -1;
      WDB[snlist + n] = -1;

      /* Read entry */

      line = ReadOFData(fp, OF_FILE_FACES);

      /* Read number of points */

      p = NextWord(line, tmpstr);
      line = &line[p];
      np = (long)atoi(tmpstr);

      /* Check type */

      if (np < 3)
        Error(ifc, "Not enough points");

      /* Store number of points */

      WDB[surf + UMSH_SURF_N_POINTS] = (double)np;

      /* Allocate memory for points */

      ptr = ReallocMem(DATA_ARRAY, np);
      WDB[surf + UMSH_SURF_PTR_POINTS] = (double)ptr;

      /* Read points */

      for (i = 0; i < np; i++)
        {
          /* Read point index */

          p = NextWord(line, tmpstr);
          line = &line[p];
          j = (long)atoi(tmpstr);

          /* Check */

          if ((j < 0) || (j > nd - 1))
            Die(FUNCTION_NAME, "Invalid point index %ld", j);

          /* Store pointer to beginning of point in point list */

          WDB[ptr++] = (double)(pointlist + j*3);

        }
    }

  /* Close file */

  fclose(fp);

#ifdef mmmaaa
  /* Loop over faces and print out */

  pointlist = (long)RDB[ifc + IFC_PTR_POINT_LIST];

  for (n = 0; n < nf; n++)
    {
      surf = (long)RDB[sslist + n];

      printf("\nFace:\n");

      /* Get number of points */

      np = (long)RDB[surf + UMSH_SURF_N_POINTS];

      /* Get pointer to surface point list */

      ptr = (long)RDB[surf + UMSH_SURF_PTR_POINTS];

      /* Loop over points */

      for (i = 0; i < np; i++)
        {
          /* Get point */
          pt = (long)RDB[ptr + i];

          /* Print point */
          printf("Point: %f %f %f\n", RDB[pt], RDB[pt + 1], RDB[pt + 2]);
        }
    }

#endif

}
