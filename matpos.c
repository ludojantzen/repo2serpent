/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : matpos.c                                       */
/*                                                                           */
/* Created:       2018/05/24 (JLe)                                           */
/* Last modified: 2019/03/29 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Prints material name at given positions                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MatPos:"

/*****************************************************************************/

void MatPos()
{
  long cell, mat, nmax, n, uni, ptr;
  unsigned long seed;
  double x, y, z, u, v, w, f, T;
  char fname[MAX_STR];
  FILE *fp;

  /* Get number of points */

  if ((nmax = (long)RDB[DATA_MATPOS_N_PTS]) == 0)
    return;

  /* Loop over OpenMP threads (this is just to avoid errors) */

  for (n = 0; n < (long)RDB[DATA_OMP_MAX_THREADS]; n++)
    {
      /* Init random number sequence */

      seed = ReInitRNG(n);
      SEED[n*RNG_SZ] = seed;
    }

  /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */

  ExpandPrivateArrays();

  fprintf(outp, "Printing materials at given positions....\n\n");

  /* Check if file is given */

  if (nmax == -1)
    {
      /* Avoid compiler warning */

      fname[0] = '\0';

      /* Get file name */

      if ((long)RDB[DATA_MATPOS_PTR_FILE] < VALID_PTR)
        Die(FUNCTION_NAME, "No file name");
      else
        sprintf(fname, "%s", GetText(DATA_MATPOS_PTR_FILE));

      /* Open file */

      if ((fp = fopen(fname, "r")) == NULL)
        Error(0, "File \"%s\" does not exist", fname);

      /* Loop over file to count number of points */

      nmax = 0;
      while(fscanf(fp, "%lf", &x) != EOF)
        nmax++;

      /* Check */

      if (nmax % 3)
        Error(0, "Invalid number of coordinates given");

      /* Allocate memory for points */

      ptr = ReallocMem(DATA_ARRAY, nmax);
      WDB[DATA_MATPOS_PTR_COORD] = (double)ptr;

      /* Rewind */

      rewind(fp);

      /* Loop over file and read data */

      while(fscanf(fp, "%lf", &WDB[ptr++]) != EOF);

      /* Close file */

      fclose(fp);

      /* Put number of points */

      nmax = nmax/3;
    }

  /* Get pointer to data */

  ptr = (long)RDB[DATA_MATPOS_PTR_COORD];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Print header */

  fprintf(outp, "           x            y            z");
  fprintf(outp, "               universe                  cell             material");
  fprintf(outp, "\n");

  /* Loop over points */

  for (n = 0; n < nmax; n++)
    {
      /* Get point */

      x = RDB[ptr++];
      y = RDB[ptr++];
      z = RDB[ptr++];

      /* Sample direction (this is necessary for STL geometries) */

      IsotropicDirection(&u, &v, &w, 0);

      /* Find position */

      if ((cell = WhereAmI(x, y, z, u, v, w, 0)) < 0)
        Error(0, "Geometry error at %E %E %E", x, y, z);

      /* Get universe */

      uni = (long)RDB[cell + CELL_PTR_UNI];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Reset density and temperature */

      f = 1.0;
      T = 0.0;

      /* Check if cell has material */

      if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
        {
          /* Get material pointer */

          mat = MatPtr(mat, 0);

          /* Get point from interface */

          IFCPoint(mat, &f, &T, -1.0, 0);
        }

      /* Print */

      fprintf(outp, "%12.5E %12.5E %12.5E  ", x, y, z);

      /* NOTE: STL-universumeilla ei oo nimiä? */

      if ((long)RDB[uni + UNIVERSE_PTR_NAME] < VALID_PTR)
        fprintf(outp, " %20s ", "N/A");
      else
        fprintf(outp, " %20s ", GetText(uni + UNIVERSE_PTR_NAME));


      fprintf(outp, " %20s ", GetText(cell + CELL_PTR_NAME));

      if ((long)RDB[cell + CELL_TYPE] == CELL_TYPE_OUTSIDE)
        fprintf(outp, "%20s  ", "outside");
      else if (mat < VALID_PTR)
        fprintf(outp, "%20s  ", "void");
      else
        fprintf(outp, "%20s  ", GetText(mat + MATERIAL_PTR_NAME));

      fprintf(outp, "\n");
    }

  /* Exit subroutine */

  fprintf(outp, "\nOK.\n\n");

  /* Terminate run */

  exit(-1);

  /***************************************************************************/
}

/*****************************************************************************/
