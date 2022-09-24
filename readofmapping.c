/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readofmapping.c                                */
/*                                                                           */
/* Created:       2018/01/26 (VVa)                                           */
/* Last modified: 2018/01/26 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Reads OpenFOAM multi-physics interface output mapping file.  */
/*                                                                           */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadOFMapping:"

/*****************************************************************************/

void ReadOFMapping(long ifc)
{
  char *line, fname[MAX_STR];
  long n, nc, dim[7], i, cgns, disttype;
  double value;
  FILE *fp;

  /* Check output flag */

  if ((long)RDB[ifc + IFC_CALC_OUTPUT] == NO)
    return;

  /* Get output mapping file filename */

  sprintf(fname, "%s", GetText(ifc + IFC_PTR_OF_MAPFILE));

  /* Get number of cells */

  nc = (long)RDB[ifc + IFC_NC_PARENTS];

  /* Check if file is given */

  if (strcmp(fname, "-1"))
    {
      /* Open file for reading */

      if ((fp = fopen(fname, "r")) == NULL)
        Error(ifc, "Output map file \"%s\" does not exist",
              fname);

      /* Read header data */

      ReadOFHeader(fp, &n, &i, (long *)dim, &disttype, &value);

      /* Check size */

      if (i != nc)
        Error(ifc,
              "Invalid number of entries in output map file");

      /* Loop over cells */

      cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];
      while (cgns > VALID_PTR)
        {
          /* Read index */

          line = ReadOFData(fp, OF_FILE_MAP);

          if (sscanf(line, "%ld", &n) == EOF)
            Die(FUNCTION_NAME,
                "Not enough entries in output map file");

          /* Check value */

          if ((n < 1) || (n > nc))
            Error(ifc, "Invalid entry %ld in output map file", n);

          /* Put value */

          WDB[cgns + IFC_TET_PRNT_STAT_IDX] = (double)(n - 1);

          /* Next cell */

          cgns = NextItem(cgns);
        }

      /* Close file */

      fclose(fp);
    }
  else
    {
      /* Map each cell to itself */

      cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];
      while (cgns > VALID_PTR)
        {
          /* Copy index */

          WDB[cgns + IFC_TET_PRNT_STAT_IDX] =
            RDB[cgns + IFC_TET_PRNT_IDX];

          /* Next cell */

          cgns = NextItem(cgns);
        }
    }
}
