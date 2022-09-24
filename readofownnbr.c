/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readofownnbr.c                                 */
/*                                                                           */
/* Created:       2018/01/26 (VVa)                                           */
/* Last modified: 2018/01/26 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Reads OpenFOAM multi-physics interface owner/neighbour files */
/*                                                                           */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadOFOwnNbr:"

/*****************************************************************************/

void ReadOFOwnNbr(long ifc)
{
  long n, nf, fi, nmax, dim[7], solist, snlist, i, nc, disttype;
  double value;
  char *line, ofname[MAX_STR], nfname[MAX_STR];
  FILE *fp;

  /* Get owner file filename */

  sprintf(ofname, "%s", GetText(ifc + IFC_PTR_OF_OFILE));

  /* Get neighbour file filename */

  sprintf(nfname, "%s", GetText(ifc + IFC_PTR_OF_NFILE));

  /* Check file formats */

  TestDOSFile(ofname);
  TestDOSFile(nfname);

  /* Get number of faces */

  nf = (long)RDB[ifc + IFC_NF_PARENTS];

  /* Get owner and neighbor lists */

  solist = (long)RDB[ifc + IFC_PTR_OWNR_LIST_PARENTS];
  snlist = (long)RDB[ifc + IFC_PTR_NBR_LIST_PARENTS];

  /* Reset number of cells */

  nc = -1;

  /* Repeat for owner and neighbour files */

  for (fi = 1; fi < 3; fi++)
    {
      /* Check mode */

      if (fi == 1)
        {
          /* Open file */

          if ((fp = fopen(ofname, "r"))
              == NULL)
            Error(ifc, "Owner file \"%s\" does not exist",
                  ofname);

          /* Read header data */

          ReadOFHeader(fp, &n, &nmax, (long *)dim, &disttype, &value);

          /* Check number of faces */

          if (nmax > nf)
            Die(FUNCTION_NAME, "Number of faces linked to owners in file "
                "%s is %ld, larger than the number of defined faces (%ld)\n", ofname, nmax, nf);
        }
      else
        {
          /* Open file */

          if ((fp = fopen(nfname, "r"))
              == NULL)
            Error(ifc, "Neighbour file \"%s\" does not exist",
                  nfname);

          /* Read header data */

          ReadOFHeader(fp, &n, &nmax, (long *)dim, &disttype, &value);

          /* Check number of faces */

          if (nmax > nf)
            Die(FUNCTION_NAME, "Number of faces linked to neighbours in file "
                "%s is %ld, larger than the number of defined faces (%ld)\n", nfname, nmax, nf);
        }

      /* Loop over faces */

      for (n = 0; n < nmax; n++)
        {
          /* Read cell index from owner/neighbour file */

          if (fi == 1)
            line = ReadOFData(fp, OF_FILE_OWNER);
          else
            line = ReadOFData(fp, OF_FILE_NEIGHBOUR);

          if (sscanf(line, "%ld", &i) == EOF)
            {
              if (fi == 1)
                Error(ifc, "Not enough entries in owner file");
              else
                Error(ifc, "Not enough entries in neighbour file");
            }

          /* Update number of cells */

          if (i + 1 > nc)
            nc = i + 1;

          /* Put pointer in array */

          if(fi == 1)
            WDB[solist + n] = (double)i;
          else
            WDB[snlist + n] = (double)i;

        }

      /* Close file */

      fclose(fp);
    }

  /* Put number of cells */

  WDB[ifc + IFC_NC] = (double)nc;
}
