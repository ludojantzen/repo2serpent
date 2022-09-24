/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readofheader.c                                 */
/*                                                                           */
/* Created:       2013/12/27 (JLe)                                           */
/* Last modified: 2018/01/26 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Reads header data from an OpenFOAM format file               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ReadOFHeader:"

/*****************************************************************************/

void ReadOFHeader(FILE *fp, long *type, long *sz, long *dim, long *disttype, double *value)
{
  long n, nread;
  char tmpstr[MAX_STR], tmpstr2[MAX_STR];
  double unival;
  signed int c;

  /* Check file pointer */

  if (fp == NULL)
    Die(FUNCTION_NAME, "Pointer error");

  /* Reset type */

  *type = -1;
  *disttype = OF_INTERNAL_FIELD_UNKNOWN;

  /* Reset dimension vector */

  for (n = 0; n < 7; n++)
    dim[n] = 0;

  /* Read first string */

  if (fscanf(fp, "%s", tmpstr) == EOF)
    Die(FUNCTION_NAME, "Unexpected EOF1");

  /* Check size */

  if ((n = atol(tmpstr)) > 0)
    {
      /* Put size */

      *sz = n;

      /* Exit subroutine */

      return;
    }

  /* Look for key word "object" */

  do
    {
      /* Check string */

      if (!strcmp(tmpstr, "object"))
        {
          /* Read next word */

          if (fscanf(fp, "%s", tmpstr) == EOF)
            Die(FUNCTION_NAME, "Unexpected EOF2");

          /* Check and set type */

          if (!strcmp(tmpstr, "points;"))
            *type = OF_FILE_POINTS;
          else if (!strcmp(tmpstr, "faces;"))
            *type = OF_FILE_FACES;
          else if (!strcmp(tmpstr, "owner;"))
            *type = OF_FILE_OWNER;
          else if (!strcmp(tmpstr, "neighbour;"))
            *type = OF_FILE_NEIGHBOUR;
          else if ((!strcmp(tmpstr, "rho;")) ||
                   (!strcmp(tmpstr, "rhok;")))
            *type = OF_FILE_DENSITY;
          else if (!strcmp(tmpstr, "T;"))
            *type = OF_FILE_TEMP;

          /* Break loop */

          break;
        }
    }
  while (fscanf(fp, "%s", tmpstr) != EOF);

  /* Loop until end of block */

  while ((n = fscanf(fp, "%s", tmpstr)) != EOF)
    if (!strcmp(tmpstr, "}"))
      break;

  /* Check */

  if (n == EOF)
    Die(FUNCTION_NAME, "Unexpected EOF");

  /* Read lines until size is found */

  while (1 != 2)
    {
      if (fgets(tmpstr, MAX_STR, fp) == NULL)
        Die(FUNCTION_NAME, "Unexpected EOF3");

      /* Skip comment lines */

      if ((tmpstr[0] == '/') && (tmpstr[1] == '/'))
        continue;

      /* Check dimensions array */

      if (!strncmp(tmpstr, "dimensions", 10))
        {
          /* Read dimensions */

          sscanf(tmpstr, "dimensions [%ld %ld %ld %ld %ld %ld %ld]", &dim[0],
                 &dim[1], &dim[2], &dim[3], &dim[4], &dim[5], &dim[6]);

        }
      else if (!strncmp(tmpstr, "internalField", 13))
        {
          /* Read distribution type */

          sscanf(tmpstr, "internalField %s", tmpstr2);

          /* Check uniform distribution type */

          if (!strncmp(tmpstr2, "uniform", 7))
            {
              /* An uniform distribution */

              /* Get value */

              nread = sscanf(tmpstr, "internalField %*s %lf", &unival);

              if (nread != 1)
                Die(FUNCTION_NAME, "Could not read value for unifrom internal field.\n");

              /* Store value */

              *value = unival;

              /* Store data-size */

              *sz = 1;

              /* Store distribution type */

              *disttype = OF_INTERNAL_FIELD_UNIFORM;

              /* Exit subroutine */

              return;
            }
          else if (!strncmp(tmpstr2, "nonuniform", 10))
            {
              /* A non-uniform field. Size should be on next line. */

              if (fgets(tmpstr, MAX_STR, fp) == NULL)
                Die(FUNCTION_NAME, "Unexpected EOF4");

              /* Store distribution type */

              *disttype = OF_INTERNAL_FIELD_NONUNIFORM;

              /* Convert to integer */

              n = (long)atoi(tmpstr);

              /* Check */

              if (n > 0)
                {
                  /* Put size */

                  *sz = n;

                  /* Loop until the first '(' */

                  while ((c = fgetc(fp)) != EOF)
                    if (c == '(')
                      break;

                  /* Exit subroutine */

                  return;
                }

            }
          else
            Die(FUNCTION_NAME, "Unknown type of internalField %s.", tmpstr2);
        }
      else
        {
          /* Convert to integer */

          n = (long)atoi(tmpstr);

          /* Check */

          if (n > 0)
            {
              /* Put size */

              *sz = n;

              /* Loop until the first '(' */

              while ((c = fgetc(fp)) != EOF)
                if (c == '(')
                  break;

              /* Exit subroutine */

              return;
            }
        }
    }

  /* Something wrong */

  Die(FUNCTION_NAME, "Shouldn't be here");
}

/*****************************************************************************/
