/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : testdosfile.c                                  */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2016/12/06 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Checks that an input file is not in DOS format               */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "TestDOSFile:"

/*****************************************************************************/

void TestDOSFile(char *fname)
{
  long i;
  signed int c, c0;
  FILE *fp;

  /* Open file */

  if ((fp = fopen(fname, "r")) == NULL)
    return;

  /* Reset length */

  i = 0;

  /* Reset last character */

  c0 = ' ';

  /* Loop over file */

  while((c = fgetc(fp)) != EOF)
    {
      /* Check ascii format */
      
      if ((!isascii((long)c)) && (c != -61)  && (c != -92) && (c != -124) &&
          (c != -74) && (c != -106) && (c != -91) && (c != -123) && 
          (c != -68) && (c != -100))
        {
          fprintf(errp, "\nFile:\n\n\"%s\"\n\nis not in ASCII format.\n\n",
                  fname);
          exit(-1);
        }
      
      /* Check newline */

      if (c == '\n')
        {
          /* Check dos format */

          if (c0 == '\r')
            {
              Error(0, "File:\n\n\"%s\"\n\nmust be converted from DOS to UNIX format", fname);
            }
          else
            {
              /* Close file and exit */

              fclose(fp);

              return;
            }
        }

      /* Remember last character */

      c0 = c;

      /* Add to length */

      i++;
    }

  /* Something wrong here */

  if (i == 0)
    Error(0, "File:\n\n\"%s\"\n\n has zero length",
          fname);

  /* Close file */

  fclose(fp);
}

/*****************************************************************************/
