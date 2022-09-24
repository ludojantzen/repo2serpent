/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : testasciifile.c                                */
/*                                                                           */
/* Created:       2014/03/07 (JLe)                                           */
/* Last modified: 2014/12/09 (JLe)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description: Checks that an input file is in ASCII format                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "TestASCIIFIle:"

/*****************************************************************************/

long TestASCIIFile(char *fname)
{
  char c;
  FILE *fp;

  /* Open file */

  if ((fp = fopen(fname, "r")) == NULL)
    return -1;

  /* Loop over file */

  while (fread(&c, sizeof(char), 1, fp))
    {
      /* Check ascii format */

      if ((!isascii((long)c)) && (c != -61)  && (c != -92) && (c != -124) &&
          (c != -74) && (c != -106) && (c != -91) && (c != -123) && 
          (c != -68) && (c != -100))
        {
          /* Close file */

          fclose(fp);

          /* Exit subroutine */

          return NO;
        }
    }

  /* Close file */

  fclose(fp);

  /* Exit subroutine */

  return YES;
}

/*****************************************************************************/
