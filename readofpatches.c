/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readofbatches.c                                */
/*                                                                           */
/* Created:       2014/03/30 (JLe)                                           */
/* Last modified: 2019/03/27 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads patch data from OpenFOAM format files                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadOFPatches:"

/*****************************************************************************/

void ReadOFPatches(long loc0)
{
  long n, loc1;
  char tmpstr[MAX_STR];
  signed int c;
  FILE *fp;

  /* Loop over files */

  for (n = 0; n < 6; n++)
    {
      /* Get file name */

      if (n == 0)
        sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_PFILE));
      else if (n == 1)
        sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_FFILE));
      else if (n == 2)
        sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_OFILE));
      else if (n == 3)
        sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_NFILE));
      else if (n == 4)
        {
          /* UMSH geometries do not supply a density file */

          if ((long)RDB[loc0 + IFC_PTR_OF_RFILE] < VALID_PTR)
            continue;

          sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_RFILE));
        }
      else if (n == 5)
        {
          /* UMSH geometries do not supply a temperature file */

          if ((long)RDB[loc0 + IFC_PTR_OF_TFILE] < VALID_PTR)
            continue;

          sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_TFILE));
        }
      else
        Die(FUNCTION_NAME, "Overflow");

      /* Open file for reading */

      if ((fp = fopen(tmpstr, "r")) == NULL)
        continue;

      /* Look over file and look for key word "boundaryField" */

      while (fscanf(fp, "%s", tmpstr) != EOF)
        {
          /* Check string */

          if (!strcmp(tmpstr, "boundaryField"))
            {
              /* Loop until begin marker */

              while ((c = fgetc(fp)) != EOF)
                if (c == '{')
                  break;

              /* Loop */

              while (1 != 2)
                {
                  /* Read next word */

                  if (fscanf(fp, "%s", tmpstr) == EOF)
                    break;

                  /* Check end marker */

                  if (tmpstr[0] == '}')
                    break;
                  else
                    {
                      /* Find match */

                      loc1 = (long)RDB[loc0 + IFC_PTR_OF_PATCHES];
                      while (loc1 > VALID_PTR)
                        {
                          /* Compare names */

                          if (!strcmp(GetText(loc1 + IFC_OF_PATCH_PTR_NAME),
                                      tmpstr))
                            break;

                          /* Next */

                          loc1 = NextItem(loc1);
                        }

                      /* Check pointer */

                      if (loc1 < VALID_PTR)
                        {
                          /* Add new batch */

                          loc1 = NewItem(loc0 + IFC_PTR_OF_PATCHES,
                                         IFC_OF_PATCH_BLOCK_SIZE);

                          /* Put name */

                          WDB[loc1 + IFC_OF_PATCH_PTR_NAME] =
                            (double)PutText(tmpstr);
                        }
                      /* Find type and save it */

                      while (fscanf(fp, "%s", tmpstr) != EOF)
                        {
                          /* Check for the control word "type" */

                          if (!strcmp(tmpstr, "type"))
                            {
                              /* Try to read type */

                              if (fscanf(fp, "%s", tmpstr) != EOF)
                                {
                                  tmpstr[strlen(tmpstr)-1] = 0;
                                  WDB[loc1 + IFC_OF_PATCH_PTR_TYPE] = PutText(tmpstr);
                                }
                              else
                                Error(loc1, "File ended while trying to read type for "
                                      "OpenFOAM patch %s", GetText(loc1 + IFC_OF_PATCH_PTR_NAME));

                              /* Break after reading and storing type */

                              break;
                            }
                        }

                      /* Loop until end marker */

                      while ((c = fgetc(fp)) != EOF)
                        if (c == '}')
                          break;
                    }
                }
            }
        }
    }
}

/*****************************************************************************/
