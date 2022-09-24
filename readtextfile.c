/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readtextfile.c                                 */
/*                                                                           */
/* Created:       2010/10/11 (JLe)                                           */
/* Last modified: 2018/08/30 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads text file and removes unnecessary stuff                */
/*                                                                           */
/* Comments: - Converted from Serpent 1.1.14 with a few additions            */
/*           - # no longer works as comment character (2.1.25)               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadTextFile:"

/*****************************************************************************/

char *ReadTextFile(char *file)
{
  long sz, n;
  char *txt, tmpstr[MAX_STR];
  signed int c, c0;
  FILE *fp, *outfp;

  /* Test file format */

  TestDOSFile(file);

  /* Avoid compiler warning */

  outfp = NULL;

  /* open file */

  if ((fp = fopen(file, "r")) != NULL)
    {
      /* Open output file */

      if ((long)RDB[DATA_COPY_INPUTS] == YES)
        {
          sprintf(tmpstr, "%s.input", GetText(DATA_PTR_INPUT_FNAME));
          
          if ((outfp = fopen(tmpstr, "a")) == NULL)
            Die(FUNCTION_NAME, "Unable to open file for writing");
        }

      /* get file size */

      sz = 0;
      while((c = fgetc(fp)) != EOF)
        {
          if ((c == '#') || (c == ':') || (c == '(') || (c == ')'))
            sz = sz + 3;
          else if (c == '-')
            sz = sz + 2;
          else
            sz++;

          /* Write to output */

          if ((long)RDB[DATA_COPY_INPUTS] == YES)
            fputc(c, outfp);
        }

      /* Check if inputs are copied */

      if ((long)RDB[DATA_COPY_INPUTS] == YES)
        {
          /* Print newline */
          
          fprintf(outfp, "\n%% --- End of file ---\n\n");
          
          /* Close output file */

          fclose(outfp);        
        }

      /* Allocate memory */

      txt = (char *)Mem(MEM_ALLOC, sz + 1, sizeof(char));
        
      /* Rewind to beginning of file */

      fseek(fp, 0, 0);

      /* read file one character at a time */

      c0 = ' ';
      n = 0;

      while((c = fgetc(fp)) != EOF)
        {
          /***** cleanup *****************************************************/

          /* Leading space for minus sign (for cell surface lists) */

          /* Tää ei sit toimi symbolisten nuklidinimien kanssa */

          /*
          if (c == '-')
            {
              if ((n > 0) && (txt[n - 1] != 'e') && (txt[n - 1] != 'E'))
                txt[n++] = ' ';

              txt[n++] =  c;
            }
          */

          /* alphanumeric characters etc. copied as they are. */
          
          /*else*/ if (isalnum(c) || (c == '.') || (c == '-') || (c == '+') || 
                       (c == '_') || (c == '\n'))
            txt[n++] = (char)c;

          /* Special characters for cell surface lists */

          else if ((c == '#') || (c == ':') || (c == '(') || (c == ')'))
            {
              txt[n++] = ' ';
              txt[n++] = (char)c;
              txt[n++] = ' ';
            }

          /* comment lines */

          else if (c == '%')
            {
              while(((c = fgetc(fp)) != EOF) && (c != '\n'))
                txt[n++] = ' ';

              if (c != EOF)
                txt[n++] = '\n';
            }

          /* comment sections */

          else if ((c0 == '/') && (c == '*'))
            {
              while((c = fgetc(fp)) != EOF)
                {
                  if ((c0 == '*') && (c == '/'))
                    break;

                  if (c0 == '\n')
                    txt[n++] = '\n';
                  else
                    txt[n++] = ' ';
                  
                  c0 = c;
                }
            }

          /* strings in quotation marks */

          else if (c == '\"')
            {
              txt[n++] = '\"';
              while(((c = fgetc(fp)) != EOF) && (c != '\"'))
                txt[n++] = (char)c;

              if (c == EOF)
                {
                  fprintf(errp, "%s Unbalanced quotation marks.\n\n",
                          FUNCTION_NAME);
                  exit(-1);
                }
              else
                txt[n++] = '\"';
            }

          else

            /* rest of the characters are converted to space */

            txt[n++] = ' ';

          /* Remember last character */

          c0 = c;

          /*******************************************************************/
        }

      /* Check count */

      if (n > sz)
        Die(FUNCTION_NAME, "Size error (%ld > %ld)", n, sz);
      
      /* close file */

      fclose(fp);
    }
  else
    {
      /* failed to open file */

      fprintf(errp, "%s Unable to open file \"%s\n", FUNCTION_NAME, file);
      exit(-1);
    }

  /* exit subroutine */

  return txt;
}

/*****************************************************************************/
