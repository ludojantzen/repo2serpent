/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : Error.c                                        */
/*                                                                           */
/* Created:       2010/09/23 (JLe)                                           */
/* Last modified: 2019/11/12 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Prints input error message                                   */
/*                                                                           */
/* Comments: - Converted from Serpent 1.1.12                                 */
/*           - Use this for user errors, errors in code terminate the run    */
/*             with Die()                                                    */
/*           - This and Due() are the only functions that should make        */
/*             a call to exit().                                             */
/*           - Output changed from outp to errp (outp for MPI tasks > 0 is   */
/*             set to dev/null)                                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "input_params.h"

#define FUNCTION_NAME "Error:"

/*****************************************************************************/

void Error(long ptr, ...)
{
  char param[MAX_STR], fname[MAX_STR], *format;
  long line, n;
  va_list argp;
  va_start (argp, ptr);

  /* Initialize variables to avoid compiler warning */

  *param = '\0';
  *fname = '\0';
  line = 0;

  /* Check pointer */

  if (ptr > VALID_PTR)
    {
      /* Put parameter name */

      if ((long)RDB[ptr + PARAM_PTR_NAME] > VALID_PTR)
        strcpy(param, GetText(ptr + PARAM_PTR_NAME));
      else
        ptr = 0;

      /* Check pointer */

      if (ptr > VALID_PTR)
        {
          /* Put file name */

          if ((long)RDB[ptr + PARAM_PTR_FNAME] > VALID_PTR)
            strcpy(fname, GetText(ptr + PARAM_PTR_FNAME));
          else
            Die(FUNCTION_NAME, "File name not set");

          /* Put line number */

          line = (long)RDB[ptr + PARAM_LINE];
        }
    }
  else if (ptr < 0)
    {
      /* Use argument values */

      strcpy(param, va_arg(argp, char *));
      strcpy(fname, va_arg(argp, char *));
      line = va_arg(argp, int);
    }

  /* Print error message */

  fprintf(errp, "\n***** %s\n\n", TimeStamp());

  if (ptr != 0)
    {
      if (line > -1)
        fprintf(errp, "Input error in parameter \"%s\" on line %ld ", param,
                line);
      else
        fprintf(errp, "Input error in parameter \"%s\" ", param);

      fprintf(errp, "in file \"%s\":\n\n", fname);
    }
  else
    fprintf(errp, "Input error:\n\n");

  format = va_arg(argp, char *);
  vfprintf(errp, format, argp);

  fprintf(errp, "\n\n");

  /* Don't print link for duplicate definitions or unidentified options */

  if (!strncmp(format, "Duplicate definition", 20) ||
      !strncmp(format, "Unidentified option", 19))
    ptr = 0;

  /* Print link to Serpent wiki */

  if (ptr != 0)
    {
      /* Find parameter from list of input parameters */

      n = 0;
      while (*inputwords[n] != '\0')
        {
          if (!strcasecmp(inputwords[n], param))
            break;
          else
            n++;
        }

      /* Print */

      if (strcmp(param, "set"))
        {
          fprintf(errp, "For parameter description, try Serpent Wiki:\n\n");
          fprintf(errp,
                  "http://serpent.vtt.fi/mediawiki/index.php/Input_syntax_manual#");

          if (*inputwords[n] == '\0')
            fprintf(errp, "set_%s\n\n", param);
          else
            fprintf(errp, "%s\n\n", param);
        }
    }

  /* Exit subroutine */

  exit(-1);
}

/*****************************************************************************/
