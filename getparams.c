/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getparams.c                                    */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2017/04/03 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Reads parameters for an input card from file                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "input_params.h"

#define FUNCTION_NAME "GetParams:"

/*****************************************************************************/

char **GetParams(char *w0, char *input, long *np0, long *i0, long min,
               long max, char *file)
{
  long np, n, i, stop, line, i1;
  char word[MAX_STR], **params;

  /* Get line number for error message */

  line = GetLineNumber(input, *i0);

  /* Reset stop flag */

  stop = NO;

  /* Loop over words and count*/

  np = 0;
  i1 = *i0;

  while ((i = NextWord(&input[i1], word)) > 0)
    {
      /* Check that word is not an input parameter */

      n = 0;
      while(inputwords[n][0] != '\0')
        if (!strcasecmp(inputwords[n++], word))
          stop = YES;

      if (stop == YES)
        break;

      /* Check overflow of params array */

      if (np > MAX_INPUT_PARAMS - 1)
        Die(FUNCTION_NAME, "params[] array size exceeded");

      /* Add count */

      np++;

      /* Last param can be of zero length */

      if (strlen(word) == 0)
        np--;

      /* Add to counter */

      i1 = i1 + i;
    }

  /* Allocate memory */

  params = (char **)Mem(MEM_ALLOC, np + 1, sizeof(char *));

  for(n = 0; n < np + 1; n++)
    params[n] = (char *)Mem(MEM_ALLOC, MAX_STR, sizeof(char));

  /* Reset stop flag */

  stop = NO;

  /* Loop over words */

  np = 0;

  while ((i = NextWord(&input[*i0], word)) > 0)
    {
      /* Check that word is not an input parameter */

      n = 0;
      while(inputwords[n][0] != '\0')
        if (!strcasecmp(inputwords[n++], word))
          stop = YES;

      if (stop == YES)
        break;

      /* Check overflow of params array */

      if (np > MAX_INPUT_PARAMS - 1)
        Die(FUNCTION_NAME, "params[] array size exceeded");

      /* Copy parameter */

      strcpy(params[np], word);

      np++;

      /* Last param can be of zero length */

      if (strlen(word) == 0)
        np--;

      /* Add to counter */

      *i0 = *i0 + i;
    }

  /* Check for reserved parameter name */

  if ((np == 0) && (np < min))
    {
      if (word[0] == '\0')
        Error(-1, w0, file, line, "Missing parameters after command word");
      else
        Error(-1, w0, file, line, "Value \"%s\" is a reserved command word",
              word);
    }

  /* Check number of parameters */

  if ((np < min) || (np > max))
    {
      if (min == max)
        Error(-1, w0, file, line, "Number of parameters must be %d", min);
      else if ((np == 0) && (!strcasecmp(w0, "set")) && (strlen(word) == 0))
        Error(-1, w0, file, line, "Missing parameter name");
      else if ((np == 0) && (!strcasecmp(w0, "set")) && (strlen(word) > 0))
        Error(-1, w0, file, line,
              "Parameter %s must be entered without \"set\"", word);
      else if (np < min)
        Error(-1, w0, file, line, "Not enough parameters given (min = %d)",
              min);
      else
        Error(-1, w0, file, line, "Too many parameters given (max = %d)", max);
    }

  /* Put number of parameters */

  *np0 = np;

  /* Return pointer to list */

  return params;
}

/*****************************************************************************/
