/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : countsensparams.c                              */
/*                                                                           */
/* Created:       2017/05/22 (VVa)                                           */
/* Last modified: 2017/05/22 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Counts the number of non-control-word parameters before the  */
/*              next control word.                                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "CountSensParams:"

/*****************************************************************************/

long CountSensParams(char **ptr, long maxi)
{
  long n, i, stop;
  char *inputwords[] = {
    "matlist",
    "zailist",
    "mtlist",
    "realist",
    "efunc",
    "xyzfunc",
    "\0"
  };

  stop = NO;

  for (i = 0; i < maxi; i++)
    {
      /* Check that word is not a control word */
      n = 0;
      while(inputwords[n][0] != '\0')
        if (!strcasecmp(inputwords[n++], ptr[i]))
          stop = YES;


      /* If a control word was found, break */

      if (stop == YES)
        break;
    }

  /* Return number of parameters */

  return i;
}

/*****************************************************************************/
