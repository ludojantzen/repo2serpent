/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getlinenumber.c                                */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2015/07/17 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Line number from position                                    */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetLineNumber:"

/*****************************************************************************/

long GetLineNumber(char *file, long imax)
{
  long nl, i;

  /* Reset number of lines */

  nl = (long)RDB[DATA_LINE_NUM_NL0];

  /* Loop to position */

  for (i = (long)RDB[DATA_LINE_NUM_N0]; i < imax; i++)
    {
      /* Check eof */

      if (file[i] == '\0')
        Die(FUNCTION_NAME, "End-of-file encountered");

      /* Check newline and add counter */

      else if (file[i] == '\n')
        nl++;
    }

  /* Update counters */

  WDB[DATA_LINE_NUM_NL0] = (double)nl;
  WDB[DATA_LINE_NUM_N0] = (double)imax;

  /* Return line number */

  return nl;
}

/*****************************************************************************/
