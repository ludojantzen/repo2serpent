/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : idxstr.c                                       */
/*                                                                           */
/* Created:       2014/02/13 (JLe)                                           */
/* Last modified: 2014/02/13 (JLe)                                           */
/* Version:       2.1.17                                                     */
/*                                                                           */
/* Description: Converts index to string with leading zeros.                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "IdxStr:"

/*****************************************************************************/

char *IdxStr(long idx, long max)
{
  long n;
  char fmt[MAX_STR];
  static char str[20];
  
  /* Compare value to maximum */

  if (idx > max)
    max = idx;

  /* Calculate number of digits */

  sprintf(str, "%ld", max);
  n = strlen(str);

  /* Formatted string */

  sprintf(fmt, "%%0%ld.0f", n);
  sprintf(str, fmt, (double)idx);

  /* Return string */

  return str;
}

/*****************************************************************************/
