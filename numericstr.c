/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : numericstr.c                                   */
/*                                                                           */
/* Created:       2017/09/07 (JLe)                                           */
/* Last modified: 2018/06/24 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Converts array to a numeric value for quick matching         */
/*                                                                           */
/* Comments: - Currently used for linking surfaces to cell surface lists     */
/*             and works only for integers.                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "NumericStr:"

/*****************************************************************************/

long NumericStr(char *str)
{
  long n;

  /* Check that all characters are digits */

  n = 0;
  while (str[n] != '\0')
    if (!isdigit(str[n++]))
      return -1;

  /* Check that first digit is not zero */

  if (str[0] == '0')
    return -1;

  /* Return number */

  return atol(str);
}

/*****************************************************************************/
