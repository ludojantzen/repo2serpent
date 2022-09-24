/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : endfcolf.c                                     */
/*                                                                           */
/* Created:       2010/11/20 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Reads integer from column formatted ENDF file                */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ENDFColI:"

/*****************************************************************************/

long ENDFColI(long col, char *line)
{
  char str[MAX_STR];

  /* Check column */

  if ((col < 1) || (col > 6))
    Die(FUNCTION_NAME, "Invalid column index %d", col);

  /* Copy string and set EOF */

  strcpy(str, &line[11*(col - 1)]);
  str[11] = '\0';

  /* Return value */

  return atoi(str);
}

/*****************************************************************************/
