/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : endfmt.c                                       */
/*                                                                           */
/* Created:       2016/08/12 (JLe)                                           */
/* Last modified: 2016/08/12 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Reads material number from formatted ENDF file               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ENDFMat:"

/*****************************************************************************/

long ENDFMat(char *line)
{
  char str[MAX_STR];

  /* Copy string and set EOF */

  strcpy(str, &line[66]);
  str[4] = '\0';

  /* Return value */

  return atol(str);
}

/*****************************************************************************/
