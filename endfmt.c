/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : endfmt.c                                       */
/*                                                                           */
/* Created:       2016/08/10 (JLe)                                           */
/* Last modified: 2016/08/10 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Reads MT number from formatted ENDF file                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ENDFMT:"

/*****************************************************************************/

long ENDFMT(char *line)
{
  char str[MAX_STR];

  /* Copy string and set EOF */

  strcpy(str, &line[72]);
  str[3] = '\0';

  /* Return value */

  return atol(str);
}

/*****************************************************************************/
