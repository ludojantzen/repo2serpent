/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : endfmf.c                                       */
/*                                                                           */
/* Created:       2016/08/10 (JLe)                                           */
/* Last modified: 2016/08/10 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Reads MF number from formatted ENDF file                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ENDFMF:"

/*****************************************************************************/

long ENDFMF(char *line)
{
  char str[MAX_STR];

  /* Copy string and set EOF */

  strcpy(str, &line[70]);
  str[2] = '\0';

  /* Return value */

  return atol(str);
}

/*****************************************************************************/
