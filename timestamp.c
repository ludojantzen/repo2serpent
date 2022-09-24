/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : timestamp.c                                    */
/*                                                                           */
/* Created:       2012/05/29 (JLe)                                           */
/* Last modified: 2012/05/29 (JLe)                                           */
/* Version:       2.1.6                                                      */
/*                                                                           */
/* Description: Returns current date and time in a string                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "TimeStamp:"

/*****************************************************************************/

char *TimeStamp()
{
  static char tmpstr[MAX_STR];
  time_t t0;

  /* Get time */

  time(&t0);

  /* Convert to string */

  strcpy(tmpstr, ctime(&t0));
  tmpstr[(int)strlen(tmpstr) - 1] = '\0';
 
  /* Return */

  return tmpstr;
}

/*****************************************************************************/
