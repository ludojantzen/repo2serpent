/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : timestr.c                                      */
/*                                                                           */
/* Created:       2010/11/14 (JLe)                                           */
/* Last modified: 2018/07/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Converts number of seconds to human-readable time string.    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "TimeStr:"

/*****************************************************************************/

char *TimeStr(long t)
{
  long h, min, s;
  char ss[10], mm[10];
  static char str[20];

  h = (long)((double)t/3600.0);
  t = t - h*3600;
  min = (long)((double)t/60.0);
  t = t - min*60;
  s = (long)t;

  /* Check negative */

  if (t < 0.0)
    {
      /* Not valid */

      sprintf(str, "-:--:--");

      /* Exit subroutine */

      return str;
    }

  /* Get seconds */

  if (s < 10)
    sprintf(ss, "0%ld", s);
  else
    sprintf(ss, "%ld", s);

  /* Get minutes */

  if (min < 10)
    sprintf(mm, "0%ld", min);
  else
    sprintf(mm, "%ld", min);
  
  /* Compose string */

  if (h < 1000000)
    sprintf(str, "%ld:%s:%s", h, mm, ss);
  else
    sprintf(str, "N/A");

  /* Return string */

  return str;
}

/*****************************************************************************/
