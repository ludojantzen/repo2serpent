/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : timeintervalstr.c                              */
/*                                                                           */
/* Created:       2010/09/11 (JLe)                                           */
/* Last modified: 2017/12/11 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Converts number of seconds to human-readable time interval   */
/*              string.                                                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "TimeIntervalStr;"

/*****************************************************************************/

char *TimeIntervalStr(double sec)
{                
  static char str[MAX_STR];
  long n;
  double rem;
  char unit[MAX_STR];

  /* Check interval */

  if (sec == 0.0)
    sprintf(unit, "seconds");
  else if (sec < 0.1E-6)
    {
      sec = sec*1E+9;
      sprintf(unit, "nanoseconds");
    }
  else if (sec < 0.1E-3)
    {
      sec = sec*1E+6;
      sprintf(unit, "microseconds");
    }
  else if (sec < 0.1)
    {
      sec = sec*1E+3;
      sprintf(unit, "milliseconds");
    }
  else if (sec < 60.0)
    sprintf(unit, "seconds");
  else if (sec < 60.0*60.0)
    {
      sec = sec/60.0;
      sprintf(unit, "minutes");
    }
  else if (sec < 24*60.0*60.0)
    {
      sec = sec/(60.0*60.0);
      sprintf(unit, "hours");
    }
  else if (sec < 365.0*24*60.0*60.0)
    {
      sec = sec/(24.0*60.0*60.0);
      sprintf(unit, "days");
    }
  else if (sec < 1E+6*365.0*24*60.0*60.0)
    {
      sec = sec/(365.0*24.0*60.0*60.0);
      sprintf(unit, "years");
    }
  else if (sec < 1E+9*365.0*24*60.0*60.0)
    {
      sec = sec/(1E+6*365.0*24.0*60.0*60.0);
      sprintf(unit, "million years");
    }
  else if (sec < 1E+12*365.0*24*60.0*60.0)
    {
      sec = sec/(1E+9*365.0*24.0*60.0*60.0);
      sprintf(unit, "billion years");
    }
  else
    {
      sec = sec/(365.0*24.0*60.0*60.0);
      sprintf(unit, "years");
    }

  /* Construct string */

  if (sec < 100.0)
    sprintf(str, "%1.3G %s", sec, unit);
  else if (sec < 1E+3)
    sprintf(str, "%1.0f %s", sec, unit);
  else if (sec < 1E+6)
    {
      /* Separate thousands */

      n = (long)(sec/1000.0);
      rem = sec - (double)n*1000.0;

      if (rem < 10.0)
        sprintf(str, "%ld,00%1.0f %s", n, rem, unit);
      else if (rem < 100.0)
        sprintf(str, "%ld,0%1.0f %s", n, rem, unit);
      else
        sprintf(str, "%ld,%1.0f %s", n, rem, unit);
    }
  else
    sprintf(str, "%1.3G %s", sec, unit);

  /* Return string */

  return str;
}

/*****************************************************************************/
