/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkvalue.c                                   */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Checks that value is within allowed limits                   */
/*                                                                           */
/* Comments: - Adopted from Serpent 1.1.13                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "CheckValue:"

/*****************************************************************************/

#ifdef DEBUG 

void CheckValue(char *func, char *param, char *comment, double value, 
	       double min, double max)
{
  char msg[MAX_STR];


  if ((value >= min) && (value <= max))
    return;
  else if (value > max)
    sprintf(msg, "Value %E of parameter \"%s\" above upper limit %E %s", value, 
	param, max, comment);
  else if (value < min)
    sprintf(msg, "Value %E of parameter \"%s\" below lower limit %E %s", value, 
	param, min, comment);
  else
    sprintf(msg, "Invalid value %E of parameter \"%s\" %s", value, param, 
	    comment);

  Die(func, "%s\n\nNOTE: This value check was performed because the code was compiled in the\n      debug mode and it may or may not be an indication of an actual problem.\n      The debugger mode can be switched off by recompiling the source code\n      without the -DDEBUG option (see Makefile).", msg);
}

#endif

/*****************************************************************************/
