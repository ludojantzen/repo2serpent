/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : atoi.c                                         */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Array to long conversion with type checking                  */
/*                                                                           */
/* Comments: - From Serpent 1.1.8                                            */
/*           - Tätä pitäisi kutsua pelkästään GetParam():sta.                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "AtoI:"

/*****************************************************************************/

long AtoI(char *str, char *param, char *file, long line)
{
  long n;

  for(n = 0; n < (long)strlen(str); n++)
    if (!isdigit(str[n]) && (str[n] != '-'))
      {
	/* Check mode */

	if (line > 0)
	  {
	    /* Error in input file */

	    Error(-1, param, file, line, "Invalid integer entry \"%s\"",str);
	  }
	else
	  {
	    /* Error in string */

	    Die(FUNCTION_NAME, "Invalid integer entry \"%s\"", str);
	  }
      }

  /* Return number */

  return atoi(str);
}

/*****************************************************************************/
