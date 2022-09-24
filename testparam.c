/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : testparam.c                                    */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Reads input value with type and boundary checking            */
/*                                                                           */
/* Comments: - From Serpent 1.1.11                                           */
/*           - Tämän pitäisi korvata kokonaan AtoI() ja AtoF().              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "TestParam:"

/*****************************************************************************/

double TestParam(char *pname, char *fname, long line, char *val, 
		 long type, ...)
{
  double retval, min, max;
  va_list argp;  
  va_start (argp, type);
    
  /* Avoid compiler warning */

  retval = -1.0;

  /* Check type */

  if (type == PTYPE_LOGICAL)
    {
      /* Check given value */

      if (!strcasecmp(val, "yes"))
	retval = (double)YES;
      else if (!strcasecmp(val, "no"))
	retval = (double)NO;
      else if (!strcasecmp(val, "1"))
	retval = (double)YES;
      else if  (!strcasecmp(val, "0"))
	retval = (double)NO;
      else
	Error(-1, pname, fname, line, 
	      "Invalid value \"%s\", must be 1, 0, \"yes\" or \"no\"", val);
    }
  else if (type == PTYPE_REAL)
    {
      /* Get value */

      retval = AtoF(val, pname, fname, line);

      /* Get boundaries */

      min = va_arg(argp, double);
      max = va_arg(argp, double);

      /* Check boundaries */

      if ((retval < min) || (retval > max))
	Error(-1, pname, fname, line,
	      "Invalid value \"%s\", must be between %1.3E and %1.3E", val,
	      min, max);
    }
  else if (type == PTYPE_INT)
    {
      /* Get value */

      retval = (double)AtoI(val, pname, fname, line);

      /* Get boundaries (NOTE: toi bugittaa longilla!!!) */

      min = (double)va_arg(argp, int);
      max = (double)va_arg(argp, int);

      /* Check boundaries */

      if ((retval < min) || (retval > max))
	Error(-1, pname, fname, line,
	      "Invalid value \"%s\", must be between %ld and %ld", val,
	      (long)min, (long)max);
    }
  else
    Die(FUNCTION_NAME, "Invalid type %d", type);

  /* Return value */

  return retval;
}

/*****************************************************************************/
