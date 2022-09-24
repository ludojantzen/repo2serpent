/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : endfnewline.c                                  */
/*                                                                           */
/* Created:       2012/12/28 (JLe)                                           */
/* Last modified: 2012/12/28 (JLe)                                           */
/* Version:       2.1.12                                                     */
/*                                                                           */
/* Description: Reads next line from a column formatted ENDF file            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ENDFNewLine:"

/*****************************************************************************/

void ENDFNewLine(char *func, char *line, FILE *fp)
{
  /* Check file pointer and line */

  if (fp == NULL)
    Die(FUNCTION_NAME, "Error in file pointer");
  else if (line == NULL)
    Die(FUNCTION_NAME, "Error in line pointer");

  /* Read */
  
  if (fgets(line, 82, fp) == NULL)
    Die(func, "fgets error");
}

/*****************************************************************************/
