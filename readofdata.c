/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readofdata.c                                   */
/*                                                                           */
/* Created:       2013/12/27 (JLe)                                           */
/* Last modified: 2016/12/07 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Reads data from an OpenFOAM format file                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ReadOFData:"

/*****************************************************************************/

char *ReadOFData(FILE *fp, long type)
{
  long n;
  signed int c;
  static char line[MAX_STR];

  /* Check type */

  if (type == OF_FILE_POINTS)
    {
      /***********************************************************************/

      /***** Points file *****************************************************/

      /* Loop until the first '(' */

      while ((c = fgetc(fp)) != EOF)
        if (c == '(')
          break;

      /* Read string */

      n = 0;
      while ((c = fgetc(fp)) != EOF)
        {
          if (c == ')')
            break;
          else if (c == '\n')
            line[n++] = ' ';
          else
            line[n++] = (char)c;
        }

      /* Put EOF */

      line[n] = '\0';
      
      /***********************************************************************/
    }
  else if (type == OF_FILE_FACES)
    {
      /***********************************************************************/

      /***** Faces file ******************************************************/

      /* Loop until the first '(' */

      n = 0;
      while ((c = fgetc(fp)) != EOF)
        {
          if (c == '(')
            break;
          else if (c == '\n')
            line[n++] = ' ';
          else
            line[n++] = (char)c;
        }

      /* Add space */

      line[n++] = ' ';

      /* Read rest of the string */

      while ((c = fgetc(fp)) != EOF)
        {
          if (c == ')')
            break;
          else if (c == '\n')
            line[n++] = ' ';
          else
            line[n++] = (char)c;
        }

      /* Put EOF */

      line[n] = '\0';
      
      /***********************************************************************/
    }
  else if ((type == OF_FILE_OWNER) || (type == OF_FILE_NEIGHBOUR) ||
           (type == OF_FILE_DENSITY) || (type == OF_FILE_TEMP) ||
           (type == OF_FILE_MATERIAL) || (type == OF_FILE_MAP))
    {
      /***********************************************************************/

      /***** Owner, neighbour, density, temperature, material ****************/

      /* Read string */

      n = 0;
      while ((c = fgetc(fp)) != EOF)
        {
          if ((n > 0) && (c == '\n'))
            break;
          else if ((c != ' ') && (c != '\n'))
            line[n++] = (char)c;
        }

      /* Put EOF */

      line[n] = '\0';
      
      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid type");

  /* Return line */

  return line;
}

/*****************************************************************************/



