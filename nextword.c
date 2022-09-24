/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nextword.c                                     */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2017/05/21 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Reads next white-space separated word from string            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "NextWord:"

/*****************************************************************************/

long NextWord(char *from, char *to)
{
  long n, i;

  /* skip leading white spaces */

  n = 0;
  while(((from[n] == ' ') || (from[n] == '\n')) && (from[n] != EOF))
    n++;

  /* check if string is in quotations */
  
  if (from[n] == '\"')
    {
      n++;

      /* read string */

      i = 0;

      while ((from[n] != '\"') && (from[n] != '\0'))
        {
          if (i < MAX_STR)
            to[i++] = from[n++];
          else
            Die(FUNCTION_NAME, "String length exceeded");
               }
      
      n++;
    }
  else
    {
      /* read string */

      i = 0;

      while((from[n] != ' ') && (from[n] != '\n') && (from[n] != '\0'))
        {
          if (i < MAX_STR)
            to[i++] = from[n++];
          else
            Die(FUNCTION_NAME, "String length exceeded");
               }
    }

  /* terminate string */

  if (i < MAX_STR)
    to[i] = '\0';
  else
    Die(FUNCTION_NAME, "String length exceeded");

  return n;
}

/*****************************************************************************/
