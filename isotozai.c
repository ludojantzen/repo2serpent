/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : isotozai.c                                     */
/*                                                                           */
/* Created:       2010/09/23 (JLe)                                           */
/* Last modified: 2017/02/27 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: - Converts isotope symbol to ZAI                             */
/*                                                                           */
/* Comments: - Taken from Serpent 1.1.10 with small modifications.           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "element_data.h"

#define FUNCTION_NAME "IsotoZAI:"

/*****************************************************************************/

long IsotoZAI(char *orig)
{     
  long n, i, Z, A, I;
  char sym[MAX_STR], iso[MAX_STR];

  /* Copy original string */

  strcpy(iso, orig);

  /* Try element name and symbol */

  n = 1;
  
  while(*element_symbols[n] != '\0')
    {
      if (!strcasecmp(element_names[n], iso))
        return n;
      else if (!strcasecmp(element_symbols[n], iso))
        return n;
      
      n++;
    }

  /* Separate element symbol */

  for (n = 0; n < (long)strlen(iso); n++)
    {
      /* Break string */

      if ((iso[n] == '-') || (isdigit(iso[n])))
        break;
      
      /* Copy character */
      
      sym[n] = iso[n];
    }

  /* Terminate string */

  sym[n] = '\0';

  /* Remove trailing '-' */

  if (iso[n] == '-')
    n++;

  /* Get I */

  if ((iso[strlen(iso) - 1] =='m') || (iso[strlen(iso) - 1] =='M'))
    {
      /* First isomeric state */

      I = 1;

      /* Remove trailing m for next check */

      iso[strlen(iso) - 1] = '\0';
    }
  else if ((iso[strlen(iso) - 1] =='n') || (iso[strlen(iso) - 1] =='N'))
    {
      /* Second isomeric state */

      I = 2;

      /* Remove trailing n for next check */

      iso[strlen(iso) - 1] = '\0';
    }
  else
    {
      /* Ground state */

      I = 0;
    }

  /* Check that remaining part is numerical */

  for (i = n; i < (long)strlen(iso); i++)
    if (!isdigit(iso[i]))
      return -1;

  /* Get A */

  A = atoi(&iso[n]);

  /* Check range */

  if ((A < 1) || (A > 299))
    return -1;

  /* Get Z (zero is reserved for neutron) */

  for (Z = 1; Z < NUMBER_OF_ELEMENTS; Z++)
    if (!strcasecmp(element_symbols[Z], sym))
      break;

  /* Check error */

  if (Z == NUMBER_OF_ELEMENTS)
    return -1;

  /* Return ZAI */
  
  return 10000*Z + 10*A + I;
}

/*****************************************************************************/
