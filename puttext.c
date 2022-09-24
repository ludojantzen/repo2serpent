/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : puttext.c                                      */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2014/02/28 (JLe)                                           */
/* Version:       2.1.19                                                     */
/*                                                                           */
/* Description: Stores text string in ASCII data block                       */
/*                                                                           */
/* Comments: - Adopted from Serpent 1.0.0                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PutText:"

/*****************************************************************************/

long PutText(char *str)
{
  long ptr, memsize, n, nl;
  char *tmpstr;

  /* Avoid compiler warning */

  nl = -1;
  tmpstr = NULL;

  /* Check string */

  if (str == NULL)
    Die(FUNCTION_NAME, "Null pointer");
  else if (*str == EOF)
    Die(FUNCTION_NAME, "Empty string");
  else if ((nl = strlen(str)) < 1)
    Die(FUNCTION_NAME, "Empty string");

  /* Get pointer to free data */

  memsize = (long)RDB[DATA_ASCII_DATA_SIZE];
  ptr = memsize;
  
  /* Copy str to tmpstr if str points to ASCII array. */
  /* Otherwise reallocation can destroy str. */
  
  if ((str >= ASCII) && (str < ASCII + ptr)) 
    {
      /* Allocate memory for a temporary array */

      tmpstr = (char *)Mem(MEM_ALLOC, nl + 1, sizeof(char));
      strcpy(tmpstr, str);
      str = tmpstr;
    }

  /* Allocate memory for data */

  memsize = memsize + nl + 1;
  ASCII = (char *)Mem(MEM_REALLOC, ASCII, memsize*sizeof(char));

  /* Store data */

  for (n = 0; n < nl; n++)
    ASCII[ptr + n] = str[n];

  /* Terminate string */

  ASCII[ptr + n] = '\0';

  /* Update data block size */

  WDB[DATA_ASCII_DATA_SIZE] = (double)memsize;

  /* Free temporary array */

  if (tmpstr != NULL)
    Mem(MEM_FREE, tmpstr);

  /* Calculate allocated memory in bytes */

  CalculateBytes();

  /* Return pointer */
  
  return ptr;
}

/*****************************************************************************/
