/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : gettext.c                                      */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Retrieves pointer to ASCII data block                        */
/*                                                                           */
/* Comments: - Adopted from Serpent 1.0.0                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetText:"

/*****************************************************************************/

char *GetText(long ptr)
{
  /* Check pointer */

  if (ptr < 1)
    Die(FUNCTION_NAME, "Invalid pointer");

  /* Get pointer to data array */

  ptr = (long)RDB[ptr];

  /* Check pointer */

  if ((ptr < VALID_PTR) || (ptr > (long)RDB[DATA_ASCII_DATA_SIZE]))
    Die(FUNCTION_NAME, "Invalid pointer");

  /* Return pointer to string */

  return &ASCII[ptr];
}

/*****************************************************************************/
