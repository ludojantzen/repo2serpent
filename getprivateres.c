/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getpricateres.c                                */
/*                                                                           */
/* Created:       2011/11/15 (JLe)                                           */
/* Last modified: 2012/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Returns reduced value from RES2 block                        */
/*                                                                           */
/* Comments: - Replaced by macro when not compiled in debug mode             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetPrivateRes:"

/*****************************************************************************/

#ifdef DEBUG

double GetPrivateRes(long ptr)
{
  double val;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "RES2 array not ready for access");
  
  /* Check if buffer is reduced */

  if ((long)RDB[DATA_RES2_REDUCED] == NO)
    Die(FUNCTION_NAME, "Array is not reduced");

  /* Get value */
  
  val = RES2[ptr];
  
  /* Return value */
  
  return val;
}

#endif

/*****************************************************************************/
