/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getprivatedata.c                               */
/*                                                                           */
/* Created:       2011/11/11 (JLe)                                           */
/* Last modified: 2014/04/04 (JLe)                                           */
/* Version:       2.1.20                                                     */
/*                                                                           */
/* Description: Gets value from PRIVA data block                             */
/*                                                                           */
/* Comments: - Replaced by macro when not compiled in debug mode             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetPrivateData:"

/*****************************************************************************/

#ifdef DEBUG 

double GetPrivateData(long ptr, long id)
{
  long sz;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "PRIVA array not ready for access");

  /* Get size of data block */
  
  sz = (long)RDB[DATA_REAL_PRIVA_SIZE];
  
  /* Check pointer */

  if ((ptr < 0) || (ptr > sz - 1))
    Die(FUNCTION_NAME, "Pointer error");

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");
  
  /* Get value */
  
  return PRIVA[ptr + id*sz];
}

#endif

/*****************************************************************************/
