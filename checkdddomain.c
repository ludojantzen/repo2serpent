/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkdddomain.c                                */
/*                                                                           */
/* Created:       2019/03/29 (JLe)                                           */
/* Last modified: 2019/03/29 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Check that material belongs to current domain.               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckDDDomain:"

/*****************************************************************************/

long CheckDDDomain(long mat)
{
  long id;
  
  /* Check domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    return YES;
  
  /* Check material pointer */

  if (mat < VALID_PTR)
    return YES;

  /* Check if material is assigned with a domain ID */

  if ((id = (long)RDB[mat + MATERIAL_MPI_ID]) < 0)
    return YES;
  else if (id == mpiid)
    return YES;
  else
    return NO;
}

/*****************************************************************************/
