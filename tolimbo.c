/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : tolimbo.c                                      */
/*                                                                           */
/* Created:       2018/03/23 (JLe)                                           */
/* Last modified: 2018/04/01 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Puts particle in limbo                                       */
/*                                                                           */
/* Comments: - Used only with domain decomposition                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ToLimbo:"

/*****************************************************************************/

void ToLimbo(long ptr, long id)
{
  long loc0;
  
  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Check if domain decomposition is on */

  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    Die(FUNCTION_NAME, "Domain decomposition is not on");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check MPI id */

  if ((long)RDB[ptr + PARTICLE_MPI_ID] != mpiid)
    Die(FUNCTION_NAME, "Error in particle MPI id");

  /* Add to balance */

  loc0 = (long)RDB[RES_DD_BALA_IN];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
  AddBuf(1.0, 1.0, loc0, id, -1, mpiid, (long)RDB[ptr + PARTICLE_TYPE] - 1);
  
  /* Add item in list */

  AddItem(OMPPtr(DATA_PART_PTR_LIMBO, id), ptr);
}

/*****************************************************************************/

