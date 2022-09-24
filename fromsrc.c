/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fromsrc.c                                      */
/*                                                                           */
/* Created:       2011/03/10 (JLe)                                           */
/* Last modified: 2015/04/08 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Retrieves neutron from source                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FromSrc:"

/*****************************************************************************/

long FromSrc(long id)
{
  long ptr, pts;
  unsigned long seed;

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Put barrier */
      
#ifdef OPEN_MP
#pragma omp critical
#endif
  
  {
    /* Get pointer to source distribution */
  
    ptr = (long)RDB[DATA_PART_PTR_SOURCE];
    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

    /* Get pointer to last item */
    
    ptr = LastItem(ptr);
    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

    /* Check type */
    
    if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_DUMMY)
      ptr = -1;
    else
      {
        /* Remove particle from source */

        RemoveItem(ptr);
      }
  }

  /* Check if source is empty */

  if (ptr < VALID_PTR)
    return -1;

  /* Check time interval */

  if (RDB[ptr + PARTICLE_T] < RDB[DATA_TIME_CUT_TMIN])
    Die(FUNCTION_NAME, "Error in time");
  else if (RDB[ptr + PARTICLE_T] >= RDB[DATA_TIME_CUT_TMAX])
    {
      /* Bank particle */

      ToBank(ptr, id);         
    }
  else
    {
      /* Init random number sequence */

      seed = ReInitRNG((long)RDB[ptr + PARTICLE_RNG_IDX]);
      SEED[id*RNG_SZ] = seed;
          
      /* Put neutron in que */
          
      ToQue(ptr, id);   
          
      /* Score initial source weight */
      
      if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT) &&
          ((long)RDB[ptr + PARTICLE_MPI_ID] == mpiid))
        {
          pts = (long)RDB[RES_INI_SRC_WGT];
          CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);
          AddBuf1D(1.0, RDB[ptr + PARTICLE_WGT], pts, id, 0);
        }  
    }

  /* Return pointer */
  
  return ptr;
}

/*****************************************************************************/
