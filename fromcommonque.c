/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fromcommonque.c                                */
/*                                                                           */
/* Created:       2019/02/21 (JLe)                                           */
/* Last modified: 2019/02/21 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Retrieves neutron / photon from common que                   */
/*                                                                           */
/* Comments: Used for load balancing, e.g. production of secondary photons.  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FromCommonQue:"

/*****************************************************************************/

long FromCommonQue(long id)
{
  long ptr, new, nt, n;

  /* Get pointer */

  ptr = (long)RDB[DATA_PART_PTR_COMMON_QUE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get number of OpenMP threads */

  nt = (long)RDB[DATA_OMP_MAX_THREADS];
  CheckValue(FUNCTION_NAME, "nt", "", nt, 1, 100000000);

#ifdef OPEN_MP
#pragma omp critical (que)
#endif
  {
    /* Get last particle */
    
    ptr = LastItem(ptr);
    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
    
    /* Check type */
    
    if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_DUMMY)
      {
        /* Que is empty, return null */
        
        ptr = -1;
      }
    else
      {
        /* Get multiplicity */

        n = (long)RDB[ptr + PARTICLE_MULTIPLICITY];

        /* Check multiplicity */
    
        if ((nt > 1) && (n > nt))
          {
            /* Subtract value */
            
            WDB[ptr + PARTICLE_MULTIPLICITY] = 
              RDB[ptr + PARTICLE_MULTIPLICITY] - nt;
            
            /* Duplicate */
            
            new = DuplicateParticle(ptr, id);
            WDB[new + PARTICLE_MULTIPLICITY] = nt - 2;

            /* Add item in list */

            AddItem(OMPPtr(DATA_PART_PTR_QUE, id), new);

            /* Duplicate */

            ptr = DuplicateParticle(new, id);
                  WDB[ptr + PARTICLE_MULTIPLICITY] = 0.0;
          }
        else if (n > 0)
          {
            /* Subtract value */
            
            WDB[ptr + PARTICLE_MULTIPLICITY] = 
              RDB[ptr + PARTICLE_MULTIPLICITY] - 1.0;
            
            /* Duplicate */
            
            ptr = DuplicateParticle(ptr, id);
            WDB[ptr + PARTICLE_MULTIPLICITY] = 0.0;      
          }
        else
          {
            /* Remove particle from que */
            
            RemoveItem(ptr);
          }
      }
  }

  /* Return pointer */
  
  return ptr;
}

/*****************************************************************************/
