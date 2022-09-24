/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : distributeddlimbo.c                            */
/*                                                                           */
/* Created:       2018/06/13 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sends the particles going to other domains in DD mode.       */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DistributeDDLimbo:"

/*****************************************************************************/

void DistributeDDLimbo()
{
  long id, part, mat, mpiidto, ptr;
  DDSend *send;

  /* Check if domain decomposition is on */

  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    return;

  /* Start overhead timer */

  StartTimer(TIMER_DD_OVERHEAD);
  
  /* Loop over OpenMP threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get pointer to limbo */

      part = (long)RDB[OMPPtr(DATA_PART_PTR_LIMBO, id)];
      CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

      /* Get pointer to first item */

      part = FirstItem(part);
      CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

      /* Check particle type */

      if ((long)RDB[part + PARTICLE_TYPE] != PARTICLE_TYPE_DUMMY)
        Die(FUNCTION_NAME, "Error in list");

      /* Get pointer to first particle */

      part = NextItem(part);

      /* Loop over particles in this thread's limbo */

      while (part > VALID_PTR)
        {
          /* Get material pointer */

          mat = (long)RDB[part + PARTICLE_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
          
          /* Get MPI id */
          
          mpiidto = (long)RDB[mat + MATERIAL_MPI_ID];
          CheckValue(FUNCTION_NAME, "mpiidto", "", mpiidto, 0, mpitasks - 1);
          
          /* Check MPI id */
          
          if (mpiidto == mpiid)
            Die(FUNCTION_NAME, "Error in index");
          
          /* Get the DDSend struct */
          
          send = &dd_part_sends[mpiidto];
          
          /* Send particle */

          SendDDParticle(send, part);
          
          /* Copy pointer */
          
          ptr = part;
          
          /* Pointer to next */
          
          part = NextItem(part);
          
          /* Remove particle from limbo */
          
          RemoveItem(ptr);
          
          /* Put particle into stack */
          
          ToStack(ptr, id);                    
        }
    }
  
  /* Check that limbos are empty */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get pointer to limbo */

      part = (long)RDB[OMPPtr(DATA_PART_PTR_LIMBO, id)];
      CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
      
      /* Check size */
      
      if (ListSize(part) != 1)
        Die(FUNCTION_NAME, "Limbo is not empty after sending particles");
      
    }
  
  /* Send all half-full particle buffers */

  FlushDDParticles();
  
  /* Stop overhead timer */

  StopTimer(TIMER_DD_OVERHEAD);
}

/*****************************************************************************/
