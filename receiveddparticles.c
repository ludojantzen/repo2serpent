/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : receiveddparticles.c                           */
/*                                                                           */
/* Created:       2018/06/13 (MGa)                                           */
/* Last modified: 2019/03/28 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Receive the particles coming from other domains in DD mode.  */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReceiveDDParticles:"

/*****************************************************************************/

void ReceiveDDParticles(long *no_particles_left)
{
  
#ifdef MPI
  
  long part, ptr, id, mat, i, j, k, np, type;
  int flag, n;
  double wgt;
  DDRecv *recv;
  
  /* Check if domain decomposition is on */
  
  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    {
      *no_particles_left = 1;
      return;
    }
  
  /* Reset the flag to signal if there are particles left to track */
  
  *no_particles_left = 1;
  
  /* Loop over domains checking if particles are coming */
  
  for (i = 0; i < mpitasks; i++)
    {
      /* Get the DDRecv struct */
      
      recv = &dd_part_recvs[i];
      
      /* Check if new neutrons have arrived from this domain */
      
      CheckDDRecv(recv, 0, &flag, &n);
      
      /* Process the message, if detected */

      if (flag)
        {
          /* The new particles have to be tracked */
          
          *no_particles_left = 0;
          
          /* Get the number of particles */
          
          np = n / dd_part_size;
          
          /* Count this particles */
          
          dd_part_count -= np;
          
          /* Reset OpenMP thread id */
          
          id = 0;
          
          /* Loop over particles */
          
          for (j = 0; j < np; j++)
            {
              /* Get particle type */

              k = j*dd_part_size + PARTICLE_TYPE - LIST_DATA_SIZE;
              type = recv->buff[k];
              
              /* This should be neutron now (this check also serves to check */
              /* that the transfer and indexing was done correctly) */
              
              if (type != PARTICLE_TYPE_NEUTRON)
                Die(FUNCTION_NAME, "Error in type (remove this check later)");
              
              /* Get particle from stack */
              
              part = FromStack(type, id);
              
              /* Copy data to particle structure */
              
              memcpy(&WDB[part + LIST_DATA_SIZE], &recv->buff[j*dd_part_size], 
                     dd_part_size*sizeof(double));

              /* Get material pointer */
              
              mat = (long)RDB[part + PARTICLE_PTR_MAT];
              CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
              
              /* Check material MPI id */
              
              if ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid)
                Die(FUNCTION_NAME, "Error in MPI id");
              
              /* Check particle MPI id */
              
              if ((long)RDB[part + PARTICLE_MPI_ID] == mpiid)
                Die(FUNCTION_NAME, "Error in MPI id");

              /* Get particle weight */
              
              wgt = RDB[part + PARTICLE_WGT];
              
              /* Add to balance */

              ptr = (long)RDB[RES_DD_BALA_OUT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(1.0, 1.0, ptr, id, -1, mpiid,
                     (long)RDB[part + PARTICLE_TYPE] - 1);

              /* Put particle to que */
              
              ToQue(part, id);
              
              /* Update thread id */
              
              if (id++ == (long)RDB[DATA_OMP_MAX_THREADS]-1)
                id = 0;
            }

          /* Post a new receive for this domain */
          
          PostDDRecv(recv);
        }
    }

#endif
  
}

/*****************************************************************************/
