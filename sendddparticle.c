/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sendddparticle.c                               */
/*                                                                           */
/* Created:       2018/12/18 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sends a particle to a given domain.                          */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*           - Particles are buffered and sent when the buffers are full.    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SendDDParticle:"

/*****************************************************************************/

void SendDDParticle(DDSend *send, long part)
{

#ifdef MPI
  
  /* Allocate new data structures if this is the first particle of a buffer */
  
  if (send->j == 0)
    ReallocDDSend(send);
  
  /* Add the particle to the buffer */
  
  memcpy(&send->buffs[send->i][send->j], &RDB[part + LIST_DATA_SIZE], 
         dd_part_size * sizeof(double));
  
  /* Count this particle */
  
  dd_part_count++;
  
  /* Increase the counter to the next position in the buffer */
  
  send->j += dd_part_size;

  /* Send the message if the buffer is full */
  
  if (send->j == dd_buff_size)
    {
      /* Send the message */
      
      PostDDSend(send);
      
      /* Clean up the DDSend struct */
      
      CleanUpDDSend(send, 0, 0);
    }

#endif
  
}

/*****************************************************************************/
