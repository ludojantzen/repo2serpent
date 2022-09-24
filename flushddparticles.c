/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : flushddparticles.c                             */
/*                                                                           */
/* Created:       2018/12/19 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sends all half-full particle buffers.                        */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FlushDDParticles:"

/*****************************************************************************/

void FlushDDParticles()
{
  
#ifdef MPI

  long i;
  DDSend *send;
  
  /* Send the half-full buffers for all domains */
  
  for (i = 0; i < mpitasks; i++)
    {
      /* Get the DDSend struct */

      send = &dd_part_sends[i];
      
      /* Send the message if there is at least one particle in the buffer */
      
      if (send->j > 0)
        {
          /* Send the message */

          PostDDSend(send);
          
          /* Clean up the DDSend struct */
          
          CleanUpDDSend(send, 0, 0);
        }
    }

#endif
}

/*****************************************************************************/
