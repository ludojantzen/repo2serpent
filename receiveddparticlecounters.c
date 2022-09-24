/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : receiveddparticlecounters.c                    */
/*                                                                           */
/* Created:       2018/06/13 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Receives the particle counters from the children MPI tasks   */
/*              during the asynchronous estimation of the tracking           */
/*              termination.                                                 */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReceiveDDParticleCounters:"

/*****************************************************************************/

void ReceiveDDParticleCounters()
{

#ifdef MPI
  
  long i;
  int flag, n;
  DDRecv *recv;
  
  /* Receive messages coming from each child task */

  for (i = 0; i < 2; i++)
    {
      /* Get the DDRecv struct */
      
      recv = &dd_part_count_recvs[i];
      
      /* Loop until there's no message coming */

      while (1)
        {
          /* Check if a message has arrived from this task */

          CheckDDRecv(recv, 0, &flag, &n);
          
          /* Break the loop if there's no message */
          
          if (!flag)
            break;
          
          /* Get the particle counter */

          dd_part_count_children[i] = (long)recv->buff[0];
          
          /* Post a new receive for this task */
          
          PostDDRecv(recv);
        }
    }

#endif
  
}

/*****************************************************************************/
