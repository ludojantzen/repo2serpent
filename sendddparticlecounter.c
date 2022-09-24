/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sendddparticlecounter.c                        */
/*                                                                           */
/* Created:       2018/06/13 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sends the particle counter to the parent MPI task during the */
/*              asynchronous estimation of the tracking termination.         */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SendDDParticleCounter:"

/*****************************************************************************/

void SendDDParticleCounter()
{
  
#ifdef MPI
  
  long dd_part_count_total;
  DDSend *send;

  /* Return if there's no parent task */
  
  if (dd_mpiid_parent == -1)
    return;
  
  /* Get the DDSend struct */
  
  send = &dd_part_count_send;
  
  /* Get the total particle counter (local + children) */
  
  dd_part_count_total = dd_part_count + dd_part_count_children[0] + 
    dd_part_count_children[1];
  
  /* Return if the counter hasn't changed since the last send */

  if (send->i > 0)
    if (dd_part_count_total == (long)(send->buffs[send->i-1][0]))
      return;
  
  /* Allocate new data structures */
  
  ReallocDDSend(send);
  
  /* Get the total particle counter (local + children) so that the global */
  /* sum moves up the binary tree */
  
  send->buffs[send->i][0] = (double)dd_part_count_total;
  
  /* Increase the counter to the next position in the buffer */

  send->j = 1;
  
  /* Send the total particle counter */

  PostDDSend(send);

#endif
  
}

/*****************************************************************************/
