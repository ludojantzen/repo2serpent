/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkfinishddasynch.c                          */
/*                                                                           */
/* Created:       2018/06/13 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Checks for termination of the tracking loop in DD mode using */
/*              asynchronous MPI communications.                             */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*           - Estimates if maybe all particles sent have been received.     */
/*           - Depending on this estimation the synchronous check is done.   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckFinishDDAsynch:"

/*****************************************************************************/

void CheckFinishDDAsynch(long *maybe_done)
{
  long i, dd_part_count_global;
  int flag, n;
  DDRecv *recv;
  DDSend *send;

  /* Check if domain decomposition is on */

  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    {
      *maybe_done = 1;
      return;
    }

  /* Send the particle message counter to the parent task */

  SendDDParticleCounter();

  /* Receive the particle message counters from the children tasks */

  ReceiveDDParticleCounters();

  /* Check the particle balance in task 0 */

  if (mpiid == 0)
    {

      /* Get the global particle balance */

      dd_part_count_global = dd_part_count + dd_part_count_children[0] + 
        dd_part_count_children[1];

      /* Check if all particles sent might have been received */

      *maybe_done = dd_part_count_global == 0;

    }

  /* Check for synchronization requests in the other tasks */

  else
    {

      /* Get the DDRecv struct */

      recv = &dd_synch_req_recv;

      /* Check if a synchronization request message has been received */

      CheckDDRecv(recv, 0, &flag, &n);

      /* The tracking might be over if the parent task says so */

      *maybe_done = flag;

      /* Post a new asynchronous receive if a message was received */

      if (flag)
        PostDDRecv(recv);

    }

  /* Send the synchronization request to the children tasks */

  if (*maybe_done)
    for (i = 0; i < 2; i++)
      {

        /* Get the DDSend struct */

        send = &dd_synch_req_sends[i];

        /* Send the synchronization request */

        ReallocDDSend(send);
        PostDDSend(send);

      }

}

/*****************************************************************************/
