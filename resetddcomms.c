/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : resetddcomms.c                                 */
/*                                                                           */
/* Created:       2018/06/13 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Resets the data used for the communications in DD mode.      */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*           - This is done before every transport cycle.                    */
/*           - Asynchronous receives (MPI_Irecv()) are posted for all        */
/*             posible incoming messages                                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ResetDDComms:"

/*****************************************************************************/

void ResetDDComms()
{
  long i;

  /* Check if domain decomposition is on */

  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    return;

  /* Reset the particle counters */

  dd_part_count = 0;
  for (i = 0; i < 2; i++)
    dd_part_count_children[i] = 0;

  /* Reset the DDSend structs to manage particle messages */
  /* Note: the i-th DDSend sends particles to domain i, with tag mpiid. */

  for (i = 0; i < mpitasks; i++)
    ResetDDSend(&dd_part_sends[i], i, mpiid, dd_buff_size);

  /* Post asynchronous receives for particle messages */

  for (i = 0; i < mpitasks; i++)
    PostDDRecv(&dd_part_recvs[i]);

  /* Reset the DDSend struct to manage particle counter messages */
  /* Note: this DDSend sends particle counter messages to the parent in the */
  /* binary tree, with tag mpitasks + mpiid. */

  ResetDDSend(&dd_part_count_send, dd_mpiid_parent, mpitasks + mpiid, 1);

  /* Post asynchronous receives for particle counter messages */

  for (i = 0; i < 2; i++)
    PostDDRecv(&dd_part_count_recvs[i]);

  /* Reset the DDSend structs to manage synchronization request messages */
  /* Note: the i-th DDSend sends synchronization request messages to */
  /* children i in the binary tree, with tag 2*mpitasks + mpiid. */

  for (i = 0; i < 2; i++)
    ResetDDSend(&dd_synch_req_sends[i], dd_mpiid_children[i], 
      2*mpitasks + mpiid, 0);

  /* Post an asynchronous receive for synchronization request messages */

  PostDDRecv(&dd_synch_req_recv);

}

/*****************************************************************************/
