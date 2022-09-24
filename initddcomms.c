/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : initddcomms.c                                  */
/*                                                                           */
/* Created:       2018/06/13 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Initializes the data used for the communications in DD mode. */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*           - The memory allocated here is freed in FreeDDComms().          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitDDComms:"

/*****************************************************************************/

void InitDDComms()
{
  long i;

  /* Check if domain decomposition is on */

  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    return;

  /* Get the parent in the binary tree for asynchronous termination check */
  /* Note: -1 means that the parent does not exist. */

  if (mpiid > 0)
    dd_mpiid_parent = (mpiid-1) / 2;
  else
    dd_mpiid_parent = -1;

  /* Get the children in the binary tree for asynchronous termination check */
  /* Note: -1 means that the child does not exist. */

  for (i = 0; i < 2; i++)
    {
      dd_mpiid_children[i] = 2*mpiid + i + 1;
      if (dd_mpiid_children[i] > mpitasks-1) dd_mpiid_children[i] = -1;
    }

  /* Set the size of the data block for each particle */
  /* Note: history data needed for ifp is excluded for now */
  /* (DATA_IFP_CHAIN_LENGTH set to -1). */

  dd_part_size = PARTICLE_BLOCK_SIZE - LIST_DATA_SIZE;

  /* Set the size of the particle buffers */

  dd_buff_size = (long)RDB[DATA_DD_BUFF_NUM_PARTS] * dd_part_size;

  /* Allocate the DDSend structs to manage particle messages */

  dd_part_sends = (DDSend *)Mem(MEM_ALLOC, mpitasks, sizeof(DDSend));

  /* Allocate and initialize the DDRecv structs to manage particle messages */
  /* Note: the i-th DDRecv receives particles from domain i, with tag i. */

  dd_part_recvs = (DDRecv *)Mem(MEM_ALLOC, mpitasks, sizeof(DDRecv));
  for (i = 0; i < mpitasks; i++)
    InitDDRecv(&dd_part_recvs[i], i, i, dd_buff_size);

  /* Initialize the DDRecv structs to manage particle counter messages */
  /* Note: the i-th DDRecv receives particle counter messages from children */
  /* i in the binary tree, with tag mpitasks + dd_mpiid_children[i]. */

  for (i = 0; i < 2; i++)
    InitDDRecv(&dd_part_count_recvs[i], dd_mpiid_children[i], 
      mpitasks + dd_mpiid_children[i], 1);

  /* Initialize the DDRecv struct to manage synchronization request messages */
  /* Note: this DDRecv receives synchronization request (empty) messages */
  /* from the parent in the binary tree, with tag 2*mpitasks + */
  /* dd_mpiid_parent. */

  InitDDRecv(&dd_synch_req_recv, dd_mpiid_parent, 
    2*mpitasks + dd_mpiid_parent, 0);

}

/*****************************************************************************/
