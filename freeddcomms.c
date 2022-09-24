/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : freeddcomms.c                                  */
/*                                                                           */
/* Created:       2018/06/13 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Frees the data used for the communications in DD mode.       */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*           - Any allocated receive buffers are freed here.                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FreeDDComms:"

/*****************************************************************************/

void FreeDDComms()
{
#ifdef MPI

  long i;

  /* Check if domain decomposition is on */
  
  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    return;
  
  /* Free the DDSend structs to manage particle messages */
  
  Mem(MEM_FREE, dd_part_sends);
  
  /* Free the DDRecv structs to manage particle messages */
  
  for (i = 0; i < mpitasks; i++)
    if (dd_part_recvs[i].buff != NULL)
      Mem(MEM_FREE, dd_part_recvs[i].buff);

  Mem(MEM_FREE, dd_part_recvs);

  /* Free the DDRecv structs to manage particle counter messages */
  
  for (i = 0; i < 2; i++)
    if (dd_part_count_recvs[i].buff != NULL)
      Mem(MEM_FREE, dd_part_count_recvs[i].buff);
  
  /* Free the DDRecv struct to manage synchronization request messages */
  
  if (dd_synch_req_recv.buff != NULL)
    Mem(MEM_FREE, dd_synch_req_recv.buff);

#endif
  
}

/*****************************************************************************/
