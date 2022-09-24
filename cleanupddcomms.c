/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : cleanupddcomms.c                               */
/*                                                                           */
/* Created:       2018/06/13 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Cleans up the data used for the communications in DD mode    */
/*              allocated during each transport cycle, after the tracking is */
/*              finished.                                                    */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CleanUpDDComms:"

/*****************************************************************************/

void CleanUpDDComms()
{
  long i;

  /* Check if domain decomposition is on */

  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    return;

  /* Clean up the DDSend structs to manage particle messages */
  /* Note: all sends must be completed by now. */

  for (i = 0; i < mpitasks; i++)
    CleanUpDDSend(&dd_part_sends[i], 0, 1);

  /* Clean up the DDRecv structs to manage particle messages */
  /* Note: all receives must be completed by now. */

  for (i = 0; i < mpitasks; i++)
    CleanUpDDRecv(&dd_part_recvs[i], 1);

  /* Send an empty message counter to make sure all sends are completed */
  /* Note: this is needed to avoid resource leakage and messages being */
  /* mistakenly received in the next cycle. */

  ReallocDDSend(&dd_part_count_send);
  PostDDSend(&dd_part_count_send);

  /* Clean up the DDRecv structs to manage particle counter messages */
  /* Note: all pending receives are completed here. */

  for (i = 0; i < 2; i++)
    CleanUpDDRecv(&dd_part_count_recvs[i], 0);

  /* Clean up the DDSend struct to manage particle counter messages */
  /* Note: all pending sends are completed here. */

  CleanUpDDSend(&dd_part_count_send, 1, 1);

  /* Clean up the DDSend structs to manage synchronization request messages */
  /* Note: all sends must be completed by now. */

  for (i = 0; i < 2; i++)
    CleanUpDDSend(&dd_synch_req_sends[i], 0, 1);

  /* Clean up the DDRecv struct to manage synchronization request messages */
  /* Note: all receives must be completed by now. */

  CleanUpDDRecv(&dd_synch_req_recv, 1);

}

/*****************************************************************************/
