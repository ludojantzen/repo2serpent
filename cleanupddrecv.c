/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : cleanupddrecv.c                                */
/*                                                                           */
/* Created:       2018/12/19 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Cleans up a DDRecv struct                                    */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*           - If needed, receives incoming messages until an empty message  */
/*             is received (this is used to make sure that all messages sent */
/*             with MPI_Isend() are received and no resources are leaked).   */
/*           - If no extra messages need to be received, the MPI_Irecv() is  */
/*             just canceled.                                                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CleanUpDDRecv:"

/*****************************************************************************/

void CleanUpDDRecv(DDRecv *recv, long cancel)
{
  
#ifdef MPI
  
  int flag, n;

  /* Nothing needs to be done if the MPI task is this task or -1 */
  
  if (recv->mpiidfrom == mpiid || recv->mpiidfrom == -1)
    return;
  
  /* Check if the receive just needs to be canceled */
  
  if (cancel)
    {
      /* Cancel the receive */

      MPI_Cancel(&recv->req);
      
      /* Nothing else to do */
      
      return;
    }

  /* Loop until the empty message is received */

  while (1)
    {
      /* Wait for the message and get its size */
      
      CheckDDRecv(recv, 1, &flag, &n);
      
      /* Break if the message is empty */
      
      if (n == 0)
        break;
      
      /* Post a new asynchronous receive for this message counter message */
      
      PostDDRecv(recv);
    }

  #endif
}

/*****************************************************************************/
