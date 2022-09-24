/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : cleanupddsend.c                                */
/*                                                                           */
/* Created:       2018/12/19 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Cleans up a DDSend struct                                    */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*           - Frees the data for completed sends (MPI_Isend()).             */
/*           - Checks that everything's fine when the tracking is finished.  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CleanUpDDSend:"

/*****************************************************************************/

void CleanUpDDSend(DDSend *send, long wait, long done)
{

#ifdef MPI
  
  long i;
  int flag;

  /* Nothing needs to be done if the MPI task is this task or -1 */

  if (send->mpiidto == mpiid || send->mpiidto == -1)
    return;
  
  /* Free the data structures for completed sends */
  
  for (i = 0; i < send->n; i++)
    {
      /* Check if this buffer is still allocated */

      if (send->stats[i] != NULL)
        {
          /* Check if the send has been completed or wait for it */

          if (!wait)
            MPI_Test(send->reqs[i], &flag, send->stats[i]);
          else
            {
              MPI_Wait(send->reqs[i], send->stats[i]);
              flag = 1;
            }
          
          /* If the send was completed deallocate the data */

          if (flag)
            {
              /* Free the data */

              if (send->buffsize != 0)
                Mem(MEM_FREE, send->buffs[i]);

              Mem(MEM_FREE, send->reqs[i]);
              Mem(MEM_FREE, send->stats[i]);
              send->buffs[i] = NULL;
              send->reqs[i] = NULL;
              send->stats[i] = NULL;

              /* Count this free */
              
              (send->size)--;
            }
          
          /* Terminate if the send needed to be completed */
          
          else if (done)
            Die(FUNCTION_NAME, "Message lost!");
        }
    }
  
  /* Nothing else needs to be done if the tracking is not finished */
  
  if (!done)
    return;
  
  /* Check if all the data has been freed */
  
  if (send->size != 0)
    Die(FUNCTION_NAME, "Memory management error in the DDSend structs!");
  
  /* Free the main data structures */
  
  if (send->stats != NULL)
    {
      /* Free the buffer, request and status arrays: */

      Mem(MEM_FREE, send->buffs);
      Mem(MEM_FREE, send->reqs);
      Mem(MEM_FREE, send->stats);
    }

#endif
}

/*****************************************************************************/
