/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkddrecv.c                                  */
/*                                                                           */
/* Created:       2018/12/19 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Checks a nonblocking receive (MPI_Irecv())                   */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckDDRecv:"

/*****************************************************************************/

void CheckDDRecv(DDRecv *recv, long wait, int *flag, int *n)
{

#ifdef MPI
  
  /* Nothing needs to be done if the MPI task is this task or -1 */
  
  if (recv->mpiidfrom == mpiid || recv->mpiidfrom == -1)
    {
      *flag = 0;
      return;
    }

  /* Check if the recieve has been completed or wait for it */
  
  if (!wait)
    MPI_Test(&recv->req, flag, &recv->stat);
  else
    {
      MPI_Wait(&recv->req, &recv->stat);
      *flag = 1;
    }
  
  /* Get the size of the message, if completed */
  
  if (*flag)
    MPI_Get_count(&recv->stat, MPI_DOUBLE, n);
  
#endif
}

/*****************************************************************************/
