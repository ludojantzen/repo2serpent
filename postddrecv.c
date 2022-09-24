/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : postddrecv.c                                   */
/*                                                                           */
/* Created:       2018/12/19 (MGa)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Begins a nonblocking receive (MPI_Irecv())                   */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PostDDRecv:"

/*****************************************************************************/

void PostDDRecv(DDRecv *recv)
{

#ifdef MPI

  /* Nothing needs to be done if the MPI task is this task or -1 */

  if (recv->mpiidfrom == mpiid || recv->mpiidfrom == -1)
    return;

  /* Post the receive */

  MPI_Irecv(recv->buff, recv->buffsize, MPI_DOUBLE, recv->mpiidfrom,
            recv->tag, my_comm, &recv->req);

#endif

}

/*****************************************************************************/
