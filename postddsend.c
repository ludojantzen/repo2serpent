/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : postddsend.c                                   */
/*                                                                           */
/* Created:       2018/12/19 (MGa)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Begins a nonblocking send (MPI_Isend())                      */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PostDDSend:"

/*****************************************************************************/

void PostDDSend(DDSend *send)
{

#ifdef MPI

  /* Nothing needs to be done if the MPI task is this task or -1 */

  if (send->mpiidto == mpiid || send->mpiidto == -1)
    return;

  /* Send the message */

  MPI_Isend(send->buffs[send->i], send->j, MPI_DOUBLE, send->mpiidto,
            send->tag, my_comm, send->reqs[send->i]);

  /* Increase the counter to the next buffer */

  (send->i)++;

  /* Reset the counter to the first position in the new buffer */

  send->j = 0;

#endif

}

/*****************************************************************************/
