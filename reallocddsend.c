/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reallocddsend.c                                */
/*                                                                           */
/* Created:       2018/12/19 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Allocates memory for a new MPI_Isend() in a DDSend struct    */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReallocDDSend:"

/*****************************************************************************/

void ReallocDDSend(DDSend *send)
{

#ifdef MPI
  
  /* Nothing needs to be done if the MPI task is this task or -1 */
  
  if (send->mpiidto == mpiid || send->mpiidto == -1)
    return;
  
  /* Increase the total number of allocated buffers */
  
  (send->n)++;
  
  /* Allocate a new buffer */
  
  send->buffs = (double **)Mem(MEM_REALLOC, send->buffs, 
                               (send->n) * sizeof(double *));

  if (send->buffsize != 0)
    send->buffs[send->i] = (double *)Mem(MEM_ALLOC, send->buffsize, 
                                         sizeof(double));
  else
    send->buffs[send->i] = NULL;
  
  /* Allocate a new MPI_Request */
  
  send->reqs = (MPI_Request **)Mem(MEM_REALLOC, send->reqs, 
                                   (send->n) * sizeof(MPI_Request *));
  send->reqs[send->i] = (MPI_Request *)Mem(MEM_ALLOC, 1, sizeof(MPI_Request));
  
  /* Allocate a new MPI_Status */
  
  send->stats = (MPI_Status **)Mem(MEM_REALLOC, send->stats, 
                                   (send->n) * sizeof(MPI_Status *));
  send->stats[send->i] = (MPI_Status *)Mem(MEM_ALLOC, 1, sizeof(MPI_Status));
  
  /* Count this allocation */
  
  (send->size)++;

#endif
}

/*****************************************************************************/
