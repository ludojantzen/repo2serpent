/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : initddrecv.c                                   */
/*                                                                           */
/* Created:       2018/12/20 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Initializes the data in a DDRecv struct                      */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitDDRecv:"

/*****************************************************************************/

void InitDDRecv(DDRecv *recv, long mpiidfrom, long tag, long buffsize)
{
  
#ifdef MPI
  
  /* Initialize the data */

  recv->mpiidfrom = mpiidfrom;
  if (mpiidfrom == mpiid || mpiidfrom == -1)
    {
      recv->tag = -1;
      recv->buffsize = -1;
    }
  else
    {
      recv->tag = tag;
      recv->buffsize = buffsize;
    }
  
  /* Allocate the buffer if needed */
  
  if (mpiidfrom == mpiid || mpiidfrom == -1 || buffsize == 0)
    recv->buff = NULL;
  else
    recv->buff = (double *)Mem(MEM_ALLOC, buffsize, sizeof(double));

#endif

}

/*****************************************************************************/
