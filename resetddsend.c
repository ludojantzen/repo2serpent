/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : resetddsend.c                                  */
/*                                                                           */
/* Created:       2018/12/20 (MGa)                                           */
/* Last modified: 2019/01/03 (MGa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Resets the data in a DDSend struct                           */
/*                                                                           */
/* Comments: - Used with domain decomposition.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ResetDDSend:"

/*****************************************************************************/

void ResetDDSend(DDSend *send, long mpiidto, long tag, long buffsize)
{
  
#ifdef MPI
  
  /* Reset the data */
  
  send->i = 0;
  send->j = 0;
  send->n = 0;
  send->size = 0;
  send->mpiidto = mpiidto;

  if (mpiidto == mpiid || mpiidto == -1)
    {
      send->tag = -1;
      send->buffsize = -1;
    }
  else
    {
      send->tag = tag;
      send->buffsize = buffsize;
    }
  
  send->buffs = NULL;
  send->reqs = NULL;
  send->stats = NULL;

#endif
  
}

/*****************************************************************************/
