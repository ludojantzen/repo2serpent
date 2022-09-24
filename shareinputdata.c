/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : shareinputdata.c                               */
/*                                                                           */
/* Created:       2010/11/23 (JLe)                                           */
/* Last modified: 2019/03/26 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Broadcasts input data to parallel MPI tasks                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ShareInputData:"

/*****************************************************************************/

void ShareInputData()
{

#ifdef MPI

  long sz, ompthreads;

  /* Check number of tasks */

  if (mpitasks == 1)
    return;

  fprintf(outp, "Sharing input data to MPI tasks...\n");

  /* Start timers */

  StartTimer(TIMER_MPI_OVERHEAD);
  StartTimer(TIMER_MPI_OVERHEAD_TOTAL);

  /***************************************************************************/

  /***** Main Data array *****************************************************/

  /* Get size of data block */

  if (mpiid == 0)
    sz = (long)RDB[DATA_REAL_MAIN_SIZE];

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Broadcast data size */

  if (MPI_Bcast(&sz, 1, MPI_LONG, 0, my_comm) != MPI_SUCCESS)
    Die(FUNCTION_NAME, "MPI Error");

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Allocate memory for data block in other tasks */

  if (mpiid > 0)
    {
      if ((WDB = Mem(MEM_ALLOC, sz, sizeof(double))) == NULL)
        Die(FUNCTION_NAME, "Cannot initialize main data array for task %ld",
            mpiid);
      else
        {
          /* Put read-only pointer */

          RDB = (const double *)WDB;
        }
    }

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Broadcast batch size */

  if (MPI_Bcast(&WDB[DATA_OPTI_MPI_BATCH_SIZE], 1, MPI_DOUBLE, 0,
                my_comm) != MPI_SUCCESS)
    Die(FUNCTION_NAME, "MPI Error");

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Transfer data */

  MPITransfer(WDB, NULL, sz, 0, MPI_METH_BC);

  /* Synchronise */

  MPI_Barrier(my_comm);

  /***************************************************************************/

  /***** ACE Data array ******************************************************/

  /* Get size of data block */

  sz = (long)RDB[DATA_REAL_ACE_SIZE];

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Allocate memory for data block in other tasks */

  if (mpiid > 0)
    if ((ACE = Mem(MEM_ALLOC, sz, sizeof(double))) == NULL)
      Die(FUNCTION_NAME, "Cannot initialize ACE data array for task %ld",
          mpiid);

  /* Transfer data */

  MPITransfer(ACE, NULL, sz, 0, MPI_METH_BC);

  /* Synchronise */

  MPI_Barrier(my_comm);

  /***************************************************************************/

  /***** ASCII Data array ****************************************************/

  /* Get size of data block */

  sz = (long)RDB[DATA_ASCII_DATA_SIZE];

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Allocate memory for data block in other tasks */

  if (mpiid > 0)
    if ((ASCII = Mem(MEM_ALLOC, sz, sizeof(char))) == NULL)
      Die(FUNCTION_NAME, "Cannot initialize ASCII array for task %ld",
          mpiid);

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Broadcast data block */

  if (MPI_Bcast(ASCII, sz, MPI_CHAR, 0, my_comm) != MPI_SUCCESS)
    Die(FUNCTION_NAME, "MPI Error");

  /* Synchronise */

  MPI_Barrier(my_comm);

  /***************************************************************************/

  /***** Private Data array **************************************************/

  /* Get number of OpenMP threads */

  ompthreads = (long)RDB[DATA_OMP_MAX_THREADS];

  /* Get size of data block */

  sz = (long)RDB[DATA_REAL_PRIVA_SIZE]*ompthreads;

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Allocate memory for data block in other tasks */

  if (mpiid > 0)
    if ((PRIVA = Mem(MEM_ALLOC, sz, sizeof(double))) == NULL)
      Die(FUNCTION_NAME, "Cannot initialize private data array for task %ld",
          mpiid);

  /* Transfer data */

  MPITransfer(PRIVA, NULL, sz, 0, MPI_METH_BC);

  /* Synchronise */

  MPI_Barrier(my_comm);

  /***************************************************************************/

  /***** Allocate memory for remaining ***************************************/

  /* First results array */

  sz = (long)RDB[DATA_REAL_RES1_SIZE];

  if (mpiid > 0)
    if ((RES1 = Mem(MEM_ALLOC, sz, sizeof(double))) == NULL)
      Die(FUNCTION_NAME, "Cannot initialize RES1 data array for task %ld",
          mpiid);

  memset(RES1, 0.0, sz*sizeof(double));

  /* Second results array */

  if ((long)RDB[DATA_OPTI_SHARED_RES2] == YES)
    sz = (long)RDB[DATA_REAL_RES2_SIZE];
  else
    sz = (long)RDB[DATA_REAL_RES2_SIZE]*ompthreads;

  if (mpiid > 0)
    if ((RES2 = Mem(MEM_ALLOC, sz, sizeof(double))) == NULL)
      Die(FUNCTION_NAME, "Cannot initialize RES2 data array for task %ld",
          mpiid);

  memset(RES2, 0.0, sz*sizeof(double));

  /* Scoring buffer */

  if ((long)RDB[DATA_OPTI_SHARED_BUF] == YES)
    sz = (long)RDB[DATA_REAL_BUF_SIZE];
  else
    sz = (long)RDB[DATA_REAL_BUF_SIZE]*ompthreads;

  if (mpiid > 0)
    if ((BUF = Mem(MEM_ALLOC, sz, sizeof(double))) == NULL)
      Die(FUNCTION_NAME, "Cannot initialize scoring buffer for task %ld",
          mpiid);

  memset(BUF, 0.0, sz*sizeof(double));

  /* Broadcast seed */

  if (MPI_Bcast(&parent_seed, 1, MPI_LONG, 0, my_comm) != MPI_SUCCESS)
    Die(FUNCTION_NAME, "MPI Error");

  /* Set number of OpenMP threads */

#ifdef OPEN_MP

  if (mpiid > 0)
    omp_set_num_threads((long)RDB[DATA_OMP_MAX_THREADS]);

#endif

  /* Synchronise */

  MPI_Barrier(my_comm);

  /***************************************************************************/

  fprintf(outp, "OK.\n\n");

  /* Stop timers */

  StopTimer(TIMER_MPI_OVERHEAD);
  StopTimer(TIMER_MPI_OVERHEAD_TOTAL);

#endif
}

/*****************************************************************************/
