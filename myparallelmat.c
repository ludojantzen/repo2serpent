/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : myparallelmat                                  */
/*                                                                           */
/* Created:       2012/06/04 (JLe)                                           */
/* Last modified: 2012/06/04 (JLe)                                           */
/* Version:       2.1.6                                                      */
/*                                                                           */
/* Description: Checks OpenMP and MPI id's and checks whether this material  */
/*              should be processed by this task and thread                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MyParalleMat:"

/*****************************************************************************/

long MyParallelMat(long mat, long mpi)
{
  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Check MPI task number */

  if (mpi == YES)
    if (((long)RDB[mat + MATERIAL_MPI_ID] != -1) &&
        ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid))
      return NO;

  /* Grad material to self if thread number is not set */

#ifdef OPEN_MP
#pragma omp critical
#endif
  {
    /* Grab material */
              
    if ((long)RDB[mat + MATERIAL_OMP_ID] < 0)
      WDB[mat + MATERIAL_OMP_ID] = OMP_THREAD_NUM;
  }

  /* Check thread number */
  
  if ((long)RDB[mat + MATERIAL_OMP_ID] != OMP_THREAD_NUM)
    return NO;

  /* OK to handle this material */

  return YES;
}

/*****************************************************************************/
