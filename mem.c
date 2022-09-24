/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : mem.c                                          */
/*                                                                           */
/* Created:       2011/12/15 (JLe)                                           */
/* Last modified: 2017/03/20 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Calls calloc(), realloc() or free() inside OpenMP critical   */
/*              clause.                                                      */
/*                                                                           */
/* Comments: - This function was written because some MPI implementations    */
/*             (OpenMPI) are not thread safe, and memory management routines */
/*             need to be protected manually. It was later discovered,       */
/*             however, that the use of nested omp critical pragmas causes   */
/*             deadlocks. The omp critical barrier is not used for now and   */
/*             problem must be fixed by using a thread-safe MPI              */
/*             implementation.                                               */
/*                                                                           */
/*           - None of the above is probably not true...                     */
/*                                                                           */
/*           - Pointers are not checked to enable debugging with valgrind.   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Mem:"

/*****************************************************************************/

void *Mem(long mode, ...)
{
  void *ptr0, *ptr1;
  size_t nelem, size;
  va_list argp;
  va_start (argp, mode);
  
  /* Set null return pointer */
  
  ptr1 = NULL;
  
  /* Check mode */
  
  if (mode == MEM_ALLOW)
    {
      /* Check pointer */
  
      if (RDB == NULL)
        Die(FUNCTION_NAME, "Pointer error");

      /* Check flag and set */

      if ((long)RDB[DATA_ALLOW_MEM_OP] == YES)
        Die(FUNCTION_NAME, "Memory allocation already allowed");
      else      
        WDB[DATA_ALLOW_MEM_OP] = (double)YES;
    }
  else if (mode == MEM_DENY)
    {
      /* Check pointer */
  
      if (RDB == NULL)
        Die(FUNCTION_NAME, "Pointer error");

      /* Set flag */

      if ((long)RDB[DATA_ALLOW_MEM_OP] == NO)
        Die(FUNCTION_NAME, "Memory allocation already denied");
      else      
        WDB[DATA_ALLOW_MEM_OP] = (double)NO;
    }
  else if (mode == MEM_ALLOC)
    {
      /* Get arguments */
      
      nelem = va_arg(argp, size_t);
      size = va_arg(argp, size_t);
      
      /* Check sizes */
      
      if (nelem < 1)
        Warn(FUNCTION_NAME, "Zero or negative block size (calloc)");
      if (size < 1)
        Die(FUNCTION_NAME, "Zero or negative element size (calloc)");        
      
      /* Allocate new block */
      
      if ((ptr1 = calloc(nelem, size)) == NULL)
        {
          if (RDB == NULL)
            Die(FUNCTION_NAME, "Memory allocation failed (calloc, %ld, %ld)",
                nelem, size);
          else
            Die(FUNCTION_NAME, 
                "Memory allocation failed (calloc, %ld, %ld, %1.2f)",
                nelem, size, RDB[DATA_REAL_BYTES]/MEGA);
        }
    }
  else if (mode == MEM_REALLOC)
    {
      /* Get arguments */
      
      ptr0 = va_arg(argp, void *);
      size = va_arg(argp, size_t);
      
      /* Check size */

      if (size < 1)
        Warn(FUNCTION_NAME, "Zero or negative block size (realloc)");
      
      /* Adjust size of existing block */
      
      if ((ptr1 = realloc(ptr0, size)) == NULL)
        Die(FUNCTION_NAME, "Memory allocation failed (realloc, size = %ld)",
            size);
      }
  else if (mode == MEM_FREE)
    {
      /* Get arguments */
      
      ptr0 = va_arg(argp, void *);
      
      /* Check pointer */
      
      if (ptr0 == NULL)
        Warn(FUNCTION_NAME, "Null pointer (free)");
      
      /* Free memory */
      
      free(ptr0);
    }
  else
    Die(FUNCTION_NAME, "Invalid mode");

  
  /* Return pointer */

  return ptr1;
}

/*****************************************************************************/
