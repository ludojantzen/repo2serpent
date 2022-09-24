/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : workarray.c                                    */
/*                                                                           */
/* Created:       2012/05/31 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Returns pointer to pre-allocated working array               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WorkArray:"

/*****************************************************************************/

double *WorkArray(long root, long type, long np, long id)
{
  long ptr, sz, loc0;
  double *dat;

  /* Check root pointer */

  if ((ptr = (long)RDB[root]) > VALID_PTR)
    {
      /* Existing data, check size */

      if (np > (long)RDB[ptr++])
        {
          /* Insufficient array size */

          Die(FUNCTION_NAME, "Array size exceeds pre-allocated data");
        }

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Avoid compiler warning */

      dat = NULL;

      /* Check type */

      if (type == DATA_ARRAY)
        {
          /* Check id */

          if (id > 0)
            Die(FUNCTION_NAME, "Error in thread id");

          /* Get pointer to data */

          dat = &WDB[ptr];
        }
      else if (type == PRIVA_ARRAY)
        {
          /* Check if access is allowed */

          if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
            Die(FUNCTION_NAME, "PRIVA array not ready for access");

          /* Get pointer to second pointer */

          ptr = (long)RDB[ptr];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get size of private data block */
  
          sz = (long)RDB[DATA_REAL_PRIVA_SIZE];
  
          /* Check pointer */

          if ((ptr < 0) || (ptr > sz - 1))
            Die(FUNCTION_NAME, "Pointer error");
          
          /* Check id */
          
          if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
            Die(FUNCTION_NAME, "Error in thread id");

          /* Get pointer to data */

          dat = &PRIVA[ptr + id*sz];
        }
      else
        Die(FUNCTION_NAME, "Invalid array type");

      /* Reset data */
      
      memset(dat, 0.0, np*sizeof(double));
    }
  else
    {
      /* No pre-existing data, check thread id */
      
      if (id > 0)
        Die(FUNCTION_NAME, "First call must be made with thread 0");

      /* Allocate memory for size and pointer */

      loc0 = ReallocMem(DATA_ARRAY, 2);
      WDB[root] = (double)loc0;

      /* Put size */

      WDB[loc0++] = (double)np;
      
      /* Check type */
      
      if (type == DATA_ARRAY)
        {
          /* Allocate memory for array */

          ptr = ReallocMem(DATA_ARRAY, np);
          WDB[loc0] = (double)ptr;

          /* Get pointer to data */

          dat = &WDB[ptr];
        }
      else if (type == PRIVA_ARRAY)
        {
          /* Allocate memory for array */

          ptr = AllocPrivateData(np, PRIVA_ARRAY);
          WDB[loc0] = (double)ptr;

          /* Get pointer to data (id must be zero) */

          dat = &PRIVA[ptr];
        }
      else
        Die(FUNCTION_NAME, "Invalid array type");

      /* Put null pointer to make sure that data is pre-allocated */
      /* before used */

      dat = NULL;
    }
  
  /* Return pointer to array */
  
  return dat;
}

/*****************************************************************************/
