/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : preallocmem.c                                  */
/*                                                                           */
/* Created:       2012/05/29 (JLe)                                           */
/* Last modified: 2018/03/14 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Pre-allocates memory from data blocks                        */
/*                                                                           */
/* Comments: - Used for allocating larger blocks at once to speed up the     */
/*             calculation.                                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PreallocMem:"

/*****************************************************************************/

void PreallocMem(long sz, long type)
{
  long ptr;

  /* Check array type */

  switch (type)
    {
    case DATA_ARRAY:
      {
        /* Allocate memory */

        ptr = ReallocMem(DATA_ARRAY, sz);

        /* Put pointer back to beginning of allocated region */

        WDB[DATA_ALLOC_MAIN_SIZE] = (double)ptr;  

        /* Break case */
        
        break;
      }
    case RES1_ARRAY:
      {
        /* Allocate memory */

        ptr = ReallocMem(RES1_ARRAY, sz);

        /* Put pointer back to beginning of allocated region */

        WDB[DATA_ALLOC_RES1_SIZE] = (double)ptr;  

        /* Break case */

        break;
      }
    case ACE_ARRAY:
      {
        /* Allocate memory */

        ptr = ReallocMem(ACE_ARRAY, sz);

        /* Put pointer back to beginning of allocated region */

        WDB[DATA_ALLOC_ACE_SIZE] = (double)ptr;  

        /* Break case */

        break;
      }
    case PRIVA_ARRAY:
      {
        /* Allocate memory */

        ptr = AllocPrivateData(sz, PRIVA_ARRAY);

        /* Put pointer back to beginning of allocated region */

        WDB[DATA_ALLOC_PRIVA_SIZE] = (double)ptr;  

        /* Break case */

        break;
      }
    case BUF_ARRAY: 
      {
        /* Allocate memory */
        
        ptr = AllocPrivateData(sz, BUF_ARRAY);

        /* Put pointer back to beginning of allocated region */

        WDB[DATA_ALLOC_BUF_SIZE] = (double)ptr;  

        /* Break case */

        break;
      }
    case RES2_ARRAY:
      {
        /* Allocate memory */

        ptr = AllocPrivateData(sz, RES2_ARRAY);

        /* Put pointer back to beginning of allocated region */

        WDB[DATA_ALLOC_RES2_SIZE] = (double)ptr;  

        /* Break case */

        break;
      }
    case RES3_ARRAY:
      {
        /* Allocate memory */

        ptr = AllocPrivateData(sz, RES3_ARRAY);

        /* Put pointer back to beginning of allocated region */

        WDB[DATA_ALLOC_RES3_SIZE] = (double)ptr;  

        /* Break case */

        break;
      }
    default:
      {
        /* Error */

        Die(FUNCTION_NAME, "Invalid data type");
      }
    }

  /* Calculate memory size in bytes */

  CalculateBytes();
}

/*****************************************************************************/
