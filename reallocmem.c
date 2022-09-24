/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reallocmem.c                                   */
/*                                                                           */
/* Created:       2010/11/19 (JLe)                                           */
/* Last modified: 2015/05/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Allocates memory from WDB, ACE and RES1 arrays               */
/*                                                                           */
/* Comments: - RES2, PRIVA and BUF are handled in AllocPrivateData()         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReallocMem:"

/*****************************************************************************/

long ReallocMem(long type, long sz)
{
  long memsize, realsize, block, loc0, totsize;

  /* Check allow flag */

  if (RDB != NULL)
    if ((long)RDB[DATA_ALLOW_MEM_OP] == NO)
      Die(FUNCTION_NAME, "Memory allocation not allowed");

  /* Check size */

  if (sz < 1)
    Die (FUNCTION_NAME, "Invalid data size %ld (type %ld)", sz, type);
    
  /* Avoid compiler warning */

  realsize = 0;
  memsize = 0;

  /* Allocate memory in 20MB blocks */
  
  block = (long)(20.0*MEGA/sizeof(double));

  /***************************************************************************/
  
  /***** Get size of datablock ***********************************************/
  
  if (type == DATA_ARRAY)
    {
      if (WDB != NULL)
        {
          memsize = (long)RDB[DATA_ALLOC_MAIN_SIZE];
          realsize = (long)RDB[DATA_REAL_MAIN_SIZE];
        }
      else
        {
          memsize = 0;
          realsize = 0;
        }
    }
  else if (type == RES1_ARRAY)
    {
      if (RES1 != NULL)
        {
          memsize = (long)RDB[DATA_ALLOC_RES1_SIZE];
          realsize = (long)RDB[DATA_REAL_RES1_SIZE];
        }
      else
        {
          memsize = 0;
          realsize = 0;
        }

      /* Multiply size by stat variable size */

      sz = sz*STAT_BLOCK_SIZE;
    }
  else if (type == ACE_ARRAY)
    {
      if (ACE != NULL)
        {
          memsize = (long)RDB[DATA_ALLOC_ACE_SIZE];
          realsize = (long)RDB[DATA_REAL_ACE_SIZE];
        }
      else
        {
          memsize = 0;
          realsize = 0;
        }
    }
  else
    Die(FUNCTION_NAME, "Invalid data type %ld", type);

  /* Remember pointer */

  loc0 = memsize;

  /* Increase size */
  
  memsize = memsize + sz;

  /***************************************************************************/
  
  /***** Check limit *********************************************************/

  /* Check pointer */

  if (WDB != NULL)
    { 
      /* Calculate total data size */
      
      CalculateBytes();
      
      /* Get total size */
      
      totsize = (long)RDB[DATA_TOTAL_BYTES];

      if (memsize > realsize)
        totsize = totsize + block*sizeof(double);
      
      /* Compare to limit */

      if ((RDB[DATA_CPU_MEM] > 0.0) && 
          ((double)totsize/GIGA > RDB[DATA_CPU_MEM_FRAC]*RDB[DATA_CPU_MEM]))
        Error(0, "Out of memory error (allocating %1.1fGb, limit %1.1fGb = %1.2f * %1.1fGb)", 
              (double)totsize/GIGA, RDB[DATA_CPU_MEM_FRAC]*RDB[DATA_CPU_MEM],
              RDB[DATA_CPU_MEM_FRAC], RDB[DATA_CPU_MEM]);
    }

  /***************************************************************************/
  
  /***** Allocate memory for main data array *********************************/

  if (type == DATA_ARRAY)
    {
      /* Compare to real data size */

      if (memsize > realsize)
        {
          /* Increase real memory size */

          realsize = memsize + block;

          /* Allocate more memory */
          
          WDB = (double *)Mem(MEM_REALLOC, WDB, realsize*sizeof(double));
          
          /* Initialize values */
          
          memset(&WDB[loc0], 0.0, (sz + block)*sizeof(double));

          /* Update real data size */

          WDB[DATA_REAL_MAIN_SIZE] = (double)realsize;
        }

      /* Update block size */
      
      WDB[DATA_ALLOC_MAIN_SIZE] = (double)memsize;
    }

  /***************************************************************************/
  
  /***** Allocate memory for RES1 array **************************************/

  else if (type == RES1_ARRAY)
    {
      /* Compare to real data size */

      if (memsize > realsize)
        {
          /* Increase real memory size */

          realsize = memsize + block;

          /* Allocate more memory */
          
          RES1 = (double *)Mem(MEM_REALLOC, RES1,realsize*sizeof(double));

          /* Initialize values */
          
          memset(&RES1[loc0], 0.0, (sz + block)*sizeof(double));

          /* Update real data size */

          WDB[DATA_REAL_RES1_SIZE] = (double)realsize;
        }
      
      /* Update block size */
      
      WDB[DATA_ALLOC_RES1_SIZE] = (double)memsize;
    }

  /***************************************************************************/

  /***** Allocate memory for ACE array ***************************************/

  else if (type == ACE_ARRAY)
    {
      /* Compare to real data size */

      if (memsize > realsize)
        {
          /* Increase real memory size */

          realsize = memsize + block;

          /* Allocate more memory */
          
          ACE = (double *)Mem(MEM_REALLOC, ACE, realsize*sizeof(double));

          /* Initialize values */
          
          memset(&ACE[loc0], 0.0, (sz + block)*sizeof(double));

          /* Update real data size */
          
          WDB[DATA_REAL_ACE_SIZE] = (double)realsize;
        }
      
      /* Update block size */
      
      WDB[DATA_ALLOC_ACE_SIZE] = (double)memsize;
    }

  /***************************************************************************/

  /* Put read-only pointers */

  RDB = (const double *)WDB;

  /* Calculate allocated memory in bytes */

  CalculateBytes();
  
  /* Return pointer to allocated memory */

  return loc0;
}

/*****************************************************************************/
