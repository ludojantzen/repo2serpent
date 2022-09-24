/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : redistributeques.c                             */
/*                                                                           */
/* Created:       2015/09/29 (VVa)                                           */
/* Last modified: 2015/09/29 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Redistributes queued particles between OpenMp threads        */
/*                                                                           */
/* Comments:   -Used only in dynamic simulation mode                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReDistributeQues:"

/*****************************************************************************/

long ReDistributeQues()
{
  long n, tomove, ave, loc0, loc1, ptr, id, id0, id1, min, max, sz, tot, nt;

  /* Get number of OpenMp threads */

  nt = (long)RDB[DATA_OMP_MAX_THREADS];

  /* Reset total counts */

  tot = 0;

  /* There is only one type of que */

  loc0 = DATA_PART_PTR_QUE;

  /* Check pointer */

  if ((long)RDB[loc0] < VALID_PTR)
    Die(FUNCTION_NAME, "Pointer error");
  /*
  printf("Que sizes:\n");
  */
  /* Loop over all threads and calculate total particle count */
  
  tot = 0;

  for (id = 0; id < nt; id++)
    {
      /* Get list size */
              
      loc1 = (long)RDB[OMPPtr(loc0, id)];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
      /*
      printf("%ld ", ListSize(loc1) - 1);
      */
      /* Add to total */

      tot += ListSize(loc1) - 1;
    }
  /*
  printf("\n");
  */
  /* Calculate average */
  /* This is long, so ques should have ave or ave + 1 particles */

  ave = tot / nt;

  /* Loop and redistribute */

  id0 = -1;
  id1 = -1;

  min = 100000000000;
  max = -1;
          
  while (1 != 2)
    {

      /* Find stack with lowest and highest number of particles */
          
          
      for (id = 0; id < nt; id++)
        {
          /* Get list size */
              
          loc1 = (long)RDB[OMPPtr(loc0, id)];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
          sz = ListSize(loc1) - 1;
              
          /* Compare to minimum */
              
          if ((sz < ave)  && (id0 < 0))
            {
              min = sz;
              id0 = id;
            }
              
          /* Compare to maximum */
              
          if ((sz > ave + 1) && (id1 < 0))
            {
              max = sz;
              id1 = id;
            }

          /* Check if we know que with less than required amount of particles */
          /* and a second que with more than required amount of particles */

          if ((id0 != -1) && (id1 != -1))
            break;
        }

      if ((id0 < 0) || (id1 < 0))
        break;

      /* Move particles from id1 to id0 until */
      /* one of them has ave or ave + 1 particles. */
          
      /* Move half of particles to other stack */

      if ((ave - min) < (max - ave))
        {
          tomove = ave - min;

          min = 100000000000;
          max = max - tomove;
        }
      else
        {
          tomove = max - ave;          

          max = 0;
          min = min + tomove;
        }
              
      for (n = 0; n < tomove; n++)
        {
          ptr = FromQue(id1);
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          ToQue(ptr, id0);          
        }

      /* Reset thread id for the que that now has the required number of */
      /* particles */

      if (min == 100000000000)
        id0 = -1;
      else if (max == 0)
        id1 = -1;
      else
        Die(FUNCTION_NAME, "Forgot to reset que size?");

    }

  /* Loop over ques to calculate new total size*/

  sz = 0;

  for (id = 0; id < nt; id++)
    {
      /* Get list size */
              
      loc1 = (long)RDB[OMPPtr(loc0, id)];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
      /*
      printf("%ld ", ListSize(loc1) - 1);
      */
      /* Add to total */

      sz += ListSize(loc1) - 1;
    }
  /*
  printf("\n\n");
  */
  if (sz != tot)
    Die(FUNCTION_NAME, "Lost %ld particles", tot - sz);

  return tot;

}

/*****************************************************************************/
