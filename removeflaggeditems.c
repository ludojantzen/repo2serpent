/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : removeflaggeditems.c                           */
/*                                                                           */
/* Created:       2010/09/22 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: - Removes flagged items from list                            */
/*                                                                           */
/* Comments: - Returns the number of removed items                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RemoveFlaggedItems:"

/*****************************************************************************/

long RemoveFlaggedItems(long lst, long param, long flag, long mode)
{
  long loc0, ptr, count;

  /* Check pointer */
  
  if (lst < 1)
    return 0;    
  
  /* Check if list is closed */

  if ((long)RDB[lst + LIST_PTR_DIRECT] > VALID_PTR)
    Die(FUNCTION_NAME, "Trying to remove items from a closed list");

  /* Reset count */

  count = 0;

  /* Loop over list */
  
  loc0 = FirstItem(lst);
  
  while (loc0 > VALID_PTR)
    {
      /* Copy pointer */

      ptr = loc0;

      /* Pointer to next */
      
      loc0 = NextItem(loc0);

      /* Check mode */

      if (mode == YES)
	{
	  /* Remove if flag is set */

	  if ((long)RDB[ptr + param] & flag)
	    {
	      RemoveItem(ptr);
	      count++;
	    }
	}
      else if (mode == NO)
	{
	  /* Remove if flag is not set */

	  if (!((long)RDB[ptr + param] & flag))
	    {
	      RemoveItem(ptr);
	      count++;
	    }
	}
      else
	Die(FUNCTION_NAME, "Invalid mode");
    }

  /* Return count */

  return count;
}

/*****************************************************************************/
