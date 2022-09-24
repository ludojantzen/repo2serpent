/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : storermxevent.c                                */
/*                                                                           */
/* Created:       2019/03/10 (JLe)                                           */
/* Last modified: 2019/03/10 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Handles buffer-wise scoring of RMX events.                   */
/*                                                                           */
/* Comments: - Used with MCNP-type weight window generator (mode 4).         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StoreRMXEvent:"

/*****************************************************************************/

long StoreRMXEvent(long rmx, long ptr0, long ptr1, double val, long idx, 
                   long id)
{
  long n, loc0;
  double v0;

  /* Check if value is response */

  if (ptr0 < 0)
    {
      /***********************************************************************/

      /***** Score buffered contributions ************************************/

      /* Pointer to buffer */

      loc0 = (long)RDB[rmx + RMX_PTR_EVENT_BUFF];
      CheckPointer(FUNCTION_NAME, "(loc0)", PRIVA_ARRAY, loc0);

      /* Loop over buffer */

      for (n = idx - 1; n > -1; n--)
        {
          /* Get values */

          ptr0 = (long)GetPrivateData(loc0 + 2*n, id);
          v0 = GetPrivateData(loc0 + 2*n + 1, id);

          /* Check */

          if (ptr0 < VALID_PTR)
            break;

          /* Reset */
          
          PutPrivateData(loc0 + 2*n, -1.0, id);
          PutPrivateData(loc0 + 2*n + 1, 0.0, id);
          
          /* Score */
          
          AddPrivateRes(ptr0, val*v0, id);
        }

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Store values in buffer ******************************************/

      /* Check pointers */

      CheckPointer(FUNCTION_NAME, "(ptr0)", RES2_ARRAY, ptr0);
      CheckPointer(FUNCTION_NAME, "(ptr1)", RES2_ARRAY, ptr1);

      /* Put total */

      AddPrivateRes(ptr0, val, id);
      
      /* Check index */
      
      if (idx > MAX_RMX_BUFF - 1)
        Die(FUNCTION_NAME, "RMX Buffer full");
      
      /* Pointer to buffer */
      
      loc0 = (long)RDB[rmx + RMX_PTR_EVENT_BUFF];
      CheckPointer(FUNCTION_NAME, "(loc0)", PRIVA_ARRAY, loc0);
      
      /* Put values in buffer */
      
      PutPrivateData(loc0 + 2*idx, (double)ptr1, id);
      PutPrivateData(loc0 + 2*idx + 1, val, id);
  
      /* Add to index */

      idx = idx + 1;

      /***********************************************************************/
    }

  /* Return index */

  return idx;
}

/*****************************************************************************/
