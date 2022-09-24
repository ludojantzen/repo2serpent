/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : createfinixifc.c                                */
/*                                                                           */
/* Created:       2013/03/11 (VVa)                                           */
/* Last modified: 2019/05/15 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Creates an IFC-template for the ReadInterface()              */
/*              for the initialization of the FINIX interface                */
/*                                                                           */
/* Comments:   -Tän vois periaatteessa kirjoittaa suoraan muistiinkin        */
/*              mutta lienee järkevintä luoda kaikki rajapinnat vasta        */
/*              ReadInterface():ssa                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#define FUNCTION_NAME "CreateFinixIFC:"

/*****************************************************************************/

void CreateFinixIFC()
{
  long fib, loc0, loc1, found, nu, ptr, i;
  char tmpstr[MAX_STR];
  double mem0, mem1;

  /* Get memory size before creation of the interface */

  mem0 = RDB[DATA_ALLOC_MAIN_SIZE];

  /* Create new interface */

  loc0 = NewItem(DATA_PTR_IFC0, IFC_BLOCK_SIZE);

  /* Print filename to string */

  sprintf(tmpstr,"./Finix.ifc");

  /* Put file name and some PARAM_N_COMMON things*/

  WDB[loc0 + IFC_PTR_INPUT_FNAME] = (double)PutText(tmpstr);
  WDB[loc0 + PARAM_PTR_NAME] = (double)PutText("ifc");
  WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(tmpstr);
  WDB[loc0 + PARAM_LINE] = 0;

  /* Create interface */

  ReadInterface(loc0, 0);

  /* Get memory size after creation of the interface */

  mem1 = RDB[DATA_ALLOC_MAIN_SIZE];

  /* Store ifc-memory size to interface (for MPI transfer) */

  WDB[loc0 + IFC_MEM_SIZE] = mem1 - mem0;

  /* Put pointer to interface to pins */

  fib = (long)RDB[DATA_PTR_FIN0];

  while (fib > VALID_PTR)
    {

      /* Get pointer to fuel rod block */

      loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Locate correct interface pin */

      while (loc1 > VALID_PTR)
        {
          /* Reset found flag */

          found = 0;

          /* Get number of pins in this interface pin     */
          /* might be larger than 1 due to axial segments */

          nu = RDB[loc1 + IFC_FUEP_N_UNI];

          /* Get pointer to segment list */

          ptr = (long)RDB[loc1 + IFC_FUEP_PTR_UNI_LIST];

          /* If one of the segments corresponds to this FINIX pin */
          /* Activate found flag */

          for (i=0; i < nu; i++)
            {
              if (CompareStr(ptr + i, fib + FINIX_PTR_UNI_NAME))
                found=1;
            }

          if (found==1)
            break;
          else
            loc1 = NextItem(loc1);

        }

      /* Check if found */

      if(loc1 < VALID_PTR)
        Die(FUNCTION_NAME, "Could not find universe");

      /* Add pointer to fuep */

      WDB[fib + FINIX_PTR_FUEP] = (double)loc1;

      /* Check that FINIX is not already linked */

      if ((long)RDB[loc1 + IFC_FUEP_PTR_FINIX] > VALID_PTR)
        Warn(FUNCTION_NAME,
             "Fuel behavior interface rod %s already linked to a FINIX rod",
             GetText(loc1 + IFC_FUEP_PTR_UNI_LIST));

      /* Add pointer to FINIX */

      WDB[loc1 + IFC_FUEP_PTR_FINIX] = (double)fib;

      /* Add pointer to interface filename */

      WDB[fib + FINIX_PTR_IFC_FNAME] = (double)PutText(tmpstr);

      /* Add pointer to interface*/

      WDB[fib + FINIX_PTR_IFC] = (double)loc0;

      fib = NextItem(fib);
    }

}
#endif
/*****************************************************************************/
