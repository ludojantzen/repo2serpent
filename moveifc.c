/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : moveifc.c                                      */
/*                                                                           */
/* Created:       2016/06/08 (VVa)                                           */
/* Last modified: 2016/09/16 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Moves interface TMP and DF from EOI list to BOI list         */
/*                                                                           */
/* Comments: -Actually just swaps the pointers                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MoveIFC:"

/*****************************************************************************/

void MoveIFC()
{
  long loc0, ptr, fpe, axi, ang, fin;

  /* Return if no interfaces are defined */

  if ((loc0 = (long)RDB[DATA_PTR_IFC0]) < VALID_PTR)
    return;

  /* Loop over interfaces */

  while (loc0 > VALID_PTR)
    {

      if ((RDB[loc0 + IFC_TYPE] != IFC_TYPE_FUEP) &&
          (RDB[loc0 + IFC_TYPE] != IFC_TYPE_FPIP))
        {
          /* Non-fuel behavior interfaces */

          /* Move density factor list */

          if ((ptr = (long)RDB[loc0 + IFC_PTR_DF_LIST_BOI]) > VALID_PTR)
            {
              /* Old EOI is new BOI */

              WDB[loc0 + IFC_PTR_DF_LIST_BOI] = RDB[loc0 + IFC_PTR_DF_LIST];

              /* Old BOI is now reused as EOI */

              WDB[loc0 + IFC_PTR_DF_LIST] = (double)ptr;

            }

          /* Move temperature list */

          if ((ptr = (long)RDB[loc0 + IFC_PTR_TMP_LIST_BOI]) > VALID_PTR)
            {
              /* Old EOI is new BOI */

              WDB[loc0 + IFC_PTR_TMP_LIST_BOI] = RDB[loc0 + IFC_PTR_TMP_LIST];

              /* Old BOI is now reused as EOI */

              WDB[loc0 + IFC_PTR_TMP_LIST] = (double)ptr;

            }
        }
      else
        {
          /* Fuel behavior interfaces */

          /* Loop over rod definitions */

          fpe = (long)RDB[loc0 + IFC_PTR_FUEP];

          while (fpe > VALID_PTR)
            {

              /* Loop over axial zones */

              axi = (long)RDB[fpe + IFC_FUEP_PTR_AX];

              while (axi > VALID_PTR)
                {

                  /* Loop over angular zones */

                  ang = (long)RDB[axi + IFC_FUEP_AX_PTR_ANG];

                  while (ang > VALID_PTR)
                    {

                      /* Move temperature list */

                      if ((ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_TMP_BOI]) >
                          VALID_PTR)
                        {
                          /* Old EOI is new BOI */

                          WDB[ang + IFC_FUEP_ANG_PTR_TMP_BOI] =
                            RDB[ang + IFC_FUEP_ANG_PTR_TMP];

                          /* Old BOI is now reused as EOI */

                          WDB[ang + IFC_FUEP_ANG_PTR_TMP] = (double)ptr;

                        }

                      /* Move DF list */

                      if ((ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_DF_BOI]) >
                          VALID_PTR)
                        {
                          /* Old EOI is new BOI */

                          WDB[ang + IFC_FUEP_ANG_PTR_DF_BOI] =
                            RDB[ang + IFC_FUEP_ANG_PTR_DF];

                          /* Old BOI is now reused as EOI */

                          WDB[ang + IFC_FUEP_ANG_PTR_DF] = (double)ptr;

                        }

                      /* Move HOT_R2 list */

                      if ((ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2_BOI]) >
                          VALID_PTR)
                        {
                          /* Old EOI is new BOI */

                          WDB[ang + IFC_FUEP_ANG_PTR_HOT_R2_BOI] =
                            RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2];

                          /* Old BOI is now reused as EOI */

                          WDB[ang + IFC_FUEP_ANG_PTR_HOT_R2] = (double)ptr;

                        }

                      /* Next angular zone */

                      ang = NextItem(ang);
                    }

                  /* Next axial zone */

                  axi = NextItem(axi);
                }

              /***********************************/
              /* Move FINIX results if available */
              /***********************************/

              if ((fin = (long)RDB[fpe + IFC_FUEP_PTR_FINIX]) > VALID_PTR)
                if ((ptr = (long)RDB[fin + FINIX_PTR_RESULTS_BOI]) > VALID_PTR)
                  {
                    /* Old EOI is new BOI */

                    WDB[fin + FINIX_PTR_RESULTS_BOI] = RDB[fin + FINIX_PTR_RESULTS];

                    /* Old BOI is now reused as EOI */

                    WDB[fin + FINIX_PTR_RESULTS] = (double)ptr;

                  }


              /* Next rod */

              fpe = NextItem(fpe);
            }
        }

      /* Next interface */

      loc0 = NextItem(loc0);
    }

}
