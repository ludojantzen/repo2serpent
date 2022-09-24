/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reprocess.c                                    */
/*                                                                           */
/* Created:       2012/07/09 (JLe)                                           */
/* Last modified: 2020/06/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Main reprocessing routine                                    */
/*                                                                           */
/* Comments: - Used only for testing purposes at the moment                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Reprocess:"

/*****************************************************************************/

void Reprocess(long dep)
{
  long rep, loc0, uni1, uni2;

  /* Get pointer to reprocessor */

  if ((rep = (long)RDB[dep + DEP_HIS_PTR_REPROC]) < VALID_PTR)
    return;

  /* This is disabled for now to test continuous reprocessing */

#ifdef DEBUG

  Warn(FUNCTION_NAME, "Function call disabled");

#endif

  return;

  fprintf(outp, "\nReprocessing fuel:\n\n");

  /***************************************************************************/

  /***** Universe swaps ******************************************************/

  /* Loop over swaps */

  loc0 = (long)RDB[rep + REPROC_PTR_SWAP_LIST];
  while (loc0 > VALID_PTR)
    {
      /* Get pointer to first universe */

      uni1 = (long)RDB[loc0 + REPROC_SWAP_PTR_UNI1];
      CheckPointer(FUNCTION_NAME, "(uni1)", DATA_ARRAY, uni1);

      /* Get pointer to second universe */

      uni2 = (long)RDB[loc0 + REPROC_SWAP_PTR_UNI2];
      CheckPointer(FUNCTION_NAME, "(uni2)", DATA_ARRAY, uni2);

      fprintf(outp, "Swapping universes %s and %s...\n",
              GetText(uni1 + UNIVERSE_PTR_NAME),
              GetText(uni2 + UNIVERSE_PTR_NAME));

      /* Swap */

      SwapUniverses(uni1, uni2);

      /* Next swap */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Re-create geometry **************************************************/

  /* Count number of cells */

  CellCount(-1, -1, 0, 1);

  /* Calculate material volumes */

  MaterialVolumes();

  /* Copy compositions to parent materials */

  SumDivCompositions();

  /***************************************************************************/

  /* Print newline */

  fprintf(outp, "\n");

  /* Plot */

  GeometryPlotter(YES);

  /*
  VolumesMC();
  */
}

/*****************************************************************************/
