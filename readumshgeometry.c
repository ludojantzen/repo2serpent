/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readumshgeometry.c                             */
/*                                                                           */
/* Created:       2013/11/23 (JLe)                                           */
/* Last modified: 2019/10/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Reads unstructured mesh based geometry                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadUMSHGeometry:"

/*****************************************************************************/

void ReadUMSHGeometry()
{
  long umsh, loc0, ptr;
#ifdef printoutumsh
  long mat, loc2, loc3;
#endif

  /* Check pointer */

  if ((umsh = (long)RDB[DATA_PTR_UMSH0]) < VALID_PTR)
    return;

  fprintf(outp, "Reading unstructured mesh based geometries:\n\n");

  /* Loop over definitions */

  while (umsh > VALID_PTR)
    {
      if ((long)RDB[umsh + UMSH_PTR_FNAME] > 0)
        {
          /* Will be read in readifcofmesh.c */

          umsh = NextItem(umsh);

          continue;
        }

      /* Pointer to interface structure */

      loc0 = (long)RDB[umsh + UMSH_PTR_IFC];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Read mesh structure/data */

      ReadOFMesh(loc0);

      /* Allocate memory for previous cell pointer */

      ptr = AllocPrivateData(1, PRIVA_ARRAY);
      WDB[loc0 + IFC_PTR_PREV_CELL] = (double)ptr;

      /* Allocate memory for previous collison */

      AllocValuePair(loc0 + IFC_PTR_PREV_COL_CELL);

      /* Next geometry */

      umsh = NextItem(umsh);

      /***********************************************************************/
    }

  /* Done */

  fprintf(outp, "OK.\n\n");

}

/*****************************************************************************/
