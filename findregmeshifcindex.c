/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findregmeshifcindex.c                          */
/*                                                                           */
/* Created:       2018/11/15 (VVa)                                           */
/* Last modified: 2018/11/15 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Translates coordinates to an index in the "cell" list of a   */
/*              regular mesh based interface.                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindRegMeshIFCIndex:"

/*****************************************************************************/

long FindRegMeshIFCIndex(long ifc, double x, double y, double z, long id)
{
  long type, msh, nmesh, lat, dumlong, i, n, idx;
  double xx, yy, zz;

  /* Get interface type */

  type = (long)RDB[ifc + IFC_TYPE];

  /* Get pointer to mesh */

  msh = (long)RDB[ifc + IFC_PTR_SEARCH_MESH_LIST];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Figure out total index */

  if (type != IFC_TYPE_REG_MESH_MULTILVL)
    i = MeshIndex(msh, x, y, z, -1.0);
  else
    {
      /*** Multi-level mesh ***/

      /* Get number of levels */

      nmesh = (long)RDB[ifc + IFC_N_MESH_LVL];

      /* Copy coordinates to working variables */

      xx = x;
      yy = y;
      zz = z;

      /* Loop over all levels to calculate index */

      lat = (long)RDB[ifc + IFC_PTR_LATTICE_LIST];
      i = 0;

      for (n = 0; n < nmesh; n++)
        {
          /* Find index for this lattice */

          idx = FindLatticeRegion(lat, -1, &xx, &yy, &zz,
                                  &dumlong, id);

          /* If we got a negative (not-found) index, break */

          if (idx < 0)
            {
              i = -1;
              break;
            }

          /* Increment i */

          i = i + idx;


          /* Get next level lattice */

          lat = NextItem(lat);

          /* If no next level lattice, break */

          if (lat < VALID_PTR)
            break;

          /* Otherwise, increment i to the next level/basis */

          i = i*(long)RDB[lat + LAT_NTOT];
        }

      /* Check that index is below number of mesh cells */

      if (i >= (long)RDB[ifc + IFC_NP])
        Die(FUNCTION_NAME, "WTF: %ld vs. %ld (i/NP)", i, (long)RDB[ifc + IFC_NP]);
    }

  return i;
}
