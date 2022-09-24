/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dataifcadens.c                                 */
/*                                                                           */
/* Created:       2018/06/05 (VVa)                                           */
/* Last modified: 2018/06/21 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns nuclide atomic density from data interface           */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DataIFCAdens:"

/*****************************************************************************/

double DataIFCAdens(long mat, long nuc, long id)
{
  long msh, arr, loc0, ptr, uni, ncol, idx;
  double x, y, z, adens;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Initialize some values */

  adens = 0.0;

  /* Check that the material has some data interfaces defined */

  if ((arr = (long)RDB[mat + MATERIAL_PTR_DATAIFC_ARR]) < VALID_PTR)
    {
      Die(FUNCTION_NAME, "Trying to find atomic density for %s in material %s "
          "from a data interface, but material has no data interfaces linked.",
          GetText(nuc + NUCLIDE_PTR_NAME), GetText(mat + MATERIAL_PTR_NAME));

      return adens;
    }

  /* Loop over material's data interfaces to find the one giving nuclide density */
  /* for "nuc" */

  while ((loc0 = (long)RDB[arr]) > VALID_PTR)
    {
      /* Check nuclide for current data interface */

      if (((long)RDB[loc0 + DATAIFC_PTR_NUCLIDE] == nuc)
          && ((long)RDB[loc0 + DATAIFC_DATA_TYPE] == DATA_TYPE_IFC_ADENS))
        break;

      /* Next data interface */

      arr++;
    }

  /* Check to make sure that we found the data interface */

  if ((loc0 = (long)RDB[arr]) < VALID_PTR)
    {
      Die(FUNCTION_NAME, "Could not find data interface corresponding to atomic density "
          "of %s for material %s", GetText(nuc + NUCLIDE_PTR_NAME), GetText(mat + MATERIAL_PTR_NAME));

      return adens;
    }

  /* Get mesh from the data interface*/

  msh = (long)RDB[loc0 + DATAIFC_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get current collision count */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Check if mesh index has already been stored for this collision */

  if ((idx = (long)TestValuePair(msh + MESH_PREV_COL_IDX, ncol, id)) < 0.0)
    {
      /************************************/
      /* Need to find the mesh index here */
      /************************************/

      /* Get the correct coordinates */

      if (RDB[msh + MESH_LOCAL_COORDS] == (double)YES)
        {
          /* Get local coordinates from collision universe */

          ptr = (long)RDB[DATA_PTR_COLLISION_UNI];
          uni = (long)GetPrivateData(ptr, id);
        }
      else
        {
          /* Get local coordinates from root universe */

          uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
        }

      /* Check universe pointer */

      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Get coordinates */

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_X];
      CheckPointer(FUNCTION_NAME, "(xptr)", PRIVA_ARRAY, ptr);
      x = GetPrivateData(ptr, id);

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Y];
      CheckPointer(FUNCTION_NAME, "(yptr)", PRIVA_ARRAY, ptr);
      y = GetPrivateData(ptr, id);

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Z];
      CheckPointer(FUNCTION_NAME, "(zptr)", PRIVA_ARRAY, ptr);
      z = GetPrivateData(ptr, id);

      /* Get mesh bin index */

      if ((idx = MeshIndex(msh, x, y, z, -1.0)) < 0)
        Die(FUNCTION_NAME, "Could not find mesh index for mesh %s "
            " linked to data interface %s at coordinates (%E, %E, %E)",
            GetText(msh + MESH_PTR_NAME), GetText(loc0 + DATAIFC_PTR_INPUT_FNAME),
            x, y, z);

      /* Store the value pair so that we won't need to seek through the same mesh again */

      StoreValuePair(msh + MESH_PREV_COL_IDX, (double)ncol, (double)idx, id);
    }

  /* Now we have the mesh index at hand, get pointer to data ifc data array */

  ptr = (long)RDB[loc0 + DATAIFC_PTR_DATA];

  /* Get atomic density from array */

  adens = RDB[ptr + idx];

  /* Get mesh index corresponding to collision coordinates */

  return adens;
}

/*****************************************************************************/
