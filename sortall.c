/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sortall.c                                      */
/*                                                                           */
/* Created:       2011/03/01 (JLe)                                           */
/* Last modified: 2019/04/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sorts reaction, surface etc. lists to speed up calculation   */
/*                                                                           */
/* Comments: - Reaction list sorting doesn't work in reproduceable MPI mode  */
/*           - Sorting the tet mesh lists is expensive if the search         */
/*             mesh is large (not in use at the moment).                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SortAll:"

/*****************************************************************************/

void SortAll()
{
  long mat, cell, nuc, lst, idx, uni, icm, loc0, msh, nx, ny, nz, i, j, k, ptr;
  long sz, surf;

  /* Get cycle index */

  idx = (long)RDB[DATA_CYCLE_IDX];

  /* Check condition (NOTE: estimateruntime.c olettaa että */
  /* sorttausväli on tuo.) */

  if (idx != (long)RDB[DATA_SORT_COUNT])
    return;
  else
    WDB[DATA_SORT_COUNT] = RDB[DATA_SORT_COUNT]*3.0;

  /***************************************************************************/

  /***** Reaction lists in nuclides and materials ****************************/

  /* Check MPI option */

  if ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO)
    {
      /* Loop over nuclides */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Check mode */

          if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
            {
              /* Pointer to partial reaction list */

              lst = (long)RDB[nuc + NUCLIDE_PTR_SAMPLE_REA_LIST];
              CheckPointer(FUNCTION_NAME, "(lst1)", DATA_ARRAY, lst);

              /* Sort list */

              SortList(lst, RLS_DATA_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
            }

          /* Next nuclide */

          nuc = NextItem(nuc);
        }

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while(mat > VALID_PTR)
        {
          /* Sort material-wise list of totals */

          if ((lst = (long)RDB[mat + MATERIAL_PTR_TOT_REA_LIST])
              > VALID_PTR)
            SortList(lst, RLS_DATA_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);

          if ((lst = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANT_LIST])
              > VALID_PTR)
            SortList(lst, RLS_DATA_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);

          /* Next material */

          mat = NextItem(mat);
        }
    }

  /***************************************************************************/

  /***** Intersection lists for cells ****************************************/

  /* Pointer to cell list */

  loc0 = (long)RDB[DATA_PTR_C0];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get list size */

  sz = ListSize(loc0);

  /* Loop over list */

#ifdef OPEN_MP
#pragma omp parallel for private(cell, lst)
#endif

  for (i = 0; i < sz; i++)
    {
      /* Get pointer */

      cell = ListPtr(loc0, i);
      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

      /* Pointer to intersection list */

      if ((lst = (long)RDB[cell + CELL_PTR_SURF_INSC]) > VALID_PTR)
        {
          /* Sort list */

          SortList(lst, CELL_INSC_PTR_OUT_COUNT, SORT_MODE_DESCEND_PRIVA);
        }
    }

  /***************************************************************************/

  /***** Cell lists for universes ********************************************/

  /* Pointer to universe list */

  loc0 = (long)RDB[DATA_PTR_U0];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get list size */

  sz = ListSize(loc0);

  /* Loop over list */

#ifdef OPEN_MP
#pragma omp parallel for private(uni, lst)
#endif

  for (i = 0; i < sz; i++)
    {
      /* Get pointer */

      uni = ListPtr(loc0, i);
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Pointer to cell list */

      if ((lst = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST]) > VALID_PTR)
        {
          /* Sort list */

          SortList(lst, CELL_LIST_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
        }

      /* Pointer to source cell list */

      if ((lst = (long)RDB[uni + UNIVERSE_PTR_SRC_CELL_LIST]) > VALID_PTR)
        {
          /* Sort list */

          SortList(lst, CELL_LIST_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
        }
    }

  /* Loop over universes */

  uni = (long)RDB[DATA_PTR_U0];
  while (uni > VALID_PTR)
    {
      /* Pointer to cell list */

      if ((lst = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST]) > VALID_PTR)
        {
          /* Pointer to list */

          loc0 = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          /* Get list size */

          sz = ListSize(loc0);

          /* Loop over list */

#ifdef OPEN_MP
#pragma omp parallel for private(lst, cell, ptr)
#endif

          for (i = 0; i < sz; i++)
            {
              /* Get pointer */

              lst = ListPtr(loc0, i);
              CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

              /* Pointer to cell */

              cell = (long)RDB[lst + CELL_LIST_PTR_CELL];
              CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

              /* Sort list */

              if ((ptr = (long)RDB[cell + CELL_PTR_SEARCH_LIST]) > VALID_PTR)
                SortList(ptr, CELL_LIST_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
            }
        }

      /* Next universe */

      uni = NextItem(uni);
    }

  /***************************************************************************/

  /***** Cell lists for surfaces *********************************************/

  /* NOTE: Disabled for now */

  loc0 = (long)RDB[DATA_PTR_S0];
  if (1 == 2)
    if ((loc0 = (long)RDB[DATA_PTR_S0]) > VALID_PTR)
      {
        /* Get list size */

        sz = ListSize(loc0);

        /* Loop over list */

#ifdef OPEN_MP
#pragma omp parallel for private(surf, ptr)
#endif

        for (i = 0; i < sz; i++)
          {
            /* Get pointer */

            surf = ListPtr(loc0, i);
            CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

            /* Sort list */

            if ((ptr = (long)RDB[surf + SURFACE_PTR_CELL_LIST]) > VALID_PTR)
              SortList(ptr, CELL_LIST_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
          }
      }

  /***************************************************************************/

  /***** Cell search meshes **************************************************/

  /* Loop over structures */

  loc0 = (long)RDB[DATA_PTR_CSM0];
  while (loc0 > VALID_PTR)
    {
      /* Get pointer to mesh */

      msh = (long)RDB[loc0 + CELL_MESH_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get size */

      nx = (long)RDB[msh + MESH_N0];
      ny = (long)RDB[msh + MESH_N1];
      nz = (long)RDB[msh + MESH_N2];

      /* Check size */

      CheckValue(FUNCTION_NAME, "size", "", nx*ny*nz, 1, 1000000000);

      /* Loop over mesh */

      for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
          for (k = 0; k < nz; k++)
            {
              /* Get pointer to list */

              lst = ReadMeshPtr(msh, i, j, k);
              CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

              lst = (long)RDB[lst];
              CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

              /* Sort list */

              SortList(lst, CELL_LIST_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
            }

      /* Pointer to next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** STL solid search meshes *********************************************/

  /* Loop over structures */

  loc0 = (long)RDB[DATA_PTR_STL0];
  while (loc0 > VALID_PTR)
    {
      /* Get pointer to mesh */

      msh = (long)RDB[loc0 + STL_PTR_SOLID_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get size */

      nx = (long)RDB[msh + MESH_N0];
      ny = (long)RDB[msh + MESH_N1];
      nz = (long)RDB[msh + MESH_N2];

      /* Check size */

      CheckValue(FUNCTION_NAME, "size", "", nx*ny*nz, 1, 1000000000);

      /* Loop over mesh */

      for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
          for (k = 0; k < nz; k++)
            {
              /* Get pointer to list */

              lst = ReadMeshPtr(msh, i, j, k);
              CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

              /* Sort list */

              if ((lst = (long)RDB[lst]) > VALID_PTR)
                SortList(lst, CELL_LIST_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
            }

      /* Pointer to next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Tet mesh search lists ***********************************************/

  /* Loop over interfaces */

  loc0 = (long)RDB[DATA_PTR_IFC0];
  while (loc0 > VALID_PTR)
    {
      /* Check type */

      if ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_TET_MESH)
        {
          /* Get pointer to search mesh */

          msh = (long)RDB[loc0 + IFC_PTR_SEARCH_MESH_LIST];
          CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

          /* Get size */

          nx = (long)RDB[msh + MESH_N0];
          ny = (long)RDB[msh + MESH_N1];
          nz = (long)RDB[msh + MESH_N2];

          /* Check size */

          CheckValue(FUNCTION_NAME, "size", "", nx*ny*nz, 1, 1000000000);

          /* Loop over mesh */

          for (i = 0; i < nx; i++)
            for (j = 0; j < ny; j++)
              for (k = 0; k < nz; k++)
                {
                  /* Get pointer to list */

                  lst = ReadMeshPtr(msh, i, j, k);
                  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

                  /* Check content and sort list */

                  if ((lst = (long)RDB[lst]) > VALID_PTR)
                    if (ListSize(lst) > 10000)
                      SortList(lst, SEARCH_MESH_PTR_CELL_COUNT,
                               SORT_MODE_DESCEND);
                }
        }

      /* Next interface */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Stochastic geometry search lists ************************************/

  /* Loop over geometries */

  loc0 = (long)RDB[DATA_PTR_PB0];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to search mesh */

      msh = (long)RDB[loc0 + PBED_PTR_SEARCH_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get size */

      nx = (long)RDB[msh + MESH_N0];
      ny = (long)RDB[msh + MESH_N1];
      nz = (long)RDB[msh + MESH_N2];

      /* Check size */

      CheckValue(FUNCTION_NAME, "size", "", nx*ny*nz, 1, 1000000000);

      /* Loop over mesh */

      for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
          for (k = 0; k < nz; k++)
            {
              /* Get pointer to list */

              lst = ReadMeshPtr(msh, i, j, k);
              CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

              /* Check content and sort list */

              if ((lst = (long)RDB[lst]) > VALID_PTR)
                SortList(lst, SEARCH_MESH_PTR_CELL_COUNT,
                         SORT_MODE_DESCEND_PRIVA);
            }

      /* Next geometry */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Misc. stuff *********************************************************/

  /* ICM break list */

  if ((icm = (long)RDB[DATA_PTR_ICM0]) > VALID_PTR)
    SortList(icm, ICM_BREAK_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);

  /***************************************************************************/
}

/*****************************************************************************/
