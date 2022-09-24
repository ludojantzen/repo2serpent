/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : detbin.c                                       */
/*                                                                           */
/* Created:       2011/07/11 (JLe)                                           */
/* Last modified: 2019/10/22 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Finds detector bin                                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DetBin:"

/*****************************************************************************/

long DetBin(long det, long mat, long part, double x, double y, double z,
            double E, double t, long id)
{
  long ptr, uni, lat, cell, umsh, ifc, idx, ncol, lvl0, lvl, flg, tst, n, m;
  long ebin, ubin, cbin, mbin, lbin, ibin, tbin, msh, prnt, tra;
  long ne, nu, nc, nm, nl, ni, nt, nmax;
  double x1, y1, z1, u, v, w;

  /* Reset bins */

  ebin = 0;
  ubin = 0;
  cbin = 0;
  mbin = 0;
  lbin = 0;
  ibin = 0;
  tbin = 0;

  /* Get number of bins */

  ne = (long)RDB[det + DET_N_EBINS];
  nu = (long)RDB[det + DET_N_UBINS];
  nc = (long)RDB[det + DET_N_CBINS];
  nm = (long)RDB[det + DET_N_MBINS];
  nl = (long)RDB[det + DET_N_LBINS];
  nt = (long)RDB[det + DET_N_TBINS];

  /* Mesh bins */

  if ((ptr = (long)RDB[det + DET_PTR_MESH]) > VALID_PTR)
    ni = (long)(RDB[ptr + MESH_N0]*RDB[ptr + MESH_N1]*RDB[ptr + MESH_N2]);
  else
    ni = 1;

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Check fissile option */

  if ((long)RDB[det + DET_SCORE_FISS_REG_ONLY] == YES)
    {
      /* Continue only if material is fissile */

      if (mat < VALID_PTR)
        return -1;
      else if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT))
        return -1;
    }

  /***************************************************************************/

  /***** Test flagging *******************************************************/

  /* Check particle pointer */

  if (part > VALID_PTR)
    {
      /* Get flags */

      flg = (long)RDB[part + PARTICLE_DET_FLAGS];

      /* Reset counters */

      n = 0;
      m = 0;

      /* Loop over flags */

      ptr = (long)RDB[det + DET_PTR_FLAGGING];
      while (ptr > VALID_PTR)
        {
          /* Get flag number */

          tst = (long)RDB[ptr + DET_FBIN_FLAG_NUMBER];
          CheckValue(FUNCTION_NAME, "tst", "", tst, 1, 64);

          /* Convert */

          tst = (long)(pow(2.0, (double)tst - 1.0));
          CheckValue(FUNCTION_NAME, "tst", "", tst, 1, LONG_MAX);

          /* Check test option */

          if ((long)RDB[ptr + DET_FBIN_FLAG_OPTION] ==
              DET_FLAG_OPT_TEST_UNSET)
            {
              /* Add to test count */

              m++;

              /* Check if flag is not set */

              if (!(flg & tst))
                n++;
              else if ((long)RDB[ptr + DET_FBIN_FLAG_AND_LOGIC] == YES)
                {
                  /* Reset count */

                  n = -1;

                  /* Break loop */

                  break;
                }
            }
          else if ((long)RDB[ptr + DET_FBIN_FLAG_OPTION] ==
                   DET_FLAG_OPT_TEST_SET)
            {

              /* Add to test count */

              m++;

              /* Check if flag is set */

              if (flg & tst)
                n++;
              else if ((long)RDB[ptr + DET_FBIN_FLAG_AND_LOGIC] == YES)
                {
                  /* Reset count */

                  n = -1;

                  /* Break loop */

                  break;
                }
            }

          /* Next flag */

          ptr = NextItem(ptr);
        }

      /* Check if flags are tested but count is zero */

      if ((m > 0) && (n <= 0))
        return -1;
    }

  /***************************************************************************/

  /***** Energy bin **********************************************************/

  /* Pointer to detector energy grid */

  if ((ptr = (long)RDB[det + DET_PTR_EGRID]) > VALID_PTR)
    {
      /* PulseDet() passes E = -INFTY to check that collision occurs    */
      /* in detector. The deposited energy is not binned at this point. */
      /* (JLe 30.6.2015 / 2.1.25) */

      if (E == -INFTY)
        {
          /* Set bin index */

          ebin = 0;
        }
      else if ((ebin = GridSearch(ptr, E)) < 0)
        {
          /* Out of bounds */

          return -1;
        }

      /* Check value */

      CheckValue(FUNCTION_NAME, "ebin", "", ebin, 0, ne - 1);
    }

  /***************************************************************************/

  /***** Time bin ************************************************************/

  /* Pointer to detector time binning */

  if ((ptr = (long)RDB[det + DET_PTR_TME]) > VALID_PTR)
    {
      /* Get bin */

      if ((tbin = SearchArray(&RDB[ptr], t, nt + 1)) < 0)
        {
          /* Out of bounds */

          return -1;
        }

      /* Check value */

      CheckValue(FUNCTION_NAME, "tbin", "", tbin, 0, nt - 1);
    }

  /***************************************************************************/

  /***** Universe bins *******************************************************/

  /* Check if universe bins are defined */

  if ((ptr = (long)RDB[det + DET_PTR_UBINS]) > VALID_PTR)
    {
      /* Loop over universe bins */

      while (ptr > VALID_PTR)
        {
          /* Universe pointer */

          uni = (long)RDB[ptr + DET_UBIN_PTR_UNI];
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Check universe type */

          if ((long)RDB[uni + UNIVERSE_TYPE] == UNIVERSE_TYPE_SUPER)
            {
              /* Super-imposed universe, test prevuous and point */

              if (TestValuePair(uni + UNIVERSE_COL_COUNT, (double)ncol, id)
                  == 1.0)
                break;
              else if (InSuperCell(uni, -1, x, y, z, id) == YES)
                break;
            }
          else
            {
              /* Physical universe, check for collision */

              if (TestValuePair(uni + UNIVERSE_COL_COUNT, (double)ncol, id)
                  > 0.0)
                break;
            }

          /* Update bin */

          ubin++;

          /* Next bin */

          ptr = NextItem(ptr);
        }

      /* Check pointer */

      if (ptr < VALID_PTR)
        {
          /* Not inside */

          return -1;
        }

      /* Check value */

      CheckValue(FUNCTION_NAME, "ubin", "", ubin, 0, nu - 1);
    }

  /***************************************************************************/

  /***** Lattice *************************************************************/

  /* Check if lattice bins are defined */

  if ((ptr = (long)RDB[det + DET_PTR_LBINS]) > VALID_PTR)
    {
      /* Pointer to lattice */

      lat = (long)RDB[ptr + DET_LBIN_PTR_LAT];
      CheckPointer(FUNCTION_NAME, "(lat)", DATA_ARRAY, lat);

      /* Check type */

      if ((long)RDB[ptr + DET_LBIN_TYPE] == DET_LBIN_LAT)
        {
          /* Lattice, get index */

          if ((lbin = (long)TestValuePair(lat + LAT_COL_CELL_IDX, (double)ncol,
                                          id)) < 0)
            return -1;
        }
      else if ((long)RDB[ptr + DET_LBIN_TYPE] == DET_LBIN_PBED)
        {
          /* Get pointer to unit where collision occurred */

          if ((ptr = (long)TestValuePair(lat + PBED_PTR_COL_PEBBLE,
                                         (double)ncol, id)) > VALID_PTR)
          {
            /* Get pebble index */

            lbin = (long)RDB[ptr + PEBBLE_IDX];
          }
          else
            return -1;
        }
      else
        Die(FUNCTION_NAME, "Invalid lattice type");

      /* Check value */

      CheckValue(FUNCTION_NAME, "lbin", "", lbin, 0, nl - 1);
    }

  /***************************************************************************/

  /***** Cell bins ***********************************************************/

  /* Check if cell bins are defined */

  if ((ptr = (long)RDB[det + DET_PTR_CBINS]) > VALID_PTR)
    {
      /* Check pointer to UMSH geometry */

      if ((umsh = (long)RDB[ptr + DET_CBIN_UMSH_PTR_UMSH]) > VALID_PTR)
        {
          /* Get pointer to interface */

          ifc = (long)RDB[umsh + UMSH_PTR_IFC];
          CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

          /* Get pointer to tet cell (set in whereami.c) */

          if ((ptr = (long)TestValuePair(ifc + IFC_PTR_PREV_COL_CELL,
                                         (double)ncol, id)) > VALID_PTR)
            {
              /* Get pointer to parent */

              prnt = (long)RDB[ptr + TET_PTR_PARENT];
              CheckPointer(FUNCTION_NAME, "(prnt)", DATA_ARRAY, prnt);

              /* Get index to statistics and score */

              if ((cbin = (long)RDB[prnt + IFC_TET_PRNT_STAT_IDX]) < 0)
                {
                  /* Not mapped */

                  return -1;
                }
            }
          else
            {
              /* Not in any cell */

              return -1;
            }
        }
      else if ((long)RDB[ptr + DET_CBIN_SUPER_CELL] == YES)
        {
          /* Loop over cell bins */

          while (ptr > VALID_PTR)
            {
              /* Cell pointer */

              cell = (long)RDB[ptr + DET_CBIN_PTR_CELL];
              CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

              /* Check cell type */

              if ((long)RDB[ptr + DET_CBIN_SUPER_CELL] == YES)
                {
                  /* Super-imposed cell, test point */

                  if (InCell(cell, x, y, z, NO, id) == YES)
                    break;
                }
              else
                Die(FUNCTION_NAME, "Physical cell in list");

              /* Update bin */

              cbin++;

              /* Next bin */

              ptr = NextItem(ptr);
            }

          /* Check pointer */

          if (ptr < VALID_PTR)
            {
              /* Not inside */

              return -1;
            }
        }
      else
        {
          /* Reset cell bin */

          cbin = -1;

          /* Loop over levels */

          lvl0 = (long)RDB[DATA_PTR_LVL0];
          while (lvl0 > VALID_PTR)
            {
              /* Pointer to private data */

              lvl = (long)RDB[lvl0 + LVL_PTR_PRIVATE_DATA];
              CheckPointer(FUNCTION_NAME, "(lvl)", PRIVA_ARRAY, lvl);

              /* Pointer to cell and detector bin */

              if ((cell = (long)GetPrivateData(lvl + LVL_PRIV_PTR_CELL, id))
                  > VALID_PTR)
                {
                  /* Loop over detector bin list */

                  ptr = (long)RDB[cell + CELL_PTR_DETBIN];
                  while (ptr > VALID_PTR)
                    {
                      /* Check detector pointer */

                      if (det == (long)RDB[ptr + DETBIN_PTR_DET])
                        {
                          /* Set bin */

                          cbin = (long)RDB[ptr + DETBIN_BIN];

                          /* Break loop */

                          break;
                        }

                      /* Next bin */

                      ptr = NextItem(ptr);
                    }
                }

              /* Check if last */

              if (GetPrivateData(lvl + LVL_PRIV_LAST, id) == YES)
                break;

              /* Next level */

              lvl0 = NextItem(lvl0);
            }

          /* Check cell bin */

          if (cbin < 0)
            return -1;
        }

      /* Check value */

      CheckValue(FUNCTION_NAME, "cbin", "", cbin, 0, nc - 1);
    }

  /***************************************************************************/

  /***** Material bins *******************************************************/

  /* Check if material bins are defined */

  if ((ptr = (long)RDB[det + DET_PTR_MBINS]) > VALID_PTR)
    {
      /* Check material pointer */

      if (mat < VALID_PTR)
        return -1;

      /* Reset bin */

      mbin = -1;

      /* Loop over detector bin list */

      ptr = (long)RDB[mat + MATERIAL_PTR_DETBIN];
      while (ptr > VALID_PTR)
        {
          /* Check detector pointer */

          if (det == (long)RDB[ptr + DETBIN_PTR_DET])
            {
              /* Set bin */

              mbin = (long)RDB[ptr + DETBIN_BIN];

              /* Break loop */

              break;
            }

          /* Next bin list */

          ptr = NextItem(ptr);
        }

      /* Repeat for parent */

      if ((mbin < 0) &&
          ((mat = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR))
        {
          /* Loop over detector bin list */

          ptr = (long)RDB[mat + MATERIAL_PTR_DETBIN];
          while (ptr > VALID_PTR)
            {
              /* Check detector pointer */

              if (det == (long)RDB[ptr + DETBIN_PTR_DET])
                {
                  /* Set bin */

                  mbin = (long)RDB[ptr + DETBIN_BIN];

                  /* Break loop */

                  break;
                }

              /* Next bin list */

              ptr = NextItem(ptr);
            }
        }

      /* Check bin */

      if (mbin < 0)
        return -1;

#ifdef WANHAA

      /* Loop over material bins */

      while (ptr > VALID_PTR)
        {
          /* Compare material pointer */

          if ((long)RDB[ptr + DET_MBIN_PTR_MAT] == mat)
            break;

          /* Compare to parent (if divided in burnup calculation) */

          if ((long)RDB[ptr + DET_MBIN_PTR_MAT] ==
              (long)RDB[mat + MATERIAL_DIV_PTR_PARENT])
            break;

          /* Update bin */

          mbin++;

          /* Next bin */

          ptr = NextItem(ptr);
        }

      /* Check pointer */

      if (ptr < VALID_PTR)
        {
          /* Not inside */

          return -1;
        }

#endif

      /* Check value */

      CheckValue(FUNCTION_NAME, "mbin", "", mbin, 0, nm - 1);
    }

  /***************************************************************************/

  /***** Mesh bins ***********************************************************/

  /* Get mesh bin */

  if ((msh = (long)RDB[det + DET_PTR_MESH]) > VALID_PTR)
    {
      /* Check if we should pass local coordinates */

      if (RDB[msh + MESH_LOCAL_COORDS] == (double)YES)
        {
          /* Get local coordinates from collision universe */

          /* Get collision universe */

          ptr = (long)RDB[DATA_PTR_COLLISION_UNI];
          uni = (long)GetPrivateData(ptr, id);

          /* Check pointer */

          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Get coordinates */

          ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_X];
          CheckPointer(FUNCTION_NAME, "(xptr)", PRIVA_ARRAY, ptr);
          x1 = GetPrivateData(ptr, id);

          ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Y];
          CheckPointer(FUNCTION_NAME, "(yptr)", PRIVA_ARRAY, ptr);
          y1 = GetPrivateData(ptr, id);

          ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Z];
          CheckPointer(FUNCTION_NAME, "(zptr)", PRIVA_ARRAY, ptr);
          z1 = GetPrivateData(ptr, id);

          /* Check coordinate transformation */

          if ((tra = (long)RDB[det + DET_PTR_TRANS]) > VALID_PTR)
            {
              /* Set dummy direction cosines */

              u = 1.0;
              v = 0.0;
              w = 0.0;

              /* Call transformation */

              CoordTrans(tra, &x1, &y1, &z1, &u, &v, &w, id);
            }

          /* Get bin */

          if ((ibin = MeshIndex(msh, x1, y1, z1, -1.0)) < 0)
            return -1;
        }
      else
        {
          /* Check coordinate transformation */

          if ((tra = (long)RDB[det + DET_PTR_TRANS]) > VALID_PTR)
            {
              /* Set dummy direction cosines */

              u = 1.0;
              v = 0.0;
              w = 0.0;

              /* Call transformation */

              CoordTrans(tra, &x, &y, &z, &u, &v, &w, id);
            }

          /* Get bin */

          if ((ibin = MeshIndex(msh, x, y, z, -1.0)) < 0)
            return -1;
        }
    }

  /***************************************************************************/

  /***** Calculate bin index *************************************************/

  /* Check bins */

#ifdef DEBUG

  if ((ebin < 0) || (ebin > (long)RDB[det + DET_N_EBINS]))
    Die(FUNCTION_NAME, "Invalid ebin");
  if ((ubin < 0) || (ubin > (long)RDB[det + DET_N_UBINS]))
    Die(FUNCTION_NAME, "Invalid ubin");
  if ((cbin < 0) || (cbin > (long)RDB[det + DET_N_CBINS]))
    Die(FUNCTION_NAME, "Invalid cbin");
  if ((mbin < 0) || (mbin > (long)RDB[det + DET_N_MBINS]))
    Die(FUNCTION_NAME, "Invalid mbin");
  if ((lbin < 0) || (lbin > (long)RDB[det + DET_N_LBINS]))
    Die(FUNCTION_NAME, "Invalid lbin");
  if ((tbin < 0) || (tbin > (long)RDB[det + DET_N_TBINS]))
    Die(FUNCTION_NAME, "Invalid tbin");

#endif

  /* Check if bins are combined */

  if ((long)RDB[det + DET_TYPE] == DETECTOR_TYPE_ADD_BINS)
    {
      /* Reset cell and material bins */

      cbin = 0;
      mbin = 0;
    }

  /* Calculate index */

  nmax = 1;
  idx = 0;

  idx = idx + ebin*nmax;
  nmax = nmax*ne;

  idx = idx + ubin*nmax;
  nmax = nmax*nu;

  idx = idx + cbin*nmax;
  nmax = nmax*nc;

  idx = idx + mbin*nmax;
  nmax = nmax*nm;

  idx = idx + lbin*nmax;
  nmax = nmax*nl;

  idx = idx + ibin*nmax;
  nmax = nmax*ni;

  idx = idx + tbin*nmax;
  nmax = nmax*nt;

  /* Return index */

  return idx;
}

/*****************************************************************************/
