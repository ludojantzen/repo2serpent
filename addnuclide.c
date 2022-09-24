/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addnuclide.c                                   */
/*                                                                           */
/* Created:       2010/09/10 (JLe)                                           */
/* Last modified: 2019/03/28 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Adds a new nuclide based on information in ACE directory     */
/*              file and ENDF format decay file                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddNuclide:"

/*****************************************************************************/

long AddNuclide(char *zaid, long ZAI, char *lib, double T, long type,
                long TMS)
{
  long nuc, ace, ptr, rea, Z, A;
  double mem;

  /* Find ace data for nuclide */

  if ((ace = FindNuclideData(zaid, ZAI, lib, T, type, TMS)) > -1)
    {
      /* Tässä oletetaan että jos nuklidi on lisätty jostain ketjusta, sen */
      /* zaid == NULL, jolloin asetetaan flägi */

      if ((ace > VALID_PTR) && (zaid == NULL))
        SetOption(ace + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_NEW_DAUGHTERS);

      /* Return pointer */

      return ace;
    }

  /* Convert pointer */

  ace = -ace;

  /* Reset memory size */

  mem = RDB[DATA_TOTAL_BYTES];

  /* Avoid compiler warning */

  nuc = -1;

  /* Check type */

  if ((long)ACE[ace + ACE_TYPE] == NUCLIDE_TYPE_TRANSMUXS)
    Die(FUNCTION_NAME, "Transmutation type read as independent nuclide");
  else if ((long)ACE[ace + ACE_TYPE] == NUCLIDE_TYPE_PHOTON)
    {
      /***********************************************************************/

      /***** Create element for photon interaction data **********************/

      /* Create new block */

      nuc = NewItem(DATA_PTR_NUC0, NUCLIDE_BLOCK_SIZE);

      /* Copy values */

      WDB[nuc + NUCLIDE_PTR_NAME] = ACE[ace + ACE_PTR_NAME];
      WDB[nuc + NUCLIDE_TYPE] = ACE[ace + ACE_TYPE];

      /* Check photon transport mode */

      if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == NO)
        Error(0, "Photon transport data (%s) without photon transport mode",
              GetText(nuc + NUCLIDE_PTR_NAME));

      /* ZAI */

      WDB[nuc + NUCLIDE_ZAI] = ACE[ace + ACE_ZAI];

      /* ZA */

      WDB[nuc + NUCLIDE_ZA] = ACE[ace + ACE_ZA];

      /* Z */

      WDB[nuc + NUCLIDE_Z] = ACE[ace + ACE_ZA]/1000.0;

      /* Check A */

      if (ACE[ace + ACE_ZA] - 1000.0*RDB[nuc + NUCLIDE_Z] != 0.0)
        Die(FUNCTION_NAME, "Non-zero A for photon interaction data");

      /* Library ID */

      if (lib == NULL)
        WDB[nuc + NUCLIDE_PTR_LIB_ID] = ACE[ace + ACE_PTR_LIB_ID];
      else
        WDB[nuc + NUCLIDE_PTR_LIB_ID] = (double)PutText(lib);

      /* Pointer to ACE data */

      WDB[nuc + NUCLIDE_PTR_ACE] = (double)ace;
      WDB[nuc + NUCLIDE_PTR_DECAY_ACE] = NULLPTR;

      /* Read ACE data */

      ReadACEFile(nuc);

      /* Add flags */

      SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_PHOTON_DATA);

      /* Reset toxicities */

      WDB[nuc + NUCLIDE_SPEC_ING_TOX] = -1.0;
      WDB[nuc + NUCLIDE_SPEC_INH_TOX] = -1.0;

      /***********************************************************************/
    }
  else if ((long)ACE[ace + ACE_TYPE] != NUCLIDE_TYPE_DECAY)
    {
      /***********************************************************************/

      /***** Create transport or dosimetry nuclide ***************************/

      /* Create new block */

      nuc = NewItem(DATA_PTR_NUC0, NUCLIDE_BLOCK_SIZE);

      /* Copy name */

      WDB[nuc + NUCLIDE_PTR_NAME] = ACE[ace + ACE_PTR_NAME];

      /* Set type to DBRC and set flag or copy type from data */

      if (type == NUCLIDE_TYPE_DBRC)
        {
          WDB[nuc + NUCLIDE_TYPE] = (double)type;
          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DBRC);
        }
      else
        WDB[nuc + NUCLIDE_TYPE] = ACE[ace + ACE_TYPE];

      /* Reset toxicities */

      WDB[nuc + NUCLIDE_SPEC_ING_TOX] = -1.0;
      WDB[nuc + NUCLIDE_SPEC_INH_TOX] = -1.0;

      /* Allocate memory for previous collision velocity and relative energy */

      AllocValuePair(nuc + NUCLIDE_PREV_COL_Z2);
      AllocValuePair(nuc + NUCLIDE_PREV_COL_COS);
      AllocValuePair(nuc + NUCLIDE_PREV_COL_ER);
      AllocValuePair(nuc + NUCLIDE_PREV_COL_ET);
      AllocValuePair(nuc + NUCLIDE_PREV_COL_DT);
      AllocValuePair(nuc + NUCLIDE_PREV_COL_T);
      AllocValuePair(nuc + NUCLIDE_PREV_COL_TV_Z2);
      AllocValuePair(nuc + NUCLIDE_PREV_COL_TV_COS);
      AllocValuePair(nuc + NUCLIDE_PREV_COL_TV_ET);
      AllocValuePair(nuc + NUCLIDE_PREV_COL_TV_T);

      /* Check type */

      if ((long)ACE[ace + ACE_TYPE] == NUCLIDE_TYPE_SAB)
        {
          /* ZAI */

          WDB[nuc + NUCLIDE_ZAI] = 10.0*ACE[ace + ACE_BOUND_ZA];

          /* ZA */

          WDB[nuc + NUCLIDE_ZA] = ACE[ace + ACE_BOUND_ZA];

          /* Z, A and I */

          Z = (long)(ACE[ace + ACE_BOUND_ZA]/1000.0);
          A = (long)ACE[ace + ACE_BOUND_ZA] - 1000*Z;

          WDB[nuc + NUCLIDE_Z] = (double)Z;
          WDB[nuc + NUCLIDE_A] = (double)A;
          WDB[nuc + NUCLIDE_I] = 0;
        }
      else
        {
          /* ZAI */

          WDB[nuc + NUCLIDE_ZAI] = ACE[ace + ACE_ZAI];

          /* ZA */

          WDB[nuc + NUCLIDE_ZA] = ACE[ace + ACE_ZA];

          /* Z, A and I */

          Z = (long)(ACE[ace + ACE_ZA]/1000.0);
          A = (long)ACE[ace + ACE_ZA] - 1000*Z;

          WDB[nuc + NUCLIDE_Z] = (double)Z;
          WDB[nuc + NUCLIDE_A] = (double)A;
          WDB[nuc + NUCLIDE_I] = ACE[ace + ACE_I];
        }

      /* Temperature */

      if (T < 0.0)
        WDB[nuc + NUCLIDE_TEMP] = ACE[ace + ACE_TEMP];
      else
        WDB[nuc + NUCLIDE_TEMP] = T;

      /* Set TMS flag */

      if (TMS == NUCLIDE_FLAG_TMS)
        SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_TMS);

      /* Real xs temperature */

      WDB[nuc + NUCLIDE_XS_TEMP] = ACE[ace + ACE_TEMP];

      /* Reset TMS temperatures */

      WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = INFTY;
      WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = -INFTY;

      /* Library ID */

      if (lib == NULL)
        WDB[nuc + NUCLIDE_PTR_LIB_ID] = ACE[ace + ACE_PTR_LIB_ID];
      else
        WDB[nuc + NUCLIDE_PTR_LIB_ID] = (double)PutText(lib);

      /* Pointer to ACE data */

      WDB[nuc + NUCLIDE_PTR_ACE] = (double)ace;
      WDB[nuc + NUCLIDE_PTR_DECAY_ACE] = NULLPTR;

      /* Read ACE data */

      ReadACEFile(nuc);

      /* Set DBRC type to transport */

      if (type == NUCLIDE_TYPE_DBRC)
        WDB[nuc + NUCLIDE_TYPE] = ACE[ace + ACE_TYPE];

      /* Add flags */

      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
        SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DOSIMETRY_DATA);
      else
        SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_TRANSPORT_DATA);

      /* Find decay data */

      if (((ptr = FindNuclideData(NULL, (long)RDB[nuc + NUCLIDE_ZAI],
                                  GetText(nuc + NUCLIDE_PTR_LIB_ID),
                                  RDB[nuc + NUCLIDE_TEMP],
                                  NUCLIDE_TYPE_DECAY, TMS)) < 0) &&
          ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DOSIMETRY))
        {
          /* Convert pointer */

          ace = -ptr;

          /* Check type (lambda = -1 for structural nuclides) */

          if (ACE[ace + ACE_LAMBDA] > -1.0)
            {
              /* Pointer to decay data */

              WDB[nuc + NUCLIDE_PTR_DECAY_ACE] = (double)ace;

              /* Set data flag */

              SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DECAY_DATA);

              /* Process decay data */

              ProcessDecayData(nuc);
            }
        }
      else if (ptr > 0)
        {
          /* Existing decay data, get pointer */

          if ((ace = (long)RDB[ptr + NUCLIDE_PTR_DECAY_ACE]) > 0)
            {
              /* Check lambda */

              if (ACE[ace + ACE_LAMBDA] < 0.0)
                Die(FUNCTION_NAME, "lambda < 0");

              /* Pointer to decay data */

              WDB[nuc + NUCLIDE_PTR_DECAY_ACE] = (double)ace;

              /* Set data flag */

              SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DECAY_DATA);

              /* Process decay data */

              ProcessDecayData(nuc);
            }
        }

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Create decay nuclide ********************************************/

      /* Create new block */

      nuc = NewItem(DATA_PTR_NUC0, NUCLIDE_BLOCK_SIZE);

      /* Copy values */

      WDB[nuc + NUCLIDE_PTR_NAME] = ACE[ace + ACE_PTR_NAME];
      WDB[nuc + NUCLIDE_TYPE] = (double)NUCLIDE_TYPE_DECAY;

      /* ZAI */

      WDB[nuc + NUCLIDE_ZAI] = ACE[ace + ACE_ZAI];

      /* ZA */

      WDB[nuc + NUCLIDE_ZA] = ACE[ace + ACE_ZA];

      /* Z, A and I */

      Z = (long)(ACE[ace + ACE_ZA]/1000.0);
      A = (long)ACE[ace + ACE_ZA] - 1000*Z;

      WDB[nuc + NUCLIDE_Z] = (double)Z;
      WDB[nuc + NUCLIDE_A] = (double)A;
      WDB[nuc + NUCLIDE_I] = ACE[ace + ACE_I];

      /* Temperature */

      if (T < 0.0)
        WDB[nuc + NUCLIDE_TEMP] = ACE[ace + ACE_TEMP];
      else
        WDB[nuc + NUCLIDE_TEMP] = T;

      /* Set TMS flag */

      if (TMS == NUCLIDE_FLAG_TMS)
        SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_TMS);

      /* Real xs temperature */

      WDB[nuc + NUCLIDE_XS_TEMP] = ACE[ace + ACE_TEMP];

      /* Reset TMS temperatures */

      WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = INFTY;
      WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = -INFTY;

      /* Library ID */

      if (lib == NULL)
        WDB[nuc + NUCLIDE_PTR_LIB_ID] = ACE[ace + ACE_PTR_LIB_ID];
      else
        WDB[nuc + NUCLIDE_PTR_LIB_ID] = (double)PutText(lib);

      /* Pointer to reaction data */

      WDB[nuc + NUCLIDE_PTR_ACE] = NULLPTR;

      /* Reset toxicities */

      WDB[nuc + NUCLIDE_SPEC_ING_TOX] = -1.0;
      WDB[nuc + NUCLIDE_SPEC_INH_TOX] = -1.0;

      /* Check type */

      if (ACE[ace + ACE_LAMBDA] < 0)
        {
          /* Structural type, put masses only */

          WDB[nuc + NUCLIDE_AWR] = ACE[ace + ACE_AWR];
          WDB[nuc + NUCLIDE_AW] = M_NEUTRON*ACE[ace + ACE_AWR];
        }
      else
        {
          /* Pointer to decay data */

          WDB[nuc + NUCLIDE_PTR_DECAY_ACE] = (double)ace;

          /* Process decay data */

          ProcessDecayData(nuc);

          /* Set flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DECAY_DATA);
        }

      /* Read additional transmutation cross section data */

      if ((ace = FindNuclideData(NULL, (long)RDB[nuc + NUCLIDE_ZAI], NULL,
                                 -1.0, NUCLIDE_TYPE_TRANSMUXS, TMS)) < 0)
        {
          /* Convert pointer */

          ace = -ace;

          /* Pointer to ACE data */

          WDB[nuc + NUCLIDE_PTR_ACE] = (double)ace;

          /* Add flag  */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_TRANSMU_DATA);

          /* Read ACE data */

          ReadACEFile(nuc);
        }

      /***********************************************************************/
    }

  /* Check nuclide pointer */

  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Add fission yields */

  ProcessFissionYields(nuc);

  /* Add branching */

  AddBranching(nuc);

  /* Set mt 4 type to special (set to partial in readacefile.c to */
  /* include isomaric state production from inelastic scattering. */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check mt */

      if (((long)RDB[rea + REACTION_MT] == 4) &&
          ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL))
        WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;

      /* Next */

      rea = NextItem(rea);
    }

  /* Update counters */

  ReactionCount();

  /* Set memory size */

  WDB[nuc + NUCLIDE_MEMSIZE] = RDB[DATA_TOTAL_BYTES] - mem;

  /* Print data */

  PrintNuclideData(nuc, NO);

  /* Put large number to level */

  if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON) ||
      (type == NUCLIDE_TYPE_DBRC))
    WDB[nuc + NUCLIDE_PATH_LEVEL] = 0.0;
  else
    WDB[nuc + NUCLIDE_PATH_LEVEL] = 1E+6;

  /* Tässä oletetaan että jos nuklidi on lisätty jostain ketjusta, sen */
  /* zaid == NULL, jolloin asetetaan flägi */

  if (zaid == NULL)
    SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_NEW_DAUGHTERS);

  /* Return pointer */

  return nuc;
}

/*****************************************************************************/
