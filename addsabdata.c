/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addsabdata.c                                   */
/*                                                                           */
/* Created:       2011/01/22 (JLe)                                           */
/* Last modified: 2018/05/25 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: - Adds S(a,b) channels to transport nuclide                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddSabData:"

/*****************************************************************************/

void AddSabData()
{
  long mat, loc0, loc1, loc2, iso, nuc, ptr, rea, new, yld, tgt, sab, rea0;
  long rea1;
  char tmpstr[MAX_STR];
  double f;

  /***************************************************************************/

  /***** Make duplicates and add channels ************************************/

  /* Check if S(a,b) data exists */

  if ((long)RDB[DATA_PTR_T0] < VALID_PTR)
    return;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check divisor */

      if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
        {
          /* Next material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Loop over S(a,b) data */

      loc0 = (long)RDB[mat + MATERIAL_PTR_SAB];
      while (loc0 > VALID_PTR)
        {
          /* Pointer to composition */

          iso = (long)RDB[loc0 + THERM_PTR_COMP];
          CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check if S(a,b) nuclide already created and linked to nuc. If it is, skip. */
          /* This occurs if S(a,b) data is linked to the same nuclide (for example 6000.09c)
             in multiple materials */

          /* HUOM! Jos sama S(a,b) -data on linkitetty useammalle nuklidille (esim. 6000.09c ja
             6000.12c), niin pitää luoda uusi S(a,b)-nuklidi jolla on orkkisnuklidin muut reaktiot,
             mutta esim. joitakin AllocValuePaireja ei saa tehdä uudestaan S(a,b)-datalle */

          if ((long)RDB[nuc + NUCLIDE_PTR_SAB_NUC] > VALID_PTR)
            {

              loc1 = (long)RDB[loc0 + THERM_PTR_THERM];
              sab = (long)RDB[loc1 + THERM_PTR_SAB];

              if ((long)RDB[sab + SAB_PTR_PREV_FRAC] > VALID_PTR)
                {
                  /* Replace original nuclide in composition with S(a,b) nuclide */

                  WDB[iso + COMPOSITION_PTR_NUCLIDE] =
                    RDB[nuc + NUCLIDE_PTR_SAB_NUC];

                  loc0 = NextItem(loc0);

                  continue;
                }
            }

          /* Create new nuclide for combined data */

          new = NewItem(DATA_PTR_NUC0, NUCLIDE_BLOCK_SIZE);

          /* Copy data */

          memcpy(&WDB[new + LIST_DATA_SIZE], &RDB[nuc + LIST_DATA_SIZE],
                 (NUCLIDE_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

          /* Put pointer */

          WDB[nuc + NUCLIDE_PTR_SAB_NUC] = (double)new;
          WDB[new + NUCLIDE_SAB_PTR_FREE] = (double)nuc;

          /* Reset toxicities */

          WDB[new + NUCLIDE_SPEC_ING_TOX] = -1.0;
          WDB[new + NUCLIDE_SPEC_INH_TOX] = -1.0;

          /* Allocate memory for previous collision velocity and relative */
          /* energy */

          WDB[new + NUCLIDE_PREV_COL_Z2] = NULLPTR;
          WDB[new + NUCLIDE_PREV_COL_COS] = NULLPTR;
          WDB[new + NUCLIDE_PREV_COL_ER] = NULLPTR;
          WDB[new + NUCLIDE_PREV_COL_ET] = NULLPTR;
          WDB[new + NUCLIDE_PREV_COL_DT] = NULLPTR;
          WDB[new + NUCLIDE_PREV_COL_T] = NULLPTR;
          WDB[new + NUCLIDE_PREV_COL_TV_Z2] = NULLPTR;
          WDB[new + NUCLIDE_PREV_COL_TV_COS] = NULLPTR;
          WDB[new + NUCLIDE_PREV_COL_TV_ET] = NULLPTR;
          WDB[new + NUCLIDE_PREV_COL_TV_T] = NULLPTR;

          AllocValuePair(new + NUCLIDE_PREV_COL_Z2);
          AllocValuePair(new + NUCLIDE_PREV_COL_COS);
          AllocValuePair(new + NUCLIDE_PREV_COL_ER);
          AllocValuePair(new + NUCLIDE_PREV_COL_ET);
          AllocValuePair(new + NUCLIDE_PREV_COL_DT);
          AllocValuePair(new + NUCLIDE_PREV_COL_T);
          AllocValuePair(new + NUCLIDE_PREV_COL_TV_Z2);
          AllocValuePair(new + NUCLIDE_PREV_COL_TV_COS);
          AllocValuePair(new + NUCLIDE_PREV_COL_TV_ET);
          AllocValuePair(new + NUCLIDE_PREV_COL_TV_T);

          /* Reset reaction pointer */

          WDB[new + NUCLIDE_PTR_REA] = NULLPTR;

          /* Taa on varmaan turha ja pitaisi olla new -nuklidille muutenkin ?
             TVi 2015-08-06 */

          /* Reset used-flag (NOTE: idea on että poistetaan nuklidi */
          /* muistin säästämiseksi jos sitä ei ole syntynyt minkään */
          /* muun nuklidin ketjussa. Vaatii sen että tätä rutiinia  */
          /* kutsutaan viimeiseksi, ja että noi flägit on asetettu. */
          /* Eikä tää välttämättä silti toimi, vaan saattaa johtaa  */
          /* siihen että nuklidia ei löydy MakeBurnMatrix():ssa.)   */

          if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_AP) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DP) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FP) &&
              !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_BP))
            ResetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

          /* Change name */

          sprintf(tmpstr, "%s", GetText(nuc + NUCLIDE_PTR_NAME));
          tmpstr[strlen(tmpstr) - 1] = 's';
          WDB[new + NUCLIDE_PTR_NAME] = (double)PutText(tmpstr);

          /* Change library ID */

          sprintf(tmpstr, "%s", GetText(nuc + NUCLIDE_PTR_LIB_ID));
          tmpstr[strlen(tmpstr) - 1] = 's';
          WDB[new + NUCLIDE_PTR_LIB_ID] = (double)PutText(tmpstr);

          /* Set S(a,b) flag */

          SetOption(new + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_SAB_DATA);

          /* Copy reaction data (this is needed to get the root */
          /* pointers right) */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* New reaction */

              ptr = NewItem(new + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);

              /* Copy data */

              memcpy(&WDB[ptr + LIST_DATA_SIZE], &RDB[rea + LIST_DATA_SIZE],
                     (REACTION_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

              WDB[ptr + REACTION_PTR_NUCLIDE] = (double)new;

              /* Put direct pointers */

              if ((long)RDB[ptr + REACTION_MT] == 301)
                WDB[new + NUCLIDE_PTR_HEATPRODXS] = (double)ptr;
              else if ((long)RDB[ptr + REACTION_MT] == 202)
                WDB[new + NUCLIDE_PTR_PHOTPRODXS] = (double)ptr;

              /* Next reaction */

              rea = NextItem(rea);
            }

          /* Put pointer */

          WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)new;

          /* Pointer to nuclide */

          loc1 = (long)RDB[loc0 + THERM_PTR_THERM];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* OTF-interpoloinnin tapauksessa sortataan SAB-lista ja luodaan */
          /* emänuklidille (new) yksi tai kaksi OTF-S(a,b) -reaktiota */

          if ((long)RDB[loc1 + THERM_INTERP_MODE] == THERM_INTERP_OTF)
            {
              /* Get pointer to Sab list (needed in the on-the-fly treatment) */

              sab = (long)RDB[loc1 + THERM_PTR_SAB];

              /***************************************************************/

              /* Sort SAB list according to temperature */

              SortList(sab, SAB_T, SORT_MODE_ASCEND);

              /***************************************************************/

              /* Check that the given libraries cover the whole temperature */
              /* range */

              sab = LastItem(sab);

              if (RDB[loc1 + THERM_OTF_MAX_TEMP] > RDB[sab + SAB_T])
                Error(loc1, "Maximum interpolation T above highest T found in S(a,b) data (%.2f > %.2f)",
                      RDB[loc1 + THERM_OTF_MAX_TEMP], RDB[sab + SAB_T]);

              sab = FirstItem(sab);

              if (RDB[loc1 + THERM_OTF_MIN_TEMP] < RDB[sab + SAB_T])
                Error(loc1, "Minimum interpolation T below smallest T found in S(a,b) data (%.2f < %.2f)",
                      RDB[loc1 + THERM_OTF_MIN_TEMP], RDB[sab + SAB_T]);

              /* Set pointers to first item in list */

              WDB[new + NUCLIDE_PTR_SAB] = (double)sab;
              WDB[loc1 + THERM_PTR_SAB] = (double)sab;

              /* Set SAB_EMAX based on first S(a,b) nuclide */
              /* Pitäisi olla sama kaikille lämpötiloille, voisi ehkä */
              /* tarkistaa joskus aikaisemmin */

              iso = (long)RDB[sab + SAB_PTR_ISO];
              CheckPointer(FUNCTION_NAME, "(iso0)" , DATA_ARRAY, iso);

              WDB[new + NUCLIDE_SAB_EMAX] = RDB[iso + NUCLIDE_EMAX];

              /* New: add placeholder S(a,b) reactions */

              rea0 = (long)RDB[new + NUCLIDE_PTR_REA];
              rea1 = (long)RDB[iso + NUCLIDE_PTR_REA];

              CheckPointer(FUNCTION_NAME, "(rea0)" , DATA_ARRAY, rea0);

              while(rea1 > VALID_PTR)
                {
                  /* Add new reaction */

                  rea = NewItem(new + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);

                  /* Copy data */

                  memcpy(&WDB[rea + LIST_DATA_SIZE],
                         &RDB[rea0 + LIST_DATA_SIZE],
                         (REACTION_BLOCK_SIZE - LIST_DATA_SIZE)
                         *sizeof(double));

                  /* Set MT as mt + 1000 (2004 or 2002) */

                  WDB[rea + REACTION_MT] = RDB[rea1 + REACTION_MT]+1000;
                  WDB[rea + REACTION_PTR_NUCLIDE] = (double)new;

                  /* Compare to nuclide minimum */

                  if (RDB[rea + REACTION_EMIN] < RDB[new + NUCLIDE_EMIN])
                    WDB[new + NUCLIDE_EMIN] = RDB[rea + REACTION_EMIN];

                  /* Compare to nuclide maximum */

                  if (RDB[rea + REACTION_EMAX] > RDB[new + NUCLIDE_EMAX])
                    WDB[new + NUCLIDE_EMAX] = RDB[rea + REACTION_EMAX];

                  rea1 = NextItem(rea1);
                }

              /* Pointer to sab data */

              sab = (long)RDB[new + NUCLIDE_PTR_SAB];

              /* Allocate value pairs for first item */

              /* Jos turha S(a,b) data päädytään jossain vaiheessa     */
              /* poistamaan niin tää pitää sitten vaihtaa ensimmäiseen */
              /* sab-itemiin */

              /* Check that the value pairs have not yet been allocated
                 (in case the same S(a,b) data is linked to multiple nuclides,
                 for example 6000.09c and 6000.12c) */

              if((long)RDB[sab + SAB_PTR_PREV_FRAC] < VALID_PTR){
                AllocValuePair(sab + SAB_PTR_PREV_FRAC);
                AllocValuePair(sab + SAB_PTR_PREV_SAB1);
              }

              /* Set used-flags for S(a,b) nuclides within relevant */
              /* temperature region */
              /* TODO: lämpötila-aluetta ei tarkisteta vielä */

              while(sab > VALID_PTR)
                {
                  iso = (long)RDB[sab + SAB_PTR_ISO];
                  CheckPointer(FUNCTION_NAME, "(iso1)" , DATA_ARRAY, iso);

                  SetOption(iso + NUCLIDE_OPTIONS, OPT_USED);

                  sab = NextItem(sab);
                }
            }

          /* Muissa kuin OTF-interpolointimoodeissa (tai ei interpolointia  */
          /* -moodissa) lisätään uudelle S(a,b) -nuklidille vaan reaktiot   */
          /* S(a,b)-datasta ja stokastisen miksauksen tapauksessa säädetään */
          /* näiden reaktioiden vaikutusalat myöhemmin lämpötilojen         */
          /* suhteessa */

          else
            {

              /* Check that TMS is not used together with pre-interpolation or
                 no interpolation modes */

              if((long)RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE)
                Error(mat, "Only on-the-fly thermal scattering interpolation can be used with TMS");

              /* Loop over S(a,b) nuclides */
              /* Listassa on 2 nuklidia jos stokastinen interpolointi on */
              /* käytössä. Muussa tapauksessa nuklideita pitäisi olla vain 1 */

              loc2 = (long)RDB[loc1 + THERM_PTR_SAB];
              while(loc2 > VALID_PTR)
                {
                  f = RDB[loc2 + SAB_FRAC];
                  loc1 = (long)RDB[loc2 + SAB_PTR_ISO];

                  /* Check pointer */

                  if (loc1 > VALID_PTR)
                    {
                      /* Loop over reactions */

                      rea = (long)RDB[loc1 + NUCLIDE_PTR_REA];
                      while(rea > VALID_PTR)
                        {
                          /* Add new reaction */

                          ptr = NewItem(new + NUCLIDE_PTR_REA,
                                        REACTION_BLOCK_SIZE);

                          /* Copy data */

                          memcpy(&WDB[ptr + LIST_DATA_SIZE],
                                 &RDB[rea + LIST_DATA_SIZE],
                                 (REACTION_BLOCK_SIZE - LIST_DATA_SIZE)
                                 *sizeof(double));

                          /* Set temperature adjustment factor */

                          WDB[ptr + REACTION_SAB_FRAC] = f;

                          /* Compare to nuclide minimum */

                          if (RDB[rea + REACTION_EMIN] <
                              RDB[nuc + NUCLIDE_EMIN])
                            WDB[nuc + NUCLIDE_EMIN] = RDB[rea + REACTION_EMIN];

                          /* Compare to nuclide maximum */

                          if (RDB[rea + REACTION_EMAX] >
                              RDB[nuc + NUCLIDE_EMAX])
                            WDB[nuc + NUCLIDE_EMAX] = RDB[rea + REACTION_EMAX];

                          /* Next reaction */

                          rea = NextItem(rea);
                        }
                    }

                  /* Next S(a,b) */

                  loc2=NextItem(loc2);
                }
            }

          /* Next therm */

          loc0 = NextItem(loc0);
        }


      /* Next material */

      mat = NextItem(mat);
    }

  /* Remove S(a,b) nuklidit (kokeillaan tätä nyt) */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check type and remove */

      /* Ei poisteta kuitenkaan jos OPT_USED-flag merkitty */

      if ( (long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_SAB &&
           !((long)RDB[nuc + NUCLIDE_OPTIONS] & OPT_USED) )
        {

          /* Copy pointer */

          ptr = nuc;

          /* Pointer to next */

          nuc = NextItem(nuc);

          /* Remove nuclide */

          RemoveItem(ptr);
        }
      else
        nuc = NextItem(nuc);
    }

  /* Set used-flags again for initial composition (Tää tarvitaan että */
  /* poistettu nuklidi luetaan mukaan jos sitä on toisen materiaalin  */
  /* alkuperäiskoostumuksessa.) */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Loop over composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];

          /* Set used-flag */

          SetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

          /* Next isotope */

          iso = NextItem(iso);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Reconfigure transmutation paths *************************************/

  /* Loop over materials in order to loop over initial nuclides */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check divisor */

      if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
        {
          /* Next material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check, whether material contains S(a,b) data */

      if((long)RDB[mat + MATERIAL_PTR_SAB] < VALID_PTR)
        {
          /* Next material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Pointer to composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Loop over material composition */

      while (iso > VALID_PTR)
        {

          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Check type */

              if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
                {
                  /* Fission, Get pointer to distribution */

                  yld = (long)RDB[yld + FISSION_YIELD_PTR_DISTR];
                  CheckPointer(FUNCTION_NAME, "yld", DATA_ARRAY, yld);

                  /* Loop over distribution */

                  while (yld > VALID_PTR)
                    {
                      /* Check pointer to new S(a,b) nuclide (NOTE: yield exists */
                      /* if data file is given, but pointers are not linked if   */
                      /* not in burnup mode). */

                      if ((tgt = (long)RDB[yld + FY_PTR_TGT]) > VALID_PTR)
                        if ((new = (long)RDB[tgt + NUCLIDE_PTR_SAB_NUC]) >
                            VALID_PTR)
                          WDB[yld + FY_PTR_TGT] = (double)new;

                      /* Next */

                      yld = NextItem(yld);
                    }
                }
              else if ((tgt = (long)RDB[rea + REACTION_PTR_TGT]) > VALID_PTR)
                {
                  /* Decay or transmutation, check pointer to new S(a,b) nuclide */

                  if ((new = (long)RDB[tgt + NUCLIDE_PTR_SAB_NUC]) > VALID_PTR)
                    WDB[rea + REACTION_PTR_TGT] = (double)new;

                }

              /* Next reaction */

              rea = NextItem(rea);
            }

          /* Next nuclide in composition */

          iso = NextItem(iso);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /* Update counters */

  ReactionCount();
}

/*****************************************************************************/
