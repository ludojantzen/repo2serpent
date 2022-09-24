/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : formtransmupathsc.c                            */
/*                                                                           */
/* Created:       2013/01/24 (JLe)                                           */
/* Last modified: 2020/04/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Recursive routine that generates decay and transmutation     */
/*              paths and adds missing data.                                 */
/*                                                                           */
/* Comments: - muuta nimi looptransmupathiksi                                */
/*                                                                           */
/*           - Tämä on yleinen silmukkarutiini transmutaatiopolkujen yli,    */
/*             ja sillä on 4 eri moodia:                                     */
/*                                                                           */
/*             1 - Lisätään lähtönuklidista syntyvien ketjujen data, jolloin */
/*                 katkaisu tehdään pituuden l mukaan                        */
/*                                                                           */
/*             2 - Muodostetaan polut, asetetaan pointterit ja lisätään      */
/*                 jäljellä oleva data, varataan jotain muistia, jne...      */
/*                                                                           */
/*             3 - Asetetaan nuklidien lämpötilat TMS:ää varten              */
/*                                                                           */
/*           - Jos palamalaskentamoodi on päällä, niin toi luuppaa noita     */
/*             ketjuja sellaisistakin materiaaleista mitä ei polteta (tää    */
/*             vaikuttaa mm. siihen miten nuklidien tms-minimi ja -maksimi   */
/*             lämpötila muodostuu). Ei varmaan (?) aiheuta virhettä, mutta  */
/*             toi voi joissain tilanteissa olla turhaa.                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FormTransmuPaths:"

/*****************************************************************************/

void FormTransmuPaths(long nuc, long l, double Tmin, double Tmax, long mode,
                      long TMS)
{
  long rea, ZAI, ZA, I, new, mt, yld, ptr, ba;
  double T;
  char lib[MAX_STR];

  /* Check pointer */

  if (nuc < VALID_PTR)
    Die(FUNCTION_NAME, "Pointer error");

  /* Check if nuclide is already processed and set flag */

  if ((long)RDB[nuc + NUCLIDE_OPTIONS] & OPT_USED)
    return;
  else if ((mode == 1) && (l > 5))
    return;
  else
    SetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

  /* Check chain length */

  if (l > 1000000)
    Die(FUNCTION_NAME, "Endless chain");

  /* Allocate memory for matrix index */

  if (mode == 2)
    AllocValuePair(nuc + NUCLIDE_PTR_MATRIX_IDX);

  /* Put TMS minimum and maximum temperatures */

  if (mode == 3)
    {
      /* Check temperatures */

      CheckValue(FUNCTION_NAME, "Tmin", "", Tmin, 0.0, INFTY);
      CheckValue(FUNCTION_NAME, "Tmax", "", Tmax, Tmin, INFTY);

      /* Check TMS flag */

      if (TMS == NO)
        Die(FUNCTION_NAME, "TMS flag not set");

      /* Put minimum temperature */

      if (Tmin < RDB[nuc + NUCLIDE_TMS_MIN_TEMP])
        WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = Tmin;

      /* Put maximum temperature */

      if (Tmax > RDB[nuc + NUCLIDE_TMS_MAX_TEMP])
        WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = Tmax;
    }

  /* Get temperature and library id */

  T = RDB[nuc + NUCLIDE_TEMP];
  sprintf(lib, "%s", GetText(nuc + NUCLIDE_PTR_LIB_ID));

  /* Loop over reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Get reaction mt */

      mt = (long)RDB[rea + REACTION_MT];

      /* Check reaction */

      if (UserTransmuCut(nuc, rea) == YES)
        {
          /* Cut-off based on user-defined list, pointer to next */

          rea = NextItem(rea);

          /* Cycle loop */

          continue;
        }
      else if (((ZAI = ReactionTargetZAI(rea)) == 0) ||
          (RDB[rea + REACTION_EMIN] >= RDB[DATA_NEUTRON_EMAX]))
        {
          /* Energy cut-off, pointer to next */

          rea = NextItem(rea);

          /* Cycle loop */

          continue;
        }
      else if (ZAI > 0)
        {
          /********************************************************************/

          /***** Decay and transmutation reactions ****************************/

          /* Check pointer to self */

          if (ZAI == (long)RDB[nuc + NUCLIDE_ZAI])
            Die(FUNCTION_NAME, "Path to self (%s %ld)",
                GetText(nuc + NUCLIDE_PTR_NAME), mt);

          /* Separate ZA and I */

          ZA = (long)((double)ZAI/10.0);
          I = ZAI - 10*ZA;

          /* Reset ba flag */

          ba = 0;

          /* Retry loop */

          do
            {
              /* Set target ZAI */

              WDB[rea + REACTION_TGT_ZAI] = (double)ZAI;

              /* Loop over existing nuclides and find match */

              new = (long)RDB[DATA_PTR_NUC0];
              while (new > VALID_PTR)
                {
                  /* Compare ZAI, T, TMS and lib */

                  if (ZAI == (long)RDB[new + NUCLIDE_ZAI])
                  if (T == RDB[new + NUCLIDE_TEMP])
                  if (TMS == ((long)RDB[new + NUCLIDE_TYPE_FLAGS] &
                              NUCLIDE_FLAG_TMS))
                  if (CompareStr(nuc + NUCLIDE_PTR_LIB_ID,
                                 new + NUCLIDE_PTR_LIB_ID) == YES)
                    {
                      /* Set daughters flag (NOTE: tämä tarvitaan siksi,  */
                      /* että nuklidi on voitu lisätä jostain sellaisesta */
                      /* materiaalista, joka ei osallistu palamalaskuun.  */
                      /* Jos flagia ei aseteta, ketju ei jatku eteenpäin. */

                      if ((mode == 1) || (mode == 2))
                        SetOption(new + NUCLIDE_TYPE_FLAGS,
                                  NUCLIDE_FLAG_NEW_DAUGHTERS);

                      /* Break loop */

                      break;
                    }

                  /* Next */

                  new = NextItem(new);
                }

              /* TODO: Tähän joku järkevä katkaisu energian ja ketjun */
              /* pituuden suhteen */

              /* Create new if not found */

              if ((new < VALID_PTR) && ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
                                        NUCLIDE_FLAG_NEW_DAUGHTERS))
                {
                  /* Check mode */

                  if (mode == 2)
                    {
                      /* Check type flag and add transport nuclide if */
                      /* fission product, otherwise decay nuclide.    */
                      /* NOTE: This was added 5.7.2019 / 2.1.32 as it */
                      /* was discovered that JEFF-3.3 fp data will    */
                      /* otherwise miss the cross sections of some    */
                      /* important neutron absorbers because the      */
                      /* direct yields are zero. */

                      if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
                          NUCLIDE_FLAG_FP)
                        {
                          if ((new = AddNuclide(NULL, ZAI, lib, T,
                                                NUCLIDE_TYPE_TRANSPORT, TMS))
                              < VALID_PTR)
                            new = AddNuclide(NULL, ZAI, lib, T,
                                             NUCLIDE_TYPE_DECAY, TMS);
                        }
                      else
                        {
                          /* NOTE: tätäkin pitää ehkä miettiä... */

                          new = AddNuclide(NULL, ZAI, lib, T,
                                           NUCLIDE_TYPE_DECAY, TMS);
                        }
                    }
                  else if (mode == 1)
                    {
                      /* Add transport or decay nuclide */

                      if ((new = AddNuclide(NULL, ZAI, lib, T,
                                            NUCLIDE_TYPE_TRANSPORT, TMS))
                          < VALID_PTR)
                        new = AddNuclide(NULL, ZAI, lib, T, NUCLIDE_TYPE_DECAY,
                                         TMS);
                    }
                }

              /* Check pointer */

              if (new > VALID_PTR)
                {
                  /* Check mode */

                  if (mode == 2)
                    {
                      /* Put pointer */

                      WDB[rea + REACTION_PTR_TGT] = (double)new;

                      /* Set decay, activation product or branch product flag */

                      if (mt > 20000)
                        SetOption(new + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_BP);
                      else if (mt > 10000)
                        SetOption(new + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DP);
                      else
                        SetOption(new + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_AP);

                      /* Set depletion flag */

                      SetOption(new + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DEP);
                    }

                  /* Call recursively */

                  FormTransmuPaths(new, l + 1, Tmin, Tmax, mode, TMS);
                }

              /* Check if found */

              if (new < VALID_PTR)
                {
                  /* Check if ba nuclide */

                  if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_BA)
                      && (ba == 0))
                    {
                      /* Try main chain */

                      ZAI = ZAI - 1000;

                      /* Put flag to avoid second call */

                      ba = 1;
                    }

                  /* Try with lower excited state */

                  else
                    {
                      I = I - 1;
                      ZAI = ZAI - 1;
                    }
                }
            }
          while ((new < VALID_PTR) && (I > -1));

          /* Check pointer and mode */

          if ((new < VALID_PTR) && (mode == 2))
            {
              /* Put pointer to lost data */

              WDB[rea + REACTION_PTR_TGT] = RDB[DATA_PTR_NUCLIDE_LOST];

              /* Add to count */

              WDB[DATA_N_DEAD_PATH] = RDB[DATA_N_DEAD_PATH] + 1.0;
            }

          /********************************************************************/
        }
      else if (mode != 1)
        {
          /********************************************************************/

          /***** Spontaneous and neutron-induced fission **********************/

          /* Set target ZAI (tää tarvitaan koska myöhemmin kutsuttavat  */
          /* rutiinit tunnistaa fissiot siitä että TGT_ZAI < 0, mikä on */
          /* kieltämättä vähän typerää). */

          if (mode == 2)
            WDB[rea + REACTION_TGT_ZAI] = (double)ZAI;

          /* Get pointer to yield */

          if ((yld = (long)RDB[rea +  REACTION_PTR_FISSY]) > VALID_PTR)
            {
              /* Check mt */

              if (((mt < 18) || (mt > 21)) && (mt != 38) && (mt != 10006))
                Die(FUNCTION_NAME, "Invalid fission mt %ld", mt);

              /* Get pointer to distribution */

              if ((yld = (long)RDB[yld + FISSION_YIELD_PTR_DISTR]) < VALID_PTR)
                Die(FUNCTION_NAME, "Pointer error");

              /* Loop over yield (missä cutoffi?) */

              while (yld > VALID_PTR)
                {
                  /* Get ZAI */

                  ZAI = (long)RDB[yld + FY_TGT_ZAI];

                  /* Separate ZA and I */

                  ZA = (long)((double)ZAI/10.0);
                  I = ZAI - 10*ZA;

                  /* Retry loop */

                  do
                    {
                      /* Loop over existing nuclides and find match */

                      new = (long)RDB[DATA_PTR_NUC0];
                      while (new > VALID_PTR)
                        {
                          /* Compare ZAI, T, TMS and lib */

                          if (ZAI == (long)RDB[new + NUCLIDE_ZAI])
                          if (T == RDB[new + NUCLIDE_TEMP])
                          if (T == RDB[new + NUCLIDE_TEMP])
                          if (TMS == ((long)RDB[new + NUCLIDE_TYPE_FLAGS] &
                                      NUCLIDE_FLAG_TMS))
                          if (CompareStr(nuc + NUCLIDE_PTR_LIB_ID,
                                         new + NUCLIDE_PTR_LIB_ID) == YES)
                            break;

                          /* Next */

                          new = NextItem(new);
                        }

                      /* Create new if not found */

                      if ((new < VALID_PTR) && (mode == 2) &&
                          ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
                           NUCLIDE_FLAG_NEW_DAUGHTERS))
                        {
                          /* Create nuclide */

                          new = AddNuclide(NULL, ZAI, lib, T,
                                           NUCLIDE_TYPE_DECAY, TMS);
                        }

                      /* Check pointer */

                      if (new > VALID_PTR)
                        {
                          /* Check mode */

                          if (mode == 2)
                            {
                              /* Put pointer */

                              WDB[yld + FY_PTR_TGT] = (double)new;

                              /* Set fission product flag and depletion flags */

                              SetOption(new + NUCLIDE_TYPE_FLAGS,
                                        NUCLIDE_FLAG_FP);
                              SetOption(new + NUCLIDE_TYPE_FLAGS,
                                        NUCLIDE_FLAG_DEP);
                            }

                          /* Call recursively */

                          FormTransmuPaths(new, l + 1, Tmin, Tmax, mode, TMS);
                        }

                      /* If not found, try with lower excited state */

                      if (new < VALID_PTR)
                        {
                          I = I - 1;
                          ZAI = ZAI - 1;
                        }
                    }
                  while ((new < VALID_PTR) && (I > -1));

                  /* Check pointer and mode */

                  if ((new < VALID_PTR) && (mode == 2))
                    {
                      /* Put pointer to lost data */

                      WDB[yld + FY_PTR_TGT] = RDB[DATA_PTR_NUCLIDE_LOST];

                      /* Get pointer to nuclide-wise yield */

                      ptr = (long)RDB[nuc + NUCLIDE_PTR_NFY_DATA];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                      /* Check total or first-chance fission and that */
                      /* yield is own. */

                      if (((mt == 18) || (mt == 19)) &&
                          ((long)RDB[ptr + FISSION_YIELD_PARENT_ZAI] ==
                           (long)RDB[nuc + NUCLIDE_ZAI]))
                        WDB[DATA_N_DEAD_PATH] = RDB[DATA_N_DEAD_PATH] + 1.0;
                    }

                  /* Next */

                  yld = NextItem(yld);
                }
            }
          else if (ZAI == -2)
            {
              /* Spontaneous fission without yield */

              WDB[rea + REACTION_PTR_TGT] = RDB[DATA_PTR_NUCLIDE_LOST];

              /* Update count */

              WDB[DATA_N_DEAD_PATH] = RDB[DATA_N_DEAD_PATH] + 1.0;
            }
          else if (((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_DECAY)
                   && ((long)RDB[DATA_PTR_NFYDATA_FNAME_LIST] > 0))
            Die(FUNCTION_NAME, "No fission yield distribution %s mt %ld",
                GetText(nuc + NUCLIDE_PTR_NAME), mt);

          /********************************************************************/
        }

      /* Next reaction */

      rea = NextItem(rea);
    }

  /* Link total fission to lost */

  if ((rea = (long)RDB[nuc + NUCLIDE_PTR_TOTFISS_REA]) > VALID_PTR)
    {
      /* Put pointer to lost data */

      WDB[rea + REACTION_PTR_TGT] = RDB[DATA_PTR_NUCLIDE_LOST];
    }
}

/*****************************************************************************/
