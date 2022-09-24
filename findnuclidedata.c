/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findnuclidedata.c                              */
/*                                                                           */
/* Created:       2010/09/26 (JLe)                                           */
/* Last modified: 2020/01/15 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Finds nuclide from list of existing or from data             */
/*                                                                           */
/* Comments: - Logiikka: 1) Alkuperäisessä koostumuksessa oleva nuklidi      */
/*                          identifioidaan vain zaid:n perusteella.          */
/*                       2) Jos zaid:a ei ole määritelty, järjestys on       */
/*                          seuraava:                                        */
/*                                                                           */
/*                          a) - ZAI, lib, T                                 */
/*                          b) - ZAI, T                                      */
/*                          c) - ZAI, lib (tää ja toi vika ei enää oo)       */
/*                          d) - ZAI                                         */
/*                                                                           */
/*                          Jos nuklidityypit on listassa oikeassa järjes-   */
/*                          tyksessä (transport - decay - struct), mitään    */
/*                          erillistä tyyppitarkistusta ei tarvita.          */
/*                                                                           */
/*                       3) Tärkeää: jos täydellistä matchia ei löydy, ei    */
/*                          käytetä nuklidia jolle on määritetty joko lib    */
/*                          tai T, vaan luodaan uusi.                        */
/*                                                                           */
/*                       4) Type-argumentti pakottaa hakemaan vain tietyn    */
/*                          tyyppistä dataa.                                 */
/*                                                                           */
/*           - Return value: > 0 : nuclide found in list of existing         */
/*                           < 0 : pointer to data in ACE array              */
/*                           = 0 : no data found                             */
/*                                                                           */
/*           - Tää rutiini on alkuprosessoinnin pullonkaula. ACE-arrayn      */
/*             voisi sortata ja hakua tehdä tehokkaammaksi (ei välttämättä   */
/*             kyllä vaikuta).                                               */
/*                                                                           */
/*           - Tää hakee tarpeettomasti uuden datan jokaiselle eri           */
/*             lämpötilassa olevalle nuklidille kun on-the-fly -menetelmä    */
/*             on käytössä.                                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindNuclideData:"

#define MAX_ISO 1000

/*****************************************************************************/

long FindNuclideData(char *zaid, long ZAI, char *lib, double T, long type,
                     long TMS)
{
  long nuc, ace, zptr[MAX_ISO], n;

  /* Check if zaid is given */

  if (zaid != NULL)
    {
      /**********************************************************************/

      /***** Initial composition ********************************************/

      /* Check if temperature is given (Doppler-broadening is in use) */

      if (T > 0.0)
        {
          /* Check that TMS flag is not set */

          if (TMS != NO)
            Die(FUNCTION_NAME, "TMS flag is set");

          /* Loop over existing nuclides */

          nuc = (long)RDB[DATA_PTR_NUC0];
          while (nuc > VALID_PTR)
            {
              /* Compare zaid, nuclide temperature and TMS flag */

              if (!strcmp(zaid, GetText(nuc + NUCLIDE_PTR_NAME)) &&
                  (RDB[nuc + NUCLIDE_TEMP] == T) &&
                  (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS)
                   == TMS))
                return nuc;

              /* Next nuclide */

              nuc = NextItem(nuc);
            }
        }
      else
        {
          /* Loop over existing nuclides */

          nuc = (long)RDB[DATA_PTR_NUC0];
          while (nuc > VALID_PTR)
            {
              /* Compare (XS temperature must be the same as nuclide */
              /* temperature) */

              if ((!strcmp(zaid, GetText(nuc + NUCLIDE_PTR_NAME))) &&
                  (RDB[nuc + NUCLIDE_TEMP] == RDB[nuc + NUCLIDE_XS_TEMP]) &&
                  (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS)
                   == TMS))
                return nuc;

              /* Next nuclide */

              nuc = NextItem(nuc);
            }
        }

      /* Loop over ace data */

      ace = (long)RDB[DATA_PTR_ACE0];
      while (ace > VALID_PTR)
        {
          /* NOTE: Tässä korvattiin alias namella, koska processnuclides.c:n */
          /* alussa COMPOSITION:iin tallennetaan NAME (2.1.22 / 12.7.2014) */

          WDB[DATA_DUMMY] = ACE[ace + ACE_PTR_NAME];
          if (!strcmp(GetText(DATA_DUMMY), zaid))
            return -ace;

          /* NOTE: Myös alias pitää ehkä tarkistaa S(a,b)-dataa varten? */
          /*
          WDB[DATA_DUMMY] = ACE[ace + ACE_PTR_ALIAS];
          if (!strcmp(GetText(DATA_DUMMY), zaid))
            return -ace;
          */
          /* next */

          ace = (long)ACE[ace + ACE_PTR_NEXT];
        }

      /* Not found */

      return 0;

      /***********************************************************************/
    }
  else if (type == NUCLIDE_TYPE_TRANSMUXS)
    {
      /***********************************************************************/

      /***** Supplementary transmutation cross section ***********************/

      /* Check ZAI */

      if ((ZAI < 10000) || (ZAI > 1200000))
        Die(FUNCTION_NAME, "Invalid ZAI %ld", ZAI);

      /* Loop over ace data */

      ace = (long)RDB[DATA_PTR_TRANSMU_ACE0];
      while (ace > VALID_PTR)
        {
          /* Compare */

          if (((long)ACE[ace + ACE_TYPE] != NUCLIDE_TYPE_TRANSMUXS) ||
              ((long)ACE[ace + ACE_ZAI] == ZAI))
            break;

          /* next */

          ace = (long)ACE[ace + ACE_PTR_NEXT];
        }

      /* Check if found */

      if (ace < 1)
        return 0;
      else if ((long)ACE[ace + ACE_TYPE] != NUCLIDE_TYPE_TRANSMUXS)
        return 0;

      /* Check ZAI */

      if ((long)ACE[ace + ACE_ZAI] != ZAI)
        Die(FUNCTION_NAME, "Mismatch in ZAI");

      /* Return pointer */

      return -ace;
    }

  /***************************************************************************/

  /***** Check input parameters **********************************************/

  /* Check ZAI, lib and T */

  if ((ZAI < 10000) || (ZAI > 1200000))
    Die(FUNCTION_NAME, "Invalid ZAI %ld", ZAI);

  if ((lib == NULL) && (type != NUCLIDE_TYPE_PHOTON))
    Die(FUNCTION_NAME, "Library ID not given");

  /* NOTE: initial compositionin decay nuklideille lämpötila on väkisin < 0. */
  /*       jos tää tarkistus otetaan pois, break-lausekkeet pitää myös ottaa */
  /*       pois vikasta osiosta. */

  if (((T < 0.0) || (T > 1E+5)) && (type != NUCLIDE_TYPE_DECAY) && (type > 0))
    Die(FUNCTION_NAME, "Invalid T %E", T);

  /***************************************************************************/

  /***** Check existing data ************************************************/

  /* Find matching ZAI, lib and T in existing nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Compare */

      if ((long)RDB[nuc + NUCLIDE_ZAI] == ZAI)
        if ((type < 1) || ((long)RDB[nuc + NUCLIDE_TYPE] == type))
          if (RDB[nuc + NUCLIDE_TEMP] == T)
            if ((((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS)
                 == TMS))
              {
                /* Check lib (photon data can be added from readphotondata.c */
                /* without lib) */

                if (type == NUCLIDE_TYPE_PHOTON)
                  return nuc;
                else if (lib == NULL)
                  Die(FUNCTION_NAME, "Null lib pointer (ZAI = %ld)", ZAI);
                if (!strcmp(GetText(nuc + NUCLIDE_PTR_LIB_ID), lib))
                  return nuc;
              }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Uusi 3.4.2012 / 2.1.4 - match initial decay nuclides with ZAI only */

#ifdef mmmmmmmmmmmmmmmmmmmmmm

  if (type == NUCLIDE_TYPE_DECAY)
    {
      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Compare ZAI, type and check initial flag */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_INITIAL)
            if ((long)RDB[nuc + NUCLIDE_ZAI] == ZAI)
              if ((long)RDB[nuc + NUCLIDE_TYPE] == type)
                {
                  /* Return pointer to nuclide */

                  return nuc;
                }

          /* Next nuclide */

          nuc = NextItem(nuc);
        }
    }

#endif

  /***************************************************************************/

  /***** Find nuclide in ACE data ********************************************/

  /* Pointer to beginning of data */

  if (type == NUCLIDE_TYPE_DECAY)
    ace = (long)RDB[DATA_PTR_DECAY_ACE0];
  else
    ace = (long)RDB[DATA_PTR_ACE0];

  /* Read nuclides with same ZAI to list */

  n = 0;
  while (ace > VALID_PTR)
    {
      /* Compare ZA and add to list */

      if (((long)ACE[ace + ACE_ZAI] == ZAI) &&
          ((long)ACE[ace + ACE_TYPE] != NUCLIDE_TYPE_TRANSMUXS))
        if ((type < 1) || (type == (long)ACE[ace + ACE_TYPE]))
          zptr[n++] = ace;

      /* Check size */

      if (n > MAX_ISO - 2)
        Die(FUNCTION_NAME, "Increase buffer size");

      /* Next */

      ace = (long)ACE[ace + ACE_PTR_NEXT];
    }

  /* Check if not found */

  if (n == 0)
    return 0;

  /* Put null pointer */

  zptr[n] = -1;

  /* Skip lib and T check if type is decay or photon */

  if ((type != NUCLIDE_TYPE_DECAY) && (type != NUCLIDE_TYPE_PHOTON))
    {
      /* Check lib */

      if (lib == NULL)
        Die(FUNCTION_NAME, "Null lib pointer (ZAI = %ld)", ZAI);

      /* lib and T */

      n = 0;
      while ((ace = zptr[n++]) > VALID_PTR)
        {
          /* Compare T */

          if (ACE[ace + ACE_TEMP] == T)
            {
              /* Library ID pointer to dummy variable */

              WDB[DATA_DUMMY] = ACE[ace + ACE_PTR_LIB_ID];

              /* Compare lib */

              if (!strcmp(GetText(DATA_DUMMY), lib))
                return -ace;
            }
        }

      /* lib */

      n = 0;
      while ((ace = zptr[n++]) > VALID_PTR)
        {
          /* Library ID pointer to dummy variable */

          WDB[DATA_DUMMY] = ACE[ace + ACE_PTR_LIB_ID];

          /* Compare lib */

          if (!strcmp(GetText(DATA_DUMMY), lib))
            return -ace;
        }

      /* T */

      n = 0;
      while ((ace = zptr[n++]) > VALID_PTR)
        {
          /* Compare T */

          if (ACE[ace + ACE_TEMP] == T)
            return -ace;
        }

      /* Noiden muiden kanssa voi tulla ongelmia lämpötilamoodeissa */
      /* (25.2.2012 / 2.1.13) */

      return 0;
    }

  /* Return first in list */

  return -zptr[0];
}

/*****************************************************************************/
