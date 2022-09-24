/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dopmicroxs.c                                   */
/*                                                                           */
/* Created:       2011/09/06 (JLe)                                           */
/* Last modified: 2019/12/18 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Returns temperature-corrected microscopic cross section      */
/*                                                                           */
/* Comments: - Lämpötilan saa vaikka materiaalista, ja mietitään myöhemmin   */
/*             koordinaatteihin perustuva haku (ei toimi nyt kun trackin-    */
/*             rutiini käyttää sekaisin delta- ja surface-trackingiä)        */
/*                                                                           */
/*           - Arvottu suhteellinen nopeus pitänee ottaa talteen ja käyttää  */
/*             samaa arvoa jos samalle nuklidille haetaan saman törmäyksen   */
/*             jälkeen samaa tai jotain toista vaikutusalaa (mietitään       */
/*             tarkemmin sitten kun saadaan jotain tuloksia)                 */
/*                                                                           */
/*           - Tää homma tulee todennäköisesti toimimaan vain optimointi-    */
/*             moodeissa 0, 1 ja 2 --> lisätään joku tarkistus alkuun?       */
/*                                                                           */
/*           - Toi id on OpenMP:n thread id, ja se välitetään melkein        */
/*             kaikille aliohjelmille                                        */
/*                                                                           */
/*           - Doppler-preprosessorirutiini pitää kytkeä pois niin kauan     */
/*             lämpötila otetaan MATERIAL_TEMP:stä.                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DopMicroXS:"

/*****************************************************************************/

double DopMicroXS(long mat, long rea, double E, double *Er, double T, long id)
{
  long nuc, ptr, loc0, ncol;
  double awr, xs, kT, g;
  double ycn, r1, z, z2, rnd1, rnd2, s,c, ar,x2;
  double T0, dT;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to nuclide data */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Check TMS flag */

  if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS))
    Die(FUNCTION_NAME, "Nuclide %s not flagged for TMS",
        GetText(nuc + NUCLIDE_PTR_NAME));

  /* If E > lower URES boundary of the nuclide, return 0 K cross section.    */
  /* Tee samoin jos thresholdireaktio. Myoskaan NJOY ei levenna vaikutualoja */
  /* URES-rajan ylapuolella. Jos ihan tarkkoja ollaan, NJOY lopettaa         */
  /* leventamisen kun ensimmainen seuraavista tayttyy:                       */
  /* - Ylitetaan URES-raja.                                                  */
  /* - Ylitetaan matalaenergimman threshold-reaktion kynnysenergia           */
  /* - energia > 1 MeV (tama on NJOY:n THNMAX-inputparametrin vakioarvo,     */
  /*   eli energia voi periaatteessa olla jokin muukin...)                   */

  /* Toi ensimmäinen testi vähän epäilyttää... Olisiko tuo parempi katsoa */
  /* REACTION_EMIN:n perusteella? Vai meneekö jotain sekaisin jos data ei */
  /* ala nollasta? */

  /* NOTE: Joillain luonnon materiaaleilla toi ures-raja voi olla jotain */
  /*       ihan häröä (esim. Zr-nat), minkä vuoksi tuo if-haara toteutuu */
  /*       vaikka ei pitäisi. Tämä sotkee koko reaktiosämpläyksen. Nyt   */
  /*       tuo on korjattu lisäämällä ylimääräinen testi että used-flagi */
  /*       on päällä. (2.1.32 / 18.12.2019 / JLe). */

  if ((RDB[rea + REACTION_XS_I0] > 0) || (E > 1.00) ||
      ((E > RDB[nuc + NUCLIDE_URES_EMIN]) &&
       ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_USED)))
    {
      /* Check for temperature sensitivity flag */

      if ((loc0 = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR)
        if (((long)RDB[loc0 + SENS_PERT_FLAGS] & SENS_PERT_FLAG_TEMPERATURE) > 0)
          {
            /* Get collision number */

            ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
            ncol = (long)GetPrivateData(ptr, id);

            /* Store an unphysical target energy */

            StoreValuePair(nuc + NUCLIDE_PREV_COL_ET, (double)ncol, 1e10, id);
          }

      *Er=E;
      return MicroXS(rea, E, id);
    }

  /* Check OTF S(a,b) treatment */

  if((long)RDB[rea + REACTION_MT] == 2004 ||
     (long)RDB[rea + REACTION_MT] == 2002 ) {

    return OTFSabXS(rea, E, T, id);
  }

  /* Compare to cross section temperature */

  if (T < RDB[nuc + NUCLIDE_TMS_MIN_TEMP])
    Die(FUNCTION_NAME, "Adjusted temperature below minimum (%s %s %E %E)",
        GetText(mat + MATERIAL_PTR_NAME), GetText(nuc + NUCLIDE_PTR_NAME),
        T, RDB[nuc + NUCLIDE_TMS_MIN_TEMP]);
  else if (T > RDB[nuc + NUCLIDE_TMS_MAX_TEMP])
    Die(FUNCTION_NAME, "Adjusted temperature above maximum (%s %s %E %E)",
        GetText(mat + MATERIAL_PTR_NAME), GetText(nuc + NUCLIDE_PTR_NAME),
        T, RDB[nuc + NUCLIDE_TMS_MAX_TEMP]);

  /* Adjust lower boundary */

  kT = T - RDB[nuc + NUCLIDE_TMS_MIN_TEMP];

  /* Remember these values so that they can be stored if energy is sampled */

  T0 = T;
  dT = kT;

  /***************************************************************************/
  /* Determine relative (target-at-rest) energy Er                         ***/
  /***************************************************************************/

  /* Check zero change (skip sampling if zero) */

  if (kT == 0.0)
    {

      /* Check for temperature sensitivity flag */

      if ((loc0 = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR)
        if (((long)RDB[loc0 + SENS_PERT_FLAGS] & SENS_PERT_FLAG_TEMPERATURE) > 0)
        {
          /* Get collision number */

          ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
          ncol = (long)GetPrivateData(ptr, id);

          /* Store zero target energy  */

          StoreValuePair(nuc + NUCLIDE_PREV_COL_ET, (double)ncol, 0.0, id);
        }

      *Er=E;
    }
  else {

    /* Get collision number */

    ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
    ncol = (long)GetPrivateData(ptr, id);

    /* Check previously stored relative energy */

    if ((*Er = TestValuePair(nuc + NUCLIDE_PREV_COL_ER, (double)ncol, id))
        < 0.0)
      {

        /* Convert to MeV */

        kT = kT*KELVIN;

        /* Get atomic weight ratio */

        awr = RDB[nuc + NUCLIDE_AWR];

        /*********************************************************************/
        /* Sämplää kohtioytimen nopeus kuten vapaakaasumallissa            ***/

        /* Lainattu vissiin mcnpstä */

        ar = awr/kT;
        ycn = sqrt(E*ar);

        do
          {
            if (RandF(id)*(ycn + 1.12837917) > ycn)
              {
                r1 = RandF(id);
                z2 = -log(r1*RandF(id));
              }
            else
              {
                do
                  {
                    rnd1 = RandF(id);
                    rnd2 = RandF(id);

                    r1 = rnd1*rnd1;
                    s = r1 + rnd2*rnd2;
                  }
                while (s > 1.0);

                z2 = -r1*log(s)/s - log(RandF(id));
              }

            z = sqrt(z2);
            c = 2.0*RandF(id) - 1.0;

            x2 = ycn*ycn + z2 - 2*ycn*z*c;

            rnd1 = RandF(id)*(ycn + z);
          }
        while (rnd1*rnd1 > x2);

        /* x2 on suhteellisen nopeuden neliö kertaa 1. NSE-paperissa     */
        /* mainittu gamma toiseen, z2 on targetin nopeuden neliö kertaa  */
        /* gamma toiseen ja c on suuntakosini verrattuna alkup neutronin */
        /* suuntaan. Target energy in MeV is z2/ar*awr */

        /* Kiitos näppärien symbolivalintojen, seuraava rivi laskee */
        /* suhteellista nopeutta vastaavan neutronin energian */

        *Er=x2/ar;

        /* Check for temperature sensitivity flag */

        if ((loc0 = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR)
          if (((long)RDB[loc0 + SENS_PERT_FLAGS] & SENS_PERT_FLAG_TEMPERATURE) > 0)
          {
            /* Store the target energy (MeV) */

            StoreValuePair(nuc + NUCLIDE_PREV_COL_ET, (double)ncol, z2/ar*awr, id);

            /* Store the collision temperature (absolute) */

            StoreValuePair(nuc + NUCLIDE_PREV_COL_T, (double)ncol,  T0, id);

            /* Store the collision temperature (relative to base XS) */

            StoreValuePair(nuc + NUCLIDE_PREV_COL_DT, (double)ncol, dT, id);
          }

        /* Tallennetaan arvotut targetin nopeuden nelio (kertaa gamma) */
        /* ja kosini jotta voidaan hyodyntaa samaa nopeutta elastisen  */
        /* sironnan kanssa  */

        /* Tallennettava z2 pitää vielä jakaa ar:llä (JLE 2012/10/26) */

        if (RDB[nuc + NUCLIDE_XS_TEMP] == 0)
          {
            StoreValuePair(nuc + NUCLIDE_PREV_COL_Z2, (double)ncol, z2/ar, id);
            StoreValuePair(nuc + NUCLIDE_PREV_COL_COS, (double)ncol, c, id);
          }

      }

    /* Tallennetaan suhteellinen energia */

    StoreValuePair(nuc + NUCLIDE_PREV_COL_ER, (double)ncol, *Er, id);

  }

  /* Noudetaan suht. nopeutta vastaava mikroskooppinen vaikutusala */

  xs = MicroXS(rea, *Er, id);

  /* Tehdään nyt PotCorr -korjaus jo tässä vaiheessa (TVi 2015-04-15)       */
  /* Monimutkaistaa asioita reaktiosämpläysvaiheessa, mutta yksinkertaistaa */
  /* myöhemin */

  /* Pitää käyttää alkuperäistä lampotilaa koska potcorr miinustaa */
  /* minimilämpötilan */

  g = PotCorr(nuc, E, T*KELVIN);
  xs = xs*g;

  /* If OTF S(a,b) treatment is in use, subtract contribution of */
  /* free atom elastic scattering and add bound atom scattering  */
  /* (for total and elastic xs) */

  /* Kannattaa huomata, että tätä S(a,b) -sironnan osuutta
     ei pidä kertoa potcorr-kertoimella */

  if(E < RDB[nuc + NUCLIDE_SAB_EMAX] &&
     ((long)RDB[rea + REACTION_MT] == 1 || (long)RDB[rea + REACTION_MT] == 2))
    {

      /* Energy region in which only one S(a,b) reaction is active
         [EMAXLOW, EMAX] needs a separate treatment. If E is outside this
         region, simply subtract elastic scattering and add S(a,b)
         reactions for mt=1 */

      if(RDB[nuc + NUCLIDE_SAB_EMAXLOW] == 0.0 ||
         E < RDB[nuc + NUCLIDE_SAB_EMAXLOW]){

        xs = xs - g*MicroXS((long)RDB[nuc + NUCLIDE_PTR_ELAXS], *Er, id);

        if((long)RDB[rea + REACTION_MT] == 1)
          xs = xs + OTFSabXS(rea, E, T, id);
      }

      /* Within [EMAXLOW, EMAX], subtract active S(a,b) from scattering xs for mt=2 */
      /* For mt=1 xs = xs - OTFSabXS(rea, E, T, id) + OTFSabXS(rea, E, T, id) = xs
         between EMAXLOW and EMAX */

      else if ( (long)RDB[rea + REACTION_MT] == 2){
        xs = xs - OTFSabXS(rea, E, T, id);
      }

    }

  /* Check */

  CheckValue(FUNCTION_NAME, "xs", "", xs, 0.0, INFTY);

  return xs;
}

/*****************************************************************************/
