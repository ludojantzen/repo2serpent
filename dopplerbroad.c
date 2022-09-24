/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dopplerbroad.c                                 */
/*                                                                           */
/* Created:       2011/07/28 (JLe)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Perfrorms Doppler-Broadening on cross sections               */
/*                                                                           */
/* Comments: - From Serpent 1.1.15                                           */
/*           - Original routine developed by Tuomas Viitanen                 */
/*           - Toteutus on tehty eri tavalla kuin Serpent1:ss‰, johtuen      */
/*             siit‰ ett‰ aliohjelmaa kutsuttaessa palamalaskun lis‰‰m‰t     */
/*             nuklidit eiv‰t ole osa materiaalikoostumuksia. Luuppi tehd‰‰n */
/*             nyt suoraan nuklidien yli.                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DopplerBroad:"

#define STD_LIMITS_SIZE 70000

/* Apufunktio find_limits lˆytyy DopplerBroad() -funktion j‰lkeen. */

long *find_limits(long, long, double, double, double, double *, long *, long);

/*****************************************************************************/

void DopplerBroad()
{
  long nuc, ace, ptr, loc0, n, nr, mt, uplimit, i;
  long NXS[16], JXS[32], NTR, NES, L0, L1, L, *limits, UNRES;
  double *XSS, T, Td, awr;

  /* Check flag */

  if ((long)RDB[DATA_USE_DOPPLER_PREPROCESSOR] == NO)
    return;

  fprintf(outp, "Running Doppler-broadening preprocessor:\n\n");

  limits = NULL;

  /* T‰m‰ pelk‰st‰‰n k‰‰nt‰j‰n varoituksen v‰ltt‰miseksi */

  UNRES = 0;

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check type */

      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
        {
          /* Get Doppler-broadened temperature in Kelvin */

          Td = RDB[nuc + NUCLIDE_TEMP];

          /* Get original temperature in Kelvin */

          T = RDB[nuc + NUCLIDE_XS_TEMP];

          /* If T = Td, skip nuclide */

          if(fabs(Td - T) < 1e-3)
            {
              /* Pointer to next */

              nuc = NextItem(nuc);

              /* Cycle loop */

              continue;
            }

          /* TODO: lis‰‰ johonkin tarkistus materiaalin l‰mpˆtilasta */

          else if (Td < T)
            Error(0, "Doppler temperature %1.1fK below original %1.1fK in %s",
                  Td, T, GetText(nuc + NUCLIDE_PTR_NAME));

          /* Print nuclide name */

          fprintf(outp,
                  "Adjusting nuclide %10s temperature from %1.0fK to %1.0fK...\n",
                  GetText(nuc + NUCLIDE_PTR_NAME), T, Td);

          /* Get atomic weight ratio */

          awr = RDB[nuc + NUCLIDE_AWR];

          /*******************************************************************/

          /***** Get NXS and JXS arrays **************************************/

          /* Pointer to ace block */

          ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];
          CheckPointer(FUNCTION_NAME, "(ace)", ACE_ARRAY, ace);

          /* Read data to NXS array */

          ptr = (long)ACE[ace + ACE_PTR_NXS];
          CheckPointer(FUNCTION_NAME, "(nxs)", ACE_ARRAY, ptr);

          for (n = 0; n < 16; n++)
            NXS[n] = (long)ACE[ptr++];

          /* Read data to JXS array */

          ptr = (long)ACE[ace + ACE_PTR_JXS];
          CheckPointer(FUNCTION_NAME, "(jxs)", ACE_ARRAY, ptr);

          for (n = 0; n < 32; n++)
            JXS[n] = (long)ACE[ptr++];

          /********************************************************************/

          /***** Make copies of ACE and XSS ***********************************/

          /* Pointer to original ace block */

          ptr = (long)RDB[nuc + NUCLIDE_PTR_ACE];
          CheckPointer(FUNCTION_NAME, "(ace)", ACE_ARRAY, ptr);

          /* Make a copy */

          ace = ReallocMem(ACE_ARRAY, ACE_BLOCK_SIZE);
          memcpy(&ACE[ace], &ACE[ptr], ACE_BLOCK_SIZE*sizeof(double));

          /* Set pointer */

          WDB[nuc + NUCLIDE_PTR_ACE] = (double)ace;

          /* Pointer to XSS array */

          ptr = (long)ACE[ace + ACE_PTR_XSS];
          CheckPointer(FUNCTION_NAME, "(xss)", ACE_ARRAY, ptr);

          /* Make a copy */

          loc0 = ReallocMem(ACE_ARRAY, NXS[0]);
          memcpy(&ACE[loc0], &ACE[ptr], NXS[0]*sizeof(double));

          /* Set pointer */

          ACE[ace + ACE_PTR_XSS] = (double)loc0;

          /* Pointer to data */

          XSS = &ACE[loc0];

          /********************************************************************/

          /***** T‰st‰ alkaa datan varsinainen k‰sittely *********************/

          /* Get number of reactions */

          NTR = NXS[3] + 3;
          CheckValue(FUNCTION_NAME, "NTR", "", NTR, 3, 1000);

          /* Loop over reaction channels */

          for (nr = 0; nr < NTR; nr++)
            {
              /* Get pointers to energy and cross section arrays */

              if (nr < 3)
                {
                  /***********************************************************/

                  /***** Total, absorption and elastic xs ********************/

                  /* Number of energy points */

                  NES = NXS[2];

                  /* Pointer to data (table F-4, page F-13) */

                  L0 = JXS[0] - 1;

                  /* Set reaction mt and type */

                  L1 = L0 + (nr + 1)*NES;

                  if (nr == 0)
                    {
                      /* Total */

                      mt = 1;
                    }
                  else if (nr == 1)
                    {
                      /* Absorption (sum of mt 100...199) */

                      mt = -2;
                    }
                  else
                    {
                      /* Elastic */

                      mt = 2;
                    }

                  /***********************************************************/
                }
              else
                {
                  /***********************************************************/

                  /***** Reaction cross sections *****************************/

                  /* Get reaction MT (Table F-6, page F-15). */

                  mt = (long)XSS[JXS[2] - 1 + nr - 3];

                  /* Get pointer to SIG-block (Table F-10,page F-17) */

                  L = (long)XSS[JXS[5] - 1 + nr - 3] + JXS[6] - 1;

                  /* Get number of energy points */

                  NES = (long)XSS[L];

                  /* Pointer to energy array */

                  L0 = JXS[0] - 1 + NXS[2] - NES;

                  /* Pointer to XS data */

                  L1 = L + 1;

                  /***********************************************************/
                }

              /* Check values */

              CheckValue(FUNCTION_NAME, "NES", "", NES, 1, 1E+6);

              /* Ensimm‰isell‰ reaktiokanavalla etsit‰‰n unresolved-raja  */
              /* mik‰li todn‰ktaulut lˆytyv‰t sek‰ lasketaan integrointi- */
              /* rajat */

              if (mt == 1)
                {
                  /* Ton ures-datan olemassaolo kannattaa testata, */
                  /* vaikka toi ei todenn‰kˆisesti mit‰‰n ongelmia */
                  /* aiheutakaan (9.6.2011 / 1.1.15 JLe): */

                  if ((JXS[22] > 0) && (XSS[JXS[22] + 6 - 1] > 1E-8))
                    {
                      UNRES = 1;
                      uplimit = NES-1;

                      /* Haarukoidaan alue ensin 1/1024 -osaan alkuper‰isest‰ */

                      for (i = 0; i < 10; i++)
                        {
                          if (XSS[L0 + UNRES + (long)floor((double)(uplimit - UNRES)/2.0)] < XSS[JXS[22] + 6 - 1])
                            UNRES = (long)floor((double)(uplimit - UNRES)/2.0) + UNRES;
                          else
                            uplimit = (long)floor((double)(uplimit - UNRES)/2.0) + UNRES;
                        }

                      /* Ja sitten etsit‰‰n tarkka oikea kohta */

                      while (XSS[L0 + UNRES] < XSS[JXS[22] + 6 - 1])
                        UNRES++;
                    }

                  /* Jos prob. tableja ei ole, tyydyt‰‰n integroimaan koko */
                  /* energia-alueen yli. T‰m‰ aiheuttaa l‰hes merkitykset- */
                  /* tˆm‰n lis‰n laskenta-aikaan ja virheen vaikutusalaan  */
                  /* unresolved-rajalle. */

                  else
                    {
                      /* Changed row below NXS[2]-1 -> NXS[2] (TVi 1.9.2011) */

                      UNRES = NXS[2];
                    }

                  /* Integrointirajat limits-vektoriin. */

                  limits = find_limits(L0, NXS[2], T, Td, awr, XSS, limits, UNRES);
                }

              /* Tehd‰‰n Doppler-levennys. Levennet‰‰n ainoastaan */
              /* kynnysenergialtaan alle 10 keV reaktiot, joiden  */
              /* mt on joko 2 (elastinen sironta) tai v‰lill‰ 16-117 */

              if (((mt == 2) || ((mt >= 16) && (mt <= 117))) &&
                  (NXS[2] - NES < UNRES) && (UNRES - NXS[2] + NES - 1 > 0) &&
                  (XSS[L0] < 0.01))
                {
                  /*
                    if( (mt==2 || (mt >=16 && mt <= 117)) && NXS[2]-NES<UNRES && XSS[L0] < 0.01 ){
                  */


                  /* Limits-vektorista ainoastaan soveltuvat osat          */
                  /* BCS-funktiolle (NXS[2]-NES l‰htien) NXS[2] sis‰lt‰‰   */
                  /* energiapisteiden m‰‰r‰n energiagridiss‰, kun taas NES */
                  /* sis‰lt‰‰ reaktion energiapisteiden m‰‰r‰n. Viimeinen  */
                  /* piste j‰tet‰‰n aina levent‰m‰tt‰ (3. argumentissa -1) */

                  BroadCrossSection(L0, L1, UNRES - (NXS[2] - NES) - 1, T, Td,
                                    awr, XSS, &limits[(NXS[2] - NES)*3],
                                    NXS[2] - NES, NXS[2]);
                }
            }

          /*******************************************************************/

          /* Set XS temperature */

          WDB[nuc + NUCLIDE_XS_TEMP] = Td;
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Free memory and print newline */

  if (limits != NULL)
    Mem(MEM_FREE, limits);

  fprintf(outp, "\n");
}

/*****************************************************************************/

/*find_limits laskee limits-vektoriin integrointirajat ja huolehtii
  tarvittaessa muistin varauksesta. Muistia varataan oletusarvoisesti
  STD_LIMITS_SIZE verran, mutta funktio osaa hoitaa myˆs tilanteen, jossa
  vaikutusalan energiapisteiden m‰‰r‰ on suurempi kuin STD_LIMITS_SIZE.
  (esmes nopeissa kirjastoissa t‰ysin mahdollista, jopa todenn‰kˆist‰)

  Limits-vektori sis‰lt‰‰ rajat:
  i*3+0 = plus-termin integroinnin alaraja
  i*3+1 = - '' - yl‰raja
  i*3+2 = miinus-termin yl‰raja */

long *find_limits(long L0, long NES, double T, double Td, double awr,
                  double *XSS, long *limits, long UNRES){
  long i;
  double alpha, intlim, v;

  /* Varataan muistia, jos limits ei ole alustettu tai se on liian pieni */

  if ((limits == NULL) && (NES < STD_LIMITS_SIZE))
    limits = (long *)Mem(MEM_ALLOC, STD_LIMITS_SIZE*3, sizeof(long));
  else if ((NES > STD_LIMITS_SIZE) && (limits != NULL))
    limits = (long *)Mem(MEM_REALLOC, limits, NES*3*sizeof(long));
  else if ((NES > STD_LIMITS_SIZE) && (limits == NULL))
    limits = (long *)Mem(MEM_ALLOC, NES*3, sizeof(long));

  alpha = (awr*NEUTRON_E0*1E4)/(SPD_C*SPD_C*2*KELVIN*(Td-T));

  /* Integrointialueeksi valitaan n=4 "keskihajontaa", jonka havaittiin
   olevan riitt‰v‰n suuri lukema. (>=3.5 on about ok, NJOYlla n=4) */

  intlim=4/sqrt(alpha);

  limits[0]=0;

  /* Etsit‰‰n sopivat integrointirajat jokaiselle energiapisteelle */

  for(i=0; i<UNRES; i++){

    v=sqrt(2*XSS[L0+i]*SPD_C*SPD_C/(1E4*NEUTRON_E0));

    /* T‰m‰ s‰‰st‰‰ paljon iterointiaikaa, kun etsint‰‰ ei tarvitse
       aloittaa nollasta. Alarajat kasvavat monotonisesti. */
    if(i!=0)
      limits[i*3]=limits[(i-1)*3];

    /*e^-alpha(V-v)^2 -termin integrointi v-4/sqrt(alpha) :sta v+4/sqrt(alpha):n */

    /* NES -> NES-1 koska yl‰rajan indeksi ei saa olla yht‰ suuri kuin NES (TVi 1.9.2011) */
    while(limits[i*3]<NES-1 && sqrt(2*XSS[L0+limits[i*3]+1]*SPD_C*SPD_C/(1E4*NEUTRON_E0))<v-intlim){
      limits[i*3]++;
    }

    limits[i*3+1]=limits[i*3];

    /* NES -> NES-1 koska yl‰rajan indeksi ei saa olla yht‰ suuri kuin NES (TVi 1.9.2011) */
    while(limits[i*3+1]<NES-1 && sqrt(2*XSS[L0+limits[i*3+1]]*SPD_C*SPD_C/(1E4*NEUTRON_E0))<v+intlim){
      limits[i*3+1]++;
    }

    /*e^-alpha(V+v)^2 -termin integrointi nollasta 4/sqrt(alpha)-v :hen */

    limits[i*3+2]=0;

    /* NES -> NES-1 koska yl‰rajan indeksi ei saa olla yht‰ suuri kuin NES (TVi 1.9.2011) */
    while(limits[i*3+2]<NES-1 && sqrt(2*XSS[L0+limits[i*3+2]]*SPD_C*SPD_C/(1E4*NEUTRON_E0))<intlim-v){
        limits[i*3+2]++;
    }
  }
  return limits;
}

/*****************************************************************************/
