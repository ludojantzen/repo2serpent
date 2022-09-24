/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : otfsabxs.c                                     */
/*                                                                           */
/* Created:       2015/03/20 (TVi)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Returns interpolated cross sections for on-the-fly S(a,b)    */
/*              treatment.                                                   */
/*                                                                           */
/* Comments: - JLe 1.11.2015 / 2.1.25 Poistin tuolta tallennuksen            */
/*             REACTION_PTR_PREV_XS:iin (liittyy ures-sämpläykseen)          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "OTFSabXS:"

/*****************************************************************************/

double OTFSabXS(long rea, double E, double T, long id){

  long nuc, sab1, sab2, nuc1, nuc2, mt, rea2, sab0, ptr, ncol; 
  double f, xs1, xs2, T1, T2, xs;
 
  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Avoid compiler warning */

  sab1 = -1;

  /* Get reaction mt and rea nuclide */

  mt = (long)RDB[rea + REACTION_MT];

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  if (E > RDB[nuc + NUCLIDE_SAB_EMAX])
    return MicroXS(rea, E, id);

  /* Check that data exists */

  if ((sab2 = (long)RDB[nuc + NUCLIDE_PTR_SAB]) < VALID_PTR )
    Die(FUNCTION_NAME, "S(a,b) data not available for nuclide %s", 
        GetText(nuc + NUCLIDE_PTR_NAME));
  
  /* Find correct temperature */
  
  while(sab2 > VALID_PTR)
    {      

      if (RDB[sab2 + SAB_T] > T ||
          (RDB[sab2 + SAB_T] == T && 
           NextItem(sab2) < VALID_PTR ) )
        break;

      sab2 = NextItem(sab2);
    }

  /* Check that sab was found */

  if (sab2 < VALID_PTR)
    Die(FUNCTION_NAME, "S(a,b) OTF nuclide not found for %s", 
        GetText(nuc + NUCLIDE_PTR_NAME));

  sab1 = PrevItem(sab2);

  /* Check Pointers */

  CheckPointer(FUNCTION_NAME, "(sab1)", DATA_ARRAY, sab1);
  CheckPointer(FUNCTION_NAME, "(sab2)", DATA_ARRAY, sab2);

  /* Temperatures */

  T1 = RDB[sab1 + SAB_T];
  T2 = RDB[sab2 + SAB_T];

  /* Pointers to S(a,b) nuclides */

  nuc1 = (long)RDB[sab1 + SAB_PTR_ISO];
  nuc2 = (long)RDB[sab2 + SAB_PTR_ISO];

  CheckPointer(FUNCTION_NAME, "(nuc1)", DATA_ARRAY, nuc1);
  CheckPointer(FUNCTION_NAME, "(nuc2)", DATA_ARRAY, nuc2);

  /* Find reaction pointers in S(a,b) data */
  
  /* Elastic scattering xs = total xs in case of S(a,b) nuclides 
     (which only have 1004 and 1002 reactions and total is 
     calculated over them) */
  
  if ((mt == 1) || (mt == 2))
    rea2 = (long)RDB[nuc1 + NUCLIDE_PTR_TOTXS];
  
  /* Find reaction at first temperature */

  else
    {
      /* Pointer to reaction */

      rea2 = (long)RDB[nuc1 + NUCLIDE_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea2)", DATA_ARRAY, rea2);

      /* JLe: Tässä voisi käyttää SeekList():iä */

      while (rea2 > VALID_PTR)
        {     
          if ((long)RDB[rea2 + REACTION_MT] == mt - 1000)
            break;
      
          rea2 = NextItem(rea2);
        }
    
      if (rea2 < VALID_PTR)
        Die(FUNCTION_NAME, "mt %ld not found for S(a,b) nuclide %s", mt, 
            GetText(nuc1 + NUCLIDE_PTR_NAME));
    }

  /* Get cross section */

  xs1 = MicroXS(rea2, E, id);
  CheckValue(FUNCTION_NAME, "xs1", "", xs1, 0.0, 1E+10);

  if ((mt == 1) || (mt == 2))
    rea2 = (long)RDB[nuc2 + NUCLIDE_PTR_TOTXS];

    /* Find reaction at second temperature */

  else
    {
      rea2 = (long)RDB[nuc2 + NUCLIDE_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea2)", DATA_ARRAY, rea2);

      /* JLe: Tässä voisi käyttää SeekList():iä */

      while (rea2 > VALID_PTR)
        {
          if ((long)RDB[rea2 + REACTION_MT] == mt - 1000)
            break;
      
          rea2 = NextItem(rea2);
        }
    
      if (rea2 < VALID_PTR)
        Die(FUNCTION_NAME, "mt %ld not found for S(a,b) nuclide %s", mt, 
            GetText(nuc1 + NUCLIDE_PTR_NAME));
    }

  /* Get cross section */

  xs2 = MicroXS(rea2, E, id);
  CheckValue(FUNCTION_NAME, "xs2", "", xs2, 0.0, 1E+10);

  /* Avoid compiler warning */
  
  f = -1.0;

  /* Calculate factor */

  if (T1 != T2)
    f = (T-T1)/(T2-T1);
  else
    Die(FUNCTION_NAME, "Division by zero");

  /* Check value */

  CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);

  /* Get collision number */
  
  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Pointer to first item in sab list */

  sab0 = (long)RDB[nuc + NUCLIDE_PTR_SAB];

  /* Store values */

  StoreValuePair(sab0 + SAB_PTR_PREV_FRAC, (double)ncol, f, id);
  StoreValuePair(sab0 + SAB_PTR_PREV_SAB1, (double)ncol, (double)sab1, id);

  xs = xs1 + f*(xs2-xs1);

  /* Tämä poistettiin juuri ennen 2.1.25:n jakelua, sillä se sotkee */
  /* myöhemmin käytettävän kokonaisvaikutusalan arvon (JLe / 19.2.2016). */
  
  /*
  StoreValuePair(rea + REACTION_PTR_PREV_XS, E, xs, id);
  */

  /* Näillä voisi varmaan optimoida myös tuota rutiinin alkupäätä */

  /* JLe: Jos palautettava arvo on samalla energialla aina sama, niin   */
  /* alkuun voi laittaa TestValuePair() -kutsun. Onko nuo kaksi ylempää */
  /* tosiaan tarkoitus kiinnittää törmäykseen eikä energiaan? */

  return xs;
}

/*****************************************************************************/
