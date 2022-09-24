/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : potcorr.c                                      */
/*                                                                           */
/* Created:       2011/10/12 (TVi)                                           */
/* Last modified: 2019/12/18 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Temperature correction to elastic potential scattering xs    */
/*                                                                           */
/* Comments: -  Tämä funktio on suoraa lainaa wanhasta. Vanha ElaCorr        */
/*              ilmeisesti tekisi saman asian hieman optimoidummin, mutta    */
/*              pidetään nyt alkajaisiksi simppelinä.                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PotCorr:"

/*****************************************************************************/

double PotCorr(long nuc, double E, double kT)
{
  double a, ykspera;

  /* Tää sama testi tehdään dopmicroxs.c:ssä */

  if ((((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_USED) &&
       (E > RDB[nuc + NUCLIDE_URES_EMIN])) || (E > 1.0))
    return 1.0;

  /* Check TMS and DBRC flags */

  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS)
    kT = kT - RDB[nuc + NUCLIDE_TMS_MIN_TEMP]*KELVIN;
  else if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC))
    Die(FUNCTION_NAME, "Nuclide %s is not flagged for TMS or DBRC",
        GetText(nuc + NUCLIDE_PTR_NAME));

  /* Check zero change (temperature is set to TMS minimum if point  */
  /* is outside temperature distribution, the value can be slightly */
  /* negative due to limited numerical precision) */

  CheckValue(FUNCTION_NAME, "kT", "", kT, -1E-6, INFTY);

  if (kT <= 0.0)
    return 1.0;

  /* Calculate constant & return 1 if E > 500 kT/A,  22 ~= sqrt(500). */
  /* Approksimaation suht. virhe < 0.001 */

  if( ( a = sqrt(RDB[nuc + NUCLIDE_AWR]*E/kT ) ) > 250)
    return 1.0;

  ykspera=1/a;

  /* Check */

  CheckValue(FUNCTION_NAME, "a", "", a, ZERO, INFTY);

  /* Nain isoilla arvoilla voidaan jo approksimoida. Kun a > 2.568, */
  /* approksimaation suht virhe < 0.001 */

  if ( a > 2.568 )
    return 1.0 + 0.5*ykspera*ykspera;

  /* Jos tahan mennessa ei ole tarpannyt, lasketaan vaikeimman kautta. */
  /* 0.56419 = 1 / sqrt(PI) */

  return (1.0 + 0.5*ykspera*ykspera)*erf(a) + exp(-a*a)*ykspera*0.56419;
}

/*****************************************************************************/
