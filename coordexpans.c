/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : coordexpans.c                                  */
/*                                                                           */
/* Created:       2013/04/04 (JLe)                                           */
/* Last modified: 2016/09/16 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Coordinate transformations for thermal expansion             */
/*                                                                           */
/* Comments: - Hot-to-cold -muunnosta tarvitaan laskettaessa tehojakaumaa    */
/*             ja materiaalikohtaisia reaktionopeuksia esim. palamalaskussa. */
/*             Kutsutaan nyt vielä pelkästään scoreinterfacepower.c:stä,     */
/*             ja noi muut jutut pitää vielä miettiä tarkemmin (esim. että   */
/*             tallennetaanko muunnetut koordinaatit suoraan rakenteisiin    */
/*             joista muut aliohjelmat ne lukee, tehdäänkö muunnos aina      */
/*             paikallisesti. Esim. fissioneutronit pitäisi tallentaa        */
/*             muuntamattomilla koordinaateilla, muuten lähdepiste siirtyy.  */
/*             Voi kyllä olla että koko homma menee kerralla oikein jos      */
/*             muunnoksen tekee suoraan whereami.c:stä.)                     */
/*                                                                           */
/*           - Cold-to-hot -muunnosta tarvitaan korjaamaan neutronin paikkaa */
/*             aikariippuvassa simulaatiossa kun geometria muuttuu askelten  */
/*             välissä (pitää vielä miettiä että kannattaako sitä korjata    */
/*             ollenkaan).                                                   */
/*                                                                           */
/*           - 2.1.19 added time dependence (reworked in 2.1.27)             */
/*                                                                           */
/*           - TODO: Axial and angular limits to arrays to use binary search */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CoordExpans:"

/*****************************************************************************/

void CoordExpans(long loc0, double *x, double *y, double *z, double t,
                 long mode)
{
  long loc1, ang, nr, ptr, i;
  double rc0 = 0.0, rc1, rh0 = 0.0, rh1, rh, rc, cost, sint, phi, phi2;

  /* Find axial zone */

  loc1 = (long)RDB[loc0 + IFC_FUEP_PTR_AX];

  while (loc1 > VALID_PTR)
    {
      /* Compare coordinates */

      if ((*z >= RDB[loc1 + IFC_FUEP_AX_ZMIN]) &&
          (*z < RDB[loc1 + IFC_FUEP_AX_ZMAX]))
        break;

      /* Next */

      loc1 = NextItem(loc1);
    }

  /* Check if found */

  if (loc1 < VALID_PTR)
    return;

  /* Find angular zone */

  ang = (long)RDB[loc1 + IFC_FUEP_AX_PTR_ANG];

  /* Get polar angle */

  phi = PolarAngle(*x,*y);

  while (ang > VALID_PTR)
    {
      /* Rotate if needed */
      if(phi > 2.0*PI+RDB[ang + IFC_FUEP_ANG_AMIN])
        phi2 = phi - 2.0*PI;
      else
        phi2 = phi;

      /* Compare coordinates */

      if ((phi2 >= RDB[ang + IFC_FUEP_ANG_AMIN]) &&
          (phi2 < RDB[ang + IFC_FUEP_ANG_AMAX]))
        break;

      /* Next */

      ang = NextItem(ang);
    }

  /* Check if found */

  if (ang < VALID_PTR)
    return;

  /* Get number of radial nodes */

  nr = (long)RDB[ang + IFC_FUEP_ANG_N_RAD];

  /* Check mode */

  if (mode == 1)
    {
      /***********************************************************************/

      /***** Hot-to-cold *****************************************************/

      /* First find the radial zone where the point is, use BOI geometry */

      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2_BOI];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      i = SearchArray(&RDB[ptr], *x*(*x)+*y*(*y), nr);

      /* Check if found */
      if (i < 0)
        return;

      /* Get the outer and inner radii of the zone */

      /* Hot radii */
      rh1 = sqrt(RDB[ptr + i + 1]);
      rh0 = sqrt(RDB[ptr + i]);

      /* Cold radii */
      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      rc1 = sqrt(RDB[ptr + i + 1]);
      rc0 = sqrt(RDB[ptr + i]);

      /* Calculate the hot radius of the point as well as directional */
      /* cosine and sine */

      rh = sqrt(*x*(*x) + *y*(*y));
      cost = *x/rh;
      sint = *y/rh;

      /* Calculate the cold radius via linear interpolation */
      rc = rc0 + (rc1-rc0)*(rh-rh0)/(rh1-rh0);

      /* Calculate new x and y values */

      *x = rc*cost;
      *y = rc*sint;

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Cold-to-hot *****************************************************/

      /* First find the radial zone where the point is */
      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      i = SearchArray(&RDB[ptr], *x*(*x)+*y*(*y), nr);

      /* Check if found */
      if (i < 0)
        return;

      /* Get the outer and inner radii of the zone */
      /* Cold radii */

      rc1 = sqrt(RDB[ptr + i + 1]);
      rc0 = sqrt(RDB[ptr + i]);

      /* Hot radii, use BOI geometry */

      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2_BOI];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      rh1 = sqrt(RDB[ptr + i + 1]);
      rh0 = sqrt(RDB[ptr + i]);

      /* Calculate the hot radius of the point as well as directional */
      /* cosine and sine */

      rc = sqrt(*x*(*x) + *y*(*y));
      cost = *x/rc;
      sint = *y/rc;

      /* Calculate the hot radius via linear interpolation */

      rh = rh0 + (rh1-rh0)*(rc-rc0)/(rc1-rc0);

      /* Calculate new x and y values */

      *x = rh*cost;
      *y = rh*sint;

      /************************************************************************/
    }
}

/*****************************************************************************/
