/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : uresfactor.c                                   */
/*                                                                           */
/* Created:       2011/01/08 (JLe)                                           */
/* Last modified: 2018/10/01 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns the ratio of ures-corrected to smooth cross section  */
/*                                                                           */
/* Comments: - Probability table sampling was moved to sampleptable.c        */
/*             (JLe 7.11.2015 / 2.1.25)                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UresFactor:"

/*****************************************************************************/

double UresFactor(long rea, double E, long id)
{
  long nuc, urs, ptr, mt, edep;
  double f, g, rnd, xs0, xs;

  /* Check Pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
  
  /* Check boundaries (yläraja oltava mukana) */

  if ((E < RDB[rea + REACTION_URES_EMIN]) ||
      (E >= RDB[rea + REACTION_URES_EMAX]))
    return 1.0;

  /* Get pointer to probability table data */

  urs = (long)RDB[rea + REACTION_PTR_URES];
  CheckPointer(FUNCTION_NAME, "(urs)", DATA_ARRAY, urs);

  /* Check existing data */
 
  if ((f = TestValuePair(urs + URES_PTR_PREV_FACT, E, id)) > -INFTY)
    {
      /* Use value only if it was obtained using same random number */

      if ((rnd = TestValuePair(urs + URES_PTR_RND_CHK, E, id)) > 0.0)
        if (rnd == TestValuePair(urs + URES_PTR_RND, E, id))
          return f;
    }

  /* Get reaction mt */

  mt = (long)RDB[rea + REACTION_MT];

  /* Get edep mode */

  edep = (long)RDB[DATA_EDEP_MODE];

  /* Check special */

  if ((mt == 101) || (mt == 1))
    {
      /***********************************************************************/

      /***** Microscopic total or total absorption xs ************************/

      /* Pointer to nuclide */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Check ures sampling flag (data may exist but sampling */
      /* is not in use because of dilution cut-off) */

      if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == NO)
        return 1.0;

      /* Get unadjusted total cross section */
      
      xs = UresDiluMicroXS(rea, E, id);
      CheckValue(FUNCTION_NAME, "xs", "", xs, 0.0, MAX_XS);

      /* Remember unadjusted value */

      xs0 = xs;

      /* Check mt */

      if (mt == 101)
        {
          /* Total absorption, adjust for capture */

          rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          xs = xs + MicroXS(rea, E, id);
          xs = xs - UresDiluMicroXS(rea, E, id);
        }
      else
        {
          /* Total, adjust for elastic */

          rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          xs = xs + MicroXS(rea, E, id);
          xs = xs - UresDiluMicroXS(rea, E, id);
          
          /* Adjust for capture */
          
          rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          xs = xs + MicroXS(rea, E, id);
          xs = xs - UresDiluMicroXS(rea, E, id);
          
          /* Adjust for fission */
          
          if ((rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS]) > VALID_PTR)
            {
              xs = xs + MicroXS(rea, E, id);
              xs = xs - UresDiluMicroXS(rea, E, id);
            }
        }
      
      /* Calculate factor */

      if (xs0 > 0.0)
        f = xs/xs0;
      else
        f = 1.0;
            
      /***********************************************************************/
    }
  else if (mt == 301)
    {
      /***********************************************************************/

      /***** Heating cross section *******************************************/

      /* Pointer to nuclide */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Check ures sampling flag (data may exist but sampling */
      /* is not in use because of dilution cut-off) */
      
      if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == NO)
        return 1.0;
      /*
      return 1.0;

      Die(FUNCTION_NAME, "Tää ei nyt pelaa, noi arvot on ihan kummallisia");
      */
      /* Sample probability table */

      f = SamplePTable(urs, E, mt, id);

      /* Pointer to total cross section */

      ptr = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Pointer to probability table data */

      ptr = (long)RDB[ptr + REACTION_PTR_URES];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Sample probability table */

      g = SamplePTable(ptr, E, mt, id);

      /* Check mode */
      
      if ((long)RDB[urs + URES_IFF] == 0)
        {          
          /* Calculate ures-corrected value */

          xs = f*g;
          
          if (edep < EDEP_MODE_NEUTRON_PHOTON)
            {
              CheckValue(FUNCTION_NAME, "xs", "", xs , -20.0, MAX_XS);
            }
          
          /* Get un-adjusted cross section */
          
          xs0 = UresDiluMicroXS(rea, E, id);
          
          if (edep < EDEP_MODE_NEUTRON_PHOTON)
            {          
              CheckValue(FUNCTION_NAME, "xs0", "", xs0, 0.0, MAX_XS);
            }
          
          /* Calculate factor */

          if (edep < EDEP_MODE_NEUTRON_PHOTON)
            {
              if (xs0 > 0.0)
                f = xs/xs0;
              else
                f = 1.0;
            }
          else
            {
              if (fabs(xs0) > ZERO)
                f = xs/xs0;
              else
                f = 1.0;
            }
        }
      else
        {
          /* Multiply factors */

          f = f*g;
        }

      /***********************************************************************/
    }
  else if (mt > 0)
    {
      /***********************************************************************/

      /***** Microscopic reaction cross section ******************************/

      /* Pointer to nuclide */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Check ures sampling flag (data may exist but sampling */
      /* is not in use because of dilution cut-off) */
      
      if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == NO)
        return 1.0;

      /* Sample probability table */

      f = SamplePTable(urs, E, mt, id);

      /* Check mode */
      
      if ((long)RDB[urs + URES_IFF] == 0)
        {
          /* Get un-adjusted cross section */
          
          xs0 = UresDiluMicroXS(rea, E, id);
          CheckValue(FUNCTION_NAME, "xs0", "", xs0, 0.0, MAX_XS);
          
          /* Calculate factor */
          
          if (xs0 > 0.0)
            f = f/xs0;
          else
            f = 1.0;
        }
      
      /***********************************************************************/
    }
  else 
    Die(FUNCTION_NAME, "Invalid reaction mode");

  /* Check factor */

  if ((edep != EDEP_MODE_NEUTRON_PHOTON) || (mt != 301))
    {
      CheckValue(FUNCTION_NAME, "f", "", f, -10.0, INFTY);
    }
  /* Remember value */

  StoreValuePair(urs + URES_PTR_PREV_FACT, E, f, id);

  /* Return factor */

  return f;
}

/*****************************************************************************/
