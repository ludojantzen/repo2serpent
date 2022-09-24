/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : radgammasrc.c                                  */
/*                                                                           */
/* Created:       2012/03/29 (JLe)                                           */
/* Last modified: 2020/06/29 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Gamma source from radioactive decay                          */
/*                                                                           */
/* Comments: - Käytetään nyt noita painoja (jako kokonaislähdetilavuudella)  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RadGammaSrc:"

/*****************************************************************************/

long RadGammaSrc(long src, long mat, double *E, double *wgt, long *idx,
                 long id)
{
  long mat0, mat1, nuc, loc0, ptr, ne, i, type, rad;
  double I, rnd, max, rate, tot, vol, norm;
  const double *x, *pdf, *cdf;

  /* Check source pointer */

  CheckPointer(FUNCTION_NAME, "(src)", DATA_ARRAY, src);

  /* Check if point is in void */

  if (mat < VALID_PTR)
    return -1;

  /* Get particle type */

  type = (long)RDB[src + SRC_TYPE];
  CheckValue(FUNCTION_NAME, "type", "", type, 1, 2);

  /***************************************************************************/

  /***** Step 1: Rejection sampling in material ******************************/

  /* Check if source material is defined */

  if ((mat0 = (long)RDB[src + SRC_PTR_RAD_SRC_MAT]) > VALID_PTR)
    {
      /* Check pointer */

      if (mat != mat0)
        {
          /* Check if material has parent */

          if ((mat1 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
            return -1;

          /* Check parent */

          if (mat1 != mat0)
            return -1;
        }
    }

  /* Avoid compiler warning */

  rate = -1.0;
  max = -1.0;
  tot = -1.0;
  vol = -1.0;

  /* Get source parameters */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      /* Copy data to variables */

      rate = RDB[mat + MATERIAL_NEUTRON_DEC_SRC_RATE];
      max = RDB[DATA_NEUTRON_DEC_SRC_MAX_I];
      tot = RDB[DATA_TOT_NEUTRON_DEC_SRC_RATE];
      vol = RDB[DATA_NEUTRON_DEC_SRC_VOL];

      /* Check total rate */

      if (tot == 0.0)
        Error(src, "Decay source but no neutron emission");
    }
  else if (type == PARTICLE_TYPE_GAMMA)
    {
      /* Copy data to variables */

      rate = RDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE];
      max = RDB[DATA_PHOTON_DEC_SRC_MAX_I];
      tot = RDB[DATA_TOT_PHOTON_DEC_SRC_RATE];
      vol = RDB[DATA_PHOTON_DEC_SRC_VOL];

      /* Check total rate */

      if (tot == 0.0)
        Error(src, "Decay source but no photon emission");
    }
  else
    Die(FUNCTION_NAME, "Invalid source type");

  /* Check source rate */

  if (rate == 0.0)
    return -1;

  /* Set normalization factor (tämä on kerroin jolla hiukkasen paino    */
  /* kerrotaan kun emittoituneita säteilylajeja on useita ja normeeraus */
  /* on kiinnitetty yhteen niistä). */

  if (mat0 > VALID_PTR)
    norm = rate/RDB[DATA_NORM_DECAY_SRC_RATE];
  else
    norm = tot/RDB[DATA_NORM_DECAY_SRC_RATE];

  /* Check */

  if (norm < ZERO)
    Die(FUNCTION_NAME, "Zero weight normalization factor");

  /* Check other parameters */

  if (RDB[mat + MATERIAL_VOLUME] == 0.0)
    Die(FUNCTION_NAME, "Material volume is zero");
  else if (max == 0.0)
    Die(FUNCTION_NAME, "Maximum intensity is zero");

  /* Check mode */

  if ((long)RDB[src + SRC_RAD_SRC_MODE] == RAD_SRC_MODE_ANALOG)
    {
      /* Perform rejection sampling on material */

      if (RandF(id) > rate/RDB[mat + MATERIAL_VOLUME]/max)
        return -1;
    }
  else if ((long)RDB[src + SRC_RAD_SRC_MODE] == RAD_SRC_MODE_IMPLICIT)
    {
      /* Or adjust weight */

      *wgt = *wgt*rate/tot/RDB[mat + MATERIAL_VOLUME]*vol;
    }
  else
    Die(FUNCTION_NAME, "Sampling mode is not set");

  /* Adjust weight */

  *wgt = *wgt*norm;
  CheckValue(FUNCTION_NAME, "(wgt)", "", *wgt, ZERO, INFTY);

  /***************************************************************************/

  /***** Step 2: Sample nuclide **********************************************/

  /* This is used as cut-off in processdecaysrc.c. (could just check the */
  /* pointer, but this is better for debugging) */

  if (rate/tot < 1E-19)
    return -1;

  /* Pointer to list */

  if (type == PARTICLE_TYPE_NEUTRON)
    ptr = (long)RDB[mat + MATERIAL_PTR_NEUTRON_DECAY_SRC];
  else
    ptr = (long)RDB[mat + MATERIAL_PTR_PHOTON_DECAY_SRC];

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Sample random number */

  rnd = RandF(id);

  /* Loop over list */

  while (ptr > VALID_PTR)
    {
      if (RDB[ptr + SRC_DECCAY_CUM_P] > rnd)
        break;

      /* Next */

      ptr = NextItem(ptr);
    }

 /* Check pointer */

  if (ptr < VALID_PTR)
    Die(FUNCTION_NAME, "Unable to sample nuclide %s",
        GetText(mat + MATERIAL_PTR_NAME));

  /***************************************************************************/

  /***** Step 3: Sample emission energy **************************************/

  /* Put index */

  *idx = (long)RDB[ptr + SRC_DECCAY_IDX];

  /* Get pointer to nuclide */

  nuc = (long)RDB[ptr + SRC_DECCAY_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Find spectrum */

  rad = (long)RDB[nuc + NUCLIDE_PTR_RADIATIONS];
  while (rad > VALID_PTR)
    {
      /* Check type */

      if ((long)RDB[rad + NUCLIDE_RAD_TYPE] == type)
        break;

      /* Pointer to next */

      rad = NextItem(rad);
    }

  /* Check pointer */

  if (rad < VALID_PTR)
    Die(FUNCTION_NAME, "Pointer error");

  /* Get total intensity */

  I = RDB[rad + NUCLIDE_RAD_SPEC_I];
  CheckValue(FUNCTION_NAME, "I", "", I, ZERO, INFTY);

  /* Sample fraction of intensity */

  I = I*RandF(id);

  /* Get pointer to spectra */

  loc0 = (long)RDB[rad + NUCLIDE_RAD_PTR_SPEC];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Loop over spectra */

  while (loc0 > VALID_PTR)
    {
      /* Compare to intensity */

      if ((I = I - RDB[loc0 + DECAY_SPEC_RI]) < 0.0)
        break;

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Check pointer */

  if (loc0 < VALID_PTR)
    Die(FUNCTION_NAME, "Emission sampling failed for %s",
        GetText(nuc + NUCLIDE_PTR_NAME));

  /* Check type */

  if ((long)RDB[loc0 + DECAY_SPEC_TYPE] == DECAY_SPEC_LINE)
    {
      /* Discrete line, get energy */

      *E = RDB[loc0 + DECAY_SPEC_LINE_E];
    }
  else if ((long)RDB[loc0 + DECAY_SPEC_TYPE] == DECAY_SPEC_CONT)
    {
      /* Continuous distribution, get pointer to energy array */

      ptr = (long)RDB[loc0 + DECAY_SPEC_CONT_PTR_E];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      x = &RDB[ptr];

      /* Get pointer to PDF array */

      ptr = (long)RDB[loc0 + DECAY_SPEC_CONT_PTR_PDF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      pdf = &RDB[ptr];

      /* Get pointer to CDF array */

      ptr = (long)RDB[loc0 + DECAY_SPEC_CONT_PTR_CDF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      cdf = &RDB[ptr];

      /* Get number of points and interpolation scheme */

      ne = (long)RDB[loc0 + DECAY_SPEC_CONT_NE];
      CheckValue(FUNCTION_NAME, "ne", "", ne, 2, 10000);

      i = (long)RDB[loc0 + DECAY_SPEC_CONT_INTT];
      CheckValue(FUNCTION_NAME, "ne", "", i, 1, 5);

      /* Sample energy */

      if ((*E = SampleTabular(x, pdf, cdf, ne, i, id)) < 0.0)
        return -1;
    }
  else
    Die(FUNCTION_NAME, "Invalid emission type");

  /* Check sampled energy */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);

  /* Return material pointer */

  return mat;

  /***************************************************************************/
}

/*****************************************************************************/
