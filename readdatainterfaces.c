/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readdatainterfaces.c                           */
/*                                                                           */
/* Created:       2018/06/04 (VVa)                                           */
/* Last modified: 2018/06/05 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads general data-fields that can be updated                */
/*                                                                           */
/* Comments:  -                                                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadDataInterfaces:"

/*****************************************************************************/

void ReadDataInterfaces()
{
  long loc0, ptr, nuc, n;
  double val, minval, maxval;
  char fname[MAX_STR];
  FILE *fp;

  if ((loc0 = (long)RDB[DATA_PTR_DATAIFC0]) < VALID_PTR)
    return;

  fprintf(outp, "Reading general data-interfaces...\n\n");

  /* Loop over the interfaces to read data */

  while (loc0 > VALID_PTR)
    {
      /* Print data-file name into string */

      sprintf(fname, "%s", GetText(loc0 + DATAIFC_PTR_INPUT_FNAME));

      /* Simply read all values (should be numerical) */

      fp = fopen(fname, "r");

      /* Get number of values to read */

      n = (long)RDB[loc0 + DATAIFC_DATA_MEM_SIZE];

      /* Get pointer to data-array */

      ptr = (long)RDB[loc0 + DATAIFC_PTR_DATA];

      /* Reset min and max values */

      minval =  INFTY;
      maxval = -INFTY;

      /* Read in all values */

      for (; n > 0; n--)
        {
          if (fscanf(fp, "%lf", &val) == EOF)
            Die(FUNCTION_NAME, "Unexpected EOF when reading data from %s, expected %ld values",
                fname, (long)RDB[loc0 + DATAIFC_DATA_MEM_SIZE]);

          /* Keep track of minimum and maximum */

          if (val < minval)
            minval = val;

          if (val > maxval)
            maxval = val;

          /* Store value to array */

          WDB[ptr++] = val;
        }

      /* Store minimum and maximum value */

      WDB[loc0 + DATAIFC_MIN_VAL] = minval;
      WDB[loc0 + DATAIFC_MAX_VAL] = maxval;

      /* Conduct a sanity check for the values */

      /* CheckDataInterfaceValues(loc0); */

      /* Next interface */

      loc0 = NextItem(loc0);
    }

  /******************************************************************/
  /* Reset external nuclide max concentrations for extra dtmajorant */
  /******************************************************************/

  loc0 = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];

  while (loc0 > VALID_PTR)
    {
      /* Check type and reset max adens */

      if ((long)RDB[loc0 + MAJORANT_EXTRA_PTR_NUC] > VALID_PTR)
          if ((long)RDB[loc0 + MAJORANT_EXTRA_TYPE] == MAJORANT_EXTRA_NUCLIDE_EXT)
              WDB[loc0 + MAJORANT_EXTRA_FRAC] = 0.0;

      /* Next extra XS */

      loc0 = NextItem(loc0);
    }

  /*******************************************************************/
  /* Update external nuclide max concentrations for extra dtmajorant */
  /*******************************************************************/

  loc0 = (long)RDB[DATA_PTR_DATAIFC0];

  while (loc0 > VALID_PTR)
    {
      if ((long)RDB[loc0 + DATAIFC_DATA_TYPE] == DATA_TYPE_IFC_ADENS)
        {
          if ((nuc = (long)RDB[loc0 + DATAIFC_PTR_NUCLIDE]) > VALID_PTR)
            {
              /* Find nuclide in majorant extra list */

              ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
              while (ptr > VALID_PTR)
                {
                  /* Check nuclide pointer */

                  if ((long)RDB[ptr + MAJORANT_EXTRA_PTR_NUC] == nuc)
                    break;

                  /* Next */

                  ptr = NextItem(ptr);
                }

              /* Check if found */

              if (ptr < VALID_PTR)
                Die(FUNCTION_NAME, "Extra majorant block should have been created for "
                    " %s (data interface atomic density) previously but could not be found",
                    GetText(nuc + NUCLIDE_PTR_NAME));

              /* Get max atomic density for this data interface */

              maxval = RDB[loc0 + DATAIFC_MAX_VAL];

              /* Update majorant atomic density */
              fprintf(outp, "Maybe putting majorant density of %E for %s\n",
                      val, GetText(nuc + NUCLIDE_PTR_NAME));

              /* If this is larger than the current maximum, update the maximum */

              if (maxval > RDB[ptr + MAJORANT_EXTRA_FRAC])
                WDB[ptr + MAJORANT_EXTRA_FRAC] = maxval;
            }
        }

      loc0 = NextItem(loc0);
    }

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
