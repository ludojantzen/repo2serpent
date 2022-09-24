/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : element.c                                      */
/*                                                                           */
/* Created:       2016/02/12 (JLe)                                           */
/* Last modified: 2019/10/26 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Decomposes natural elements into isotopic compositions       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "natural_elements.h"
#include "element_data.h"

#define FUNCTION_NAME "Element:"

/*****************************************************************************/

void Element(long iso, char *name, char *dens, char *id)
{
  double d, tot;
  long Z, A, n, i;
  char tmpstr[MAX_STR];

  /* Get Z */

  if ((Z = atol(name)) < 1)
    for (n = 1; n < NUMBER_OF_ELEMENTS; n++)
      if (!strcasecmp(name, element_symbols[n]))
        {
          /* Set Z */

          Z = n;

          /* Break loop */

          break;
        }

  /* Check */

  if (Z < 1)
    Error(0, "Element %s does not exist", name);
  else if (Z > NUMBER_OF_ELEMENTS)
    {
      /* Try dividing by 1000 */

      Z = (long)((double)Z/1000.0);

      if ((Z < 1) || (Z > NUMBER_OF_ELEMENTS))
        Error(0, "Element %s does not exist", name);
    }

  /* Check if natural composition exists */

  n = 0;
  while(nat_frac[n][0] > 0)
    {
      /* Check Z */

      if (nat_frac[n][0] == Z)
        break;

      /* Next */

      n++;
    }

  /* Check */

  if (nat_frac[n][0] == Z)
    {
      /* Check if printed */

      if (iso < VALID_PTR)
        fprintf(outp, "\nIsotopic composition for natural %s:\n\n",
                element_names[Z]);
    }
  else
    {
      /* Not found */

      if (iso < VALID_PTR)
        fprintf(outp, "\nElement %s has no natural isotopes.\n\n",
                element_names[Z]);

      /* Exit */

      return;
    }

  /* Reset count */

  i = 0;

  /* Check density */

  if ((d = atof(dens)) < 0.0)
    {
      /* Calculate total for mass fraction */

      tot = 0.0;

      n = 0;
      while(nat_frac[n][0] > 0)
        {
          /* Check Z and add to total */

          if (nat_frac[n][0] == Z)
            tot = tot + nat_frac[n][2]*nat_frac[n][3];

          /* Next */

          n++;
        }

      /* Loop over data */

      n = 0;
      while(nat_frac[n][0] > 0)
        {
          /* Check Z and print */

          if (nat_frac[n][0] == Z)
            {
              /* Get A */

              A = (long)nat_frac[n][1];

              /* Ta-180 is actually isomeric */

              if ((Z == 73) && (A == 180))
                A = A + 200;

              /* Check if printed */

              if (iso < VALID_PTR)
                {
                  /* Check if id is given */

                  if (id == NULL)
                    fprintf(outp, "%6ld  %1.5E\n", 1000*Z + A,
                            nat_frac[n][2]*nat_frac[n][3]*d/tot);
                  else
                    fprintf(outp, "%6ld.%s  %1.5E\n", 1000*Z + A, id,
                            nat_frac[n][2]*nat_frac[n][3]*d/tot);
                }
              else
                {
                  /* Create a new entry */

                  if (i > 0)
                    iso = DuplicateItem(iso);
                  else
                    i++;

                  /* Put name */

                  sprintf(tmpstr, "%ld.%s", 1000*Z + A, id);
                  WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)PutText(tmpstr);

                  /* Put density */

                  WDB[iso + COMPOSITION_ADENS] =
                    nat_frac[n][2]*nat_frac[n][3]*d/tot;
                }
            }

          /* Next */

          n++;
        }

    }
  else if (d > 0.0)
    {
      /* Loop over data */

      n = 0;
      while(nat_frac[n][0] > 0)
        {
          /* Check Z and print */

          if (nat_frac[n][0] == Z)
            {
              /* Get A */

              A = (long)nat_frac[n][1];

              /* Ta-180 is actually isomeric */

              if ((Z == 73) && (A == 180))
                A = A + 200;

              /* Check if printed */

              if (iso < VALID_PTR)
                {
                  /* Check if id is given */

                  if (id == NULL)
                    fprintf(outp, "%6ld  %1.5E\n", 1000*Z + A,
                            nat_frac[n][3]*d);
                  else
                    fprintf(outp, "%6ld.%s  %1.5E\n", 1000*Z + A, id,
                            nat_frac[n][3]*d);
                }
              else
                {
                  /* Create a new entry */

                  if (i > 0)
                    iso = DuplicateItem(iso);
                  else
                    i++;

                  /* Put name */

                  sprintf(tmpstr, "%ld.%s", 1000*Z + A, id);
                  WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)PutText(tmpstr);

                  /* Put density */

                  WDB[iso + COMPOSITION_ADENS] = nat_frac[n][3]*d;
                }
            }

          /* Next */

          n++;
        }
    }
  else if (iso < VALID_PTR)
    Error(0, "Zero density");

  /* Exit OK */

  if (iso < VALID_PTR)
    fprintf(outp, "\n");
}

/*****************************************************************************/
