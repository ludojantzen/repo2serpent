/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sensbufval.c                                   */
/*                                                                           */
/* Created:       2018/06/21 (VVa)                                           */
/* Last modified: 2018/06/28 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates buffer value for a specific sensitivity bin index */
/*              Handles summing up sum reaction modes, sum bins and total    */
/*              bin in a nice manner                                         */
/*                                                                           */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SensBufVal:"
#define MY_MAXMT 100

/*****************************************************************************/

double SensBufVal(long stp, long gen, long imat,
                  long izai, long irea, long iene)
{
  long mtarr, maxmt, ireaNew, idx, sens, i, mat0, zai0;
  long reaarr, matarr, zaiarr, mt, mtlist[MY_MAXMT];
  double val;

  /* Reset value */

  val = 0.0;

  /* Get pointer to the sensitivity block or return 0.0 */

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return val;

  /* Get pointer to index array */

  reaarr = (long)RDB[sens + SENS_PTR_PERT_INDICES];
  CheckPointer(FUNCTION_NAME, "(reaarr)", DATA_ARRAY, reaarr);

  /* Get pointer to index array */

  matarr = (long)RDB[sens + SENS_PTR_MAT_INDICES];
  CheckPointer(FUNCTION_NAME, "(matarr)", DATA_ARRAY, matarr);

  /* Get pointer to index array */

  zaiarr = (long)RDB[sens + SENS_PTR_ZAI_INDICES];
  CheckPointer(FUNCTION_NAME, "(zaiarr)", DATA_ARRAY, zaiarr);

  /******************************/
  /* Handle ZAI integrated data */
  /******************************/

  /* Total ZAI */

  if (((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ZAI_TOT)
      && (izai == (long)RDB[zaiarr + SENS_TOT_ZAI_IDX]))
    {
      /* Loop over all ZAI bins */

      val = 0.0;

      /* Get index to first zai that should be dealt with */

      zai0 = (long)RDB[zaiarr + SENS_NON_ZAI_IDX];

      /* Do not call recursively for total or sum */

      for (i = zai0; i < (long)RDB[sens + SENS_N_ZAI]; i++)
        val = val + SensBufVal(stp, gen, imat, i, irea, iene);

      /* Return the sum */

      return val;
    }

  /* Sum ZAI */

  if (((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ZAI_SUM)
      && (izai == (long)RDB[zaiarr + SENS_SUM_ZAI_IDX]))
    {
      /* Loop over all bins starting from the next one from the sum bin */

      val = 0.0;

      /* Get index to first zai that should be dealt with */

      zai0 = (long)RDB[zaiarr + SENS_NON_ZAI_IDX] + 1;

      for (i = zai0; i < (long)RDB[sens + SENS_N_ZAI]; i++)
        val = val + SensBufVal(stp, gen, imat, i, irea, iene);

      return val;
    }

  /***********************************/
  /* Handle material integrated data */
  /***********************************/

  /* Total material */

  if (((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_TOT)
      && (imat == (long)RDB[matarr + SENS_TOT_MAT_IDX]))
    {
      /* Loop over all MAT bins */

      val = 0.0;

      /* Get index to first material that should be dealt with */

      mat0 = (long)RDB[matarr + SENS_NON_MAT_IDX];

      /* Do not call recursively for total or sum */

      for (i = mat0; i < (long)RDB[sens + SENS_N_MAT]; i++)
        val = val + SensBufVal(stp, gen, i, izai, irea, iene);

      /* Return the sum */

      return val;
    }

  /* Sum of included materials */

  if (((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_SUM)
      && (imat == (long)RDB[matarr + SENS_SUM_MAT_IDX]))
    {
      /* Loop over all bins starting from the next one from the sum bin */

      val = 0.0;

      /* Get index to first material that should be dealt with */

      mat0 = (long)RDB[matarr + SENS_NON_MAT_IDX] + 1;

      for (i = mat0; i < (long)RDB[sens + SENS_N_MAT]; i++)
        val = val + SensBufVal(stp, gen, i, izai, irea, iene);

      return val;
    }

  /* The remaining data cannot be sum or total in mat or zai */

  /*********************************/
  /* Handle energy integrated data */
  /*********************************/

  if (((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ENE_TOT)
      && (iene == 0))
    {
      /* Loop over all energy bins */

      val = 0.0;

      /* Include also bin zero (not found) as it contains scores that didn't */
      /* fall into any energy bin */

      idx = CompressSensLabel(imat, izai, irea, iene, 1.0);
      if (gen >= 0)
        val = BufVal(stp, gen, idx);
      else
        val = BufVal(stp, idx);

      /* Do not call bin zero recursively */

      for (i = 1; i < (long)RDB[sens + SENS_N_ENE] + 1; i++)
        val = val + SensBufVal(stp, gen, imat, izai, irea, i);

      /* Return the sum */

      return val;
    }

  /*************************/
  /* Handle sum reactions  */
  /*************************/

  /* Handle total nubar */

  if (((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_NUBAR)
      && (irea == (long)RDB[reaarr + NUBAR_TOT_IDX]))
    {
      /* Calculate as a sum from the delayed and prompt nubar */
      val = 0.0;
      val += SensBufVal(stp, gen, imat, izai, (long)RDB[reaarr + NUBAR_PRO_IDX], iene);
      val += SensBufVal(stp, gen, imat, izai, (long)RDB[reaarr + NUBAR_DEL_IDX], iene);

      return val;
    }

  /* Handle total fission spectrum */

  if (((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_CHI)
      && (irea == (long)RDB[reaarr + CHI_TOT_IDX]))
    {
      /* Calculate as a sum from the delayed and prompt spectra */
      val = 0.0;
      val += SensBufVal(stp, gen, imat, izai, (long)RDB[reaarr + CHI_PRO_IDX], iene);
      val += SensBufVal(stp, gen, imat, izai, (long)RDB[reaarr + CHI_DEL_IDX], iene);

      return val;
    }

  /* Handle total XS */

  if (irea == (long)RDB[reaarr + TOT_REA_IDX])
    {
      /* Reset sum */

      val = 0.0;

      /* Sum up all partial reactions */

      if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XS)
        {
          /* From different sum reaction modes ("ela", "inl", etc.) */

          for (i = 1; i < (long)RDB[sens + SENS_MAX_XS_INDEX]; i++)
            {
              /* Check if the sum reaction mode is included in calculation */

              if ((ireaNew = (long)RDB[reaarr + i]) > 0)
                {
                  /* Compress the label and add the value from the buffer */
                  /* Values for sum reactions may be 0 */

                  idx = CompressSensLabel(imat, izai, ireaNew, iene, 1.0);
                  if (gen >= 0)
                    val = val + BufVal(stp, gen, idx);
                  else
                    val = val + BufVal(stp, idx);
                }
            }

          /* Return sum of sum reaction modes */

          return val;
        }
      if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XSMT)
        {
          /* From different partial reaction modes (MT 2, etc.) */
          /* Loop over all possible partial MTs */

          /* Get maximum MT to include */

          maxmt = (long)RDB[sens + SENS_MAX_MT];

          /* Limit maxmt from up to 117 (n,d alpha) */
          /* Usually, higher MTs are not used */
          /* NOTE: MTs 600-891 can be used to give partial reaction modes */
          /*       for many production reactions */

          if (maxmt > 117)
            maxmt = 117;

          /* Get pointer to mt index array */

          mtarr = (long)RDB[sens + SENS_PTR_MT_INDICES];

          /* Loop over the different MTs and include also 0 to which */
          /* not included MTs have been scored during tracking */

          for (i = 0; i < maxmt; i++)
            {
              /* Skip mt if it is not included for calculation */

              if ((ireaNew = (long)RDB[mtarr + i]) == 0)
                continue;

              /* Compress the label and add the value from the buffer */
              /* Values for sum reactions may be 0 */

              idx = CompressSensLabel(imat, izai, ireaNew, iene, 1.0);
              if (gen >= 0)
                val = val + BufVal(stp, gen, idx);
              else
                val = val + BufVal(stp, idx);
            }

          /* Return sum of sum reaction modes */

          return val;
        }
    }

  /* Handle MT numbers that indicate a sum reaction (may be required to sum from partials) */

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XSMT)
    {

      /* Figure out the MT corresponding to this reaction */

      /* Get maximum perturbed MT */

      maxmt = (long)RDB[sens + SENS_MAX_MT];

      /* Limit maxmt to 18 (n,f), which is currently maximum sum reaction mode */
      /* for which data is summed up here */

      if (maxmt > 18)
        maxmt = 18;

      /* Get pointer to indices */

      mtarr = (long)RDB[sens + SENS_PTR_MT_INDICES];
      CheckPointer(FUNCTION_NAME, "(mtarr)", DATA_ARRAY, mtarr);

      /* Loop over all possible MT numbers */
      /* and break when found */

      for (mt = 0; mt < maxmt + 1; mt++)
        if ((long)RDB[mtarr + mt] == irea)
          break;

      /* Set the mtlist to zero */

      memset(mtlist, 0.0, MY_MAXMT*sizeof(long));

      /* Check if we found the mt */

      if (mt <= maxmt)
        {
          /* Create partial reaction list based on mt */
          /* Include also the sum reaction mode. */
          /* NB: During neutron transport, either the sum reaction mode is sampled */
          /*     or the partials so we should be able to sum them up (one is zero) */

          if (mt == 4)
            {
              /* Inelastic scattering */

              /* Include sum reaction mode if it is present for some nuclides */

              mtlist[0] = 4;

              /* Include partials if they are present for some nuclides */
              /* MTs 51-91 */

              for (i = 51; i < 92; i++)
                mtlist[1+i-51] = i;
            }
          else if (mt == 18)
            {
              /* Fission */

              /* Include total fission */

              mtlist[0] = 18;

              /* First, second and third chance fission */

              mtlist[1] = 19;
              mtlist[2] = 20;
              mtlist[3] = 21;
            }
        }

      /* If mtlist contains some MTs to sum up, let's sum them up and return */

      if (mtlist[0] != 0)
        {
          /* Reset value to return */

          val = 0.0;

          /* Loop over the mtlist */

          for (i = 0; i < MY_MAXMT; i++)
            {
              /* Break loop if mtlist contains no more values */

              if ((mt = mtlist[i]) == 0)
                break;

              /* Break loop if the rest of the MTs are larger than the ones */
              /* scored for sensitivity */

              if (mt > (long)RDB[sens + SENS_MAX_MT])
                break;

              /* Get reaction index from mt index array and add to sum if partial */
              /* mt was scored during calculation */

              if ((ireaNew = (long)RDB[mtarr + mt]) > 0)
                {
                  idx = CompressSensLabel(imat, izai, ireaNew, iene, 1.0);
                  if (gen >= 0)
                    val = val + BufVal(stp, gen, idx);
                  else
                    val = val + BufVal(stp, idx);
                }
            }

          /* Return value*/

          return val;
        }

    }

  /* Normally just compress the label and return value from the buffer */

  idx = CompressSensLabel(imat, izai, irea, iene, 1.0);
  if (gen >= 0)
    val = BufVal(stp, gen, idx);
  else
    val = BufVal(stp, idx);

  return val;
}
