/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : allocprecdet.c                                 */
/*                                                                           */
/* Created:       2015/11/03 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Sets up some detectors for implicit treatment of delayed     */
/*              neutron emission                                             */
/*                                                                           */
/* Comments: -The group detector is allocated only at processprecdet.c       */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AllocPrecDet:"

/*****************************************************************************/


void AllocPrecDet()
{
  long loc0, loc1, ptr;
  char tmpstr[MAX_STR];


  /* Get pointer to precursor detector or return */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  if ((long)RDB[loc0 + PRECDET_PTR_OUT_FNAME] < VALID_PTR)
    return;

#ifdef DNPRINT
  fprintf(outp, "allocprecdet.c-->\n");
#endif

  /*******************************************************/
  /* Create detector for tallying number of live neutrons*/
  /*******************************************************/

  if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /* Create new item */

      loc1 = NewItem(DATA_PTR_DET0, DET_BLOCK_SIZE);

      /* Store this to precursor detector block */

      WDB[loc0 + PRECDET_PTR_LIVE_DET] = (double)loc1;
    }
  else
    {
      /* In dynamic mode we want to save live neutron population */
      /* manually at the end of the simulation */

      loc1 = NewItem(loc0 + PRECDET_PTR_LIVE_DET, DET_BLOCK_SIZE);

      /* Put total number fo bins */

      WDB[loc1 + DET_N_TOT_BINS] = (double)1;

      /* Allocate memory */

      ptr = NewStat("LiveNeutronDet_Stat", 2, 1, 1);

      /* Put pointer */

      WDB[loc1 + DET_PTR_STAT] = (double)ptr;

      /* Allocate memory for history */

      if((long)RDB[loc1 + DET_WRITE_HIS] == 1)
        AllocStatHistory(ptr);

    }

  /* Put name, file name and line number */

  WDB[loc1 + PARAM_PTR_NAME] = RDB[loc0 + PARAM_PTR_NAME];
  WDB[loc1 + PARAM_PTR_FNAME] = RDB[loc0 + PARAM_PTR_FNAME];
  WDB[loc1 + PARAM_LINE] = RDB[loc0 + PARAM_LINE];

  /* Reset volume */

  WDB[loc1 + DET_VOL] = 1.0;

  /* Reset mesh pointer */

  WDB[loc1 + DET_PTR_MESH] = NULLPTR;

  /* Detector name */

  WDB[loc1 + DET_PTR_NAME] = (double)PutText("LiveNeutronDetector");

  /* Set particle type */

  WDB[loc1 + DET_PARTICLE] = (double)PARTICLE_TYPE_NEUTRON;

  /* Detector reaction, allocate memory for structure */

  loc1 = NewItem(loc1 + DET_PTR_RBINS, DET_RBIN_BLOCK_SIZE);

  /* Store MT -15: Flux multiplied with 1/v. */
  /* Should give integral of neutron density */
  /* i.e. total number of live neutrons */

  WDB[loc1 + DET_RBIN_MT] = MT_NEUTRON_DENSITY;

  /* Store material name (void) */

  WDB[loc1 + DET_RBIN_PTR_MAT] = (double)PutText("void");

  /******************************************************/
  /* Create detector for writing live neutrons to file  */
  /******************************************************/

  /* Create new item */

  if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      loc1 = NewItem(DATA_PTR_DET0, DET_BLOCK_SIZE);

      /* Store this to precursor detector block */

      WDB[loc0 + PRECDET_PTR_FILE_DET] = (double)loc1;
    }
  else
    {
      /* In dynamic mode we want to save live neutrons manually */
      /* at the end of the simulation */

      loc1 = NewItem(loc0 + PRECDET_PTR_FILE_DET, DET_BLOCK_SIZE);

      /* Allocate memory for buffer */

      ptr = ReallocMem(DATA_ARRAY,
                       (long)RDB[DATA_SRC_FILE_BUF_SIZE]*SRC_BUF_BLOCK_SIZE);

      /* Put pointer and buffer size */

      WDB[loc1 + DET_WRITE_PTR_BUF] = (double)ptr;
      WDB[loc1 + DET_WRITE_BUF_SZ] = RDB[DATA_SRC_FILE_BUF_SIZE];

      /* Reset index */

      WDB[loc1 + DET_WRITE_BUF_IDX] = 0.0;
    }

  /* Put name, file name and line number */

  WDB[loc1 + PARAM_PTR_NAME] = RDB[loc0 + PARAM_PTR_NAME];
  WDB[loc1 + PARAM_PTR_FNAME] = RDB[loc0 + PARAM_PTR_FNAME];
  WDB[loc1 + PARAM_LINE] = RDB[loc0 + PARAM_LINE];

  /* Reset volume */

  WDB[loc1 + DET_VOL] = 1.0;

  /* Reset mesh pointer */

  WDB[loc1 + DET_PTR_MESH] = NULLPTR;

  /* Detector name */

  WDB[loc1 + DET_PTR_NAME] = (double)PutText("FileNeutronDetector");

  /* Set particle type */

  WDB[loc1 + DET_PARTICLE] = (double)PARTICLE_TYPE_NEUTRON;

  /* Print live neutron filename to tmpstr */

  sprintf(tmpstr, "%s.live", GetText(loc0 + PRECDET_PTR_OUT_FNAME));

  /* Store output file name to detector */

  WDB[loc1 + DET_WRITE_PTR_FILE] =
    (double)PutText(tmpstr);

  /* Store fraction of points to store */

  WDB[loc1 + DET_WRITE_PROB] = RDB[loc0 + PRECDET_SAVE_FRAC_LIVE];

  /* Set binary output on */

  WDB[loc1 + DET_WRITE_BINARY] = (double)YES;

  /* Remove old file */

  remove(GetText(loc1 + DET_WRITE_PTR_FILE));

  /*****************************************************/
  /* Create detector for writing precursors to file    */
  /*****************************************************/

  /* Create new detector, do not create it to DET0 list */
  /* as its scoring will be handled manually */

  loc1 = NewItem(loc0 + PRECDET_PTR_PREC_DET, DET_BLOCK_SIZE);

  /* Put name, file name and line number */

  WDB[loc1 + PARAM_PTR_NAME] = RDB[loc0 + PARAM_PTR_NAME];
  WDB[loc1 + PARAM_PTR_FNAME] = RDB[loc0 + PARAM_PTR_FNAME];
  WDB[loc1 + PARAM_LINE] = RDB[loc0 + PARAM_LINE];

  /* Reset volume */

  WDB[loc1 + DET_VOL] = 1.0;

  /* Reset mesh pointer */

  WDB[loc1 + DET_PTR_MESH] = NULLPTR;

  /* Detector name */

  WDB[loc1 + DET_PTR_NAME] = (double)PutText("FilePrecursorDetector");

  /* Set particle type */

  WDB[loc1 + DET_PARTICLE] = (double)PARTICLE_TYPE_NEUTRON;

  /* Print live neutron filename to tmpstr */

  sprintf(tmpstr, "%s.precpoints", GetText(loc0 + PRECDET_PTR_OUT_FNAME));

  /* Store output file name to detector */

  WDB[loc1 + DET_WRITE_PTR_FILE] =
    (double)PutText(tmpstr);

  /* Set binary output on */

  WDB[loc1 + DET_WRITE_BINARY] = (double)YES;

  /* Store fraction of points to store */

  WDB[loc1 + DET_WRITE_PROB] = RDB[loc0 + PRECDET_SAVE_FRAC_PREC];

  /* Allocate memory for buffer */

  ptr = ReallocMem(DATA_ARRAY,
                   (long)RDB[DATA_SRC_FILE_BUF_SIZE]*SRC_BUF_BLOCK_SIZE);

  /* Put pointer and buffer size */

  WDB[loc1 + DET_WRITE_PTR_BUF] = (double)ptr;
  WDB[loc1 + DET_WRITE_BUF_SZ] = RDB[DATA_SRC_FILE_BUF_SIZE];

  /* Reset index */

  WDB[loc1 + DET_WRITE_BUF_IDX] = 0.0;

  /* Remove old file */

  remove(GetText(loc1 + DET_WRITE_PTR_FILE));

#ifdef DNPRINT
  fprintf(outp, "<-- allocprecdet.c\n\n");
#endif

}
