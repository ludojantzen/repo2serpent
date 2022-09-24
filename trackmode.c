/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : trackmode.c                                    */
/*                                                                           */
/* Created:       2012/05/23 (JLe)                                           */
/* Last modified: 2020/04/09 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Selects delta- or surface-tracking for next path length      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TrackMode:"

/*****************************************************************************/

long TrackMode(long part, long mat, double E, double totxs, double majorant,
               long type, long id)
{
  long ptr, tmode;
  double trsh;

  /* Get tracking mode from domain decomposition */

  tmode = (long)RDB[part + PARTICLE_DD_TRACK_MODE];

  /* Reset mode */

  WDB[part + PARTICLE_DD_TRACK_MODE] = 0.0;

  /* Check if particle came from different domain */

  if (tmode == TRACK_MODE_DT)
    {
      /* Use delta-tracking for next (mpiid:tä ei resetoida jotta */
      /* MoveDT():ssä päästään suoraan rejektiosämpläykseen) */

      return TRACK_MODE_DT;
    }
  else if (tmode == TRACK_MODE_ST)
    {
      /* Check MPI id */

      if ((long)RDB[part + PARTICLE_MPI_ID] == mpiid)
        Die(FUNCTION_NAME, "Error in MPI id");

      /* Set the new MPI id */

      WDB[part + PARTICLE_MPI_ID] = mpiid;
    }
  else if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
           ((long)RDB[part + PARTICLE_MPI_ID] != mpiid))
    Die(FUNCTION_NAME, "Error in MPI id");

  /* Check forced dt by flag */

  ptr = (long)RDB[DATA_DT_ENFORCE_NEXT_TRACK];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

  if (GetPrivateData(ptr, id) == YES)
    {
      /* Reset mode */

      PutPrivateData(ptr, NO, id);

      /* Check if particle came from dd using st */

      if (tmode == TRACK_MODE_ST)
        Die(FUNCTION_NAME, "WTF?");

      /* Use DT */

      return TRACK_MODE_DT;
    }

  /* Check forced dt in material */

  if (mat > VALID_PTR)
    {
      if ((long)RDB[mat + MATERIAL_DT_MODE] == DT_MAT_FORCE)
        return TRACK_MODE_DT;
      else if ((long)RDB[mat + MATERIAL_DT_MODE] == DT_MAT_BLOCK)
        return TRACK_MODE_ST;
    }

  /* Check cross sections */

  CheckValue(FUNCTION_NAME, "totxs", "", totxs, 0.0, INFTY);
  CheckValue(FUNCTION_NAME, "majorant", "", majorant, ZERO, INFTY);

  /* Compare total to majorant (NOTE: tää voi ylittyä jos moodi on < 2 */
  /* ja grid thinning käytössä). */

  if ((mat > VALID_PTR) && (totxs/majorant - 1.0 > 1E-6))
    {
      /* Check particle type */

      if ((long)RDB[part + PARTICLE_TYPE] == PARTICLE_TYPE_NEUTRON)
        {
          /* Check reaction lists */

          CheckReaListSum(mat, PARTICLE_TYPE_NEUTRON, E, NO, id);

          /* Exit */

          Die(FUNCTION_NAME,
              "Neutron total exceeds majorant (E = %E, f = %E, mat = %s)",
              E, totxs/majorant - 1.0, GetText(mat + MATERIAL_PTR_NAME));
        }
      else
        {
          /* Check reaction lists */

          CheckReaListSum(mat, PARTICLE_TYPE_GAMMA, E, NO, id);

          /* Exit */

          Die(FUNCTION_NAME,
              "Photon total exceeds majorant (E = %E, f = %E, mat = %s)",
              E, totxs/majorant - 1.0, GetText(mat + MATERIAL_PTR_NAME));
        }
    }

  /* Get delta-tracking threshold */

  if (type == PARTICLE_TYPE_NEUTRON)
    trsh = RDB[DATA_DT_NTHRESH];
  else
    trsh = RDB[DATA_DT_PTHRESH];

  /* Check mode */

  if (trsh == 0.0)
    return TRACK_MODE_DT;
  else if (trsh == 1.0)
    return TRACK_MODE_ST;
  else if (totxs/majorant < trsh)
    return TRACK_MODE_ST;
  else
    return TRACK_MODE_DT;
}

/*****************************************************************************/
