/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : storehistorypoint.c                            */
/*                                                                           */
/* Created:       2011/05/13 (JLe)                                           */
/* Last modified: 2019/04/04 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Adds point in particle track in event array.                 */
/*                                                                           */
/* Comments: - Fixed source -laskussa pitää jotenkin huomioida se että       */
/*             lähdeneutroneilla / protoneilla ei oo historiaa               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StoreHistoryPoint:"

/*****************************************************************************/

void StoreHistoryPoint(long part, long mat, long rea, double x, double y,
                       double z, double u, double v, double w, double E,
                       double t, double wgt, double flx, long trk, long id)
{
  long ptr, flags, ok;

  /* Check if events are recorded */

  if ((flags = (long)RDB[DATA_EVENT_RECORD_FLAGS]) == 0)
    return;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Reset OK flag */

  ok = NO;

  /* Types for plotter */

  if ((trk == TRACK_END_STRT) || (trk == TRACK_END_COLL) ||
      (trk == TRACK_END_LEAK) || (trk == TRACK_END_SURF) ||
      (trk == TRACK_END_TCUT) || (trk == TRACK_END_ECUT) ||
      (trk == TRACK_END_WCUT) || (trk == TRACK_END_BC) ||
      (trk == TRACK_END_WWIN))
    if (flags & RECORD_EVENT_PLOTTER)
      ok = YES;

  /* Types for importance */

  if (flx > 0.0)
    if ((trk == TRACK_END_COLL) || (trk == TRACK_END_VIRT))
      if (flags & RECORD_EVENT_IMPORTANCE)
        ok = YES;

  /* Types for IFP (event is added in fission.c) */

  if (trk == TRACK_END_FISS)
    if ((flags & RECORD_EVENT_IFP) || (flags & RECORD_EVENT_SENS))
      ok = NO;

  /* Check ok flag */

  if (ok == NO)
    return;

  /* New event from bank */

  ptr = EventFromBank(part, id);
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Put type */

  WDB[ptr + EVENT_TYPE] = (double)trk;

  /* Put coordinates */

  CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);

  WDB[ptr + EVENT_X] = x;
  WDB[ptr + EVENT_Y] = y;
  WDB[ptr + EVENT_Z] = z;

  /* Put energy and time */

  CheckValue(FUNCTION_NAME, "E", "", E, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "t", "", t, 0.0, INFTY);

  WDB[ptr + EVENT_E] = E;
  WDB[ptr + EVENT_T] = t;

  /* Weight and flux */

  WDB[ptr + EVENT_WGT] = wgt;
  WDB[ptr + EVENT_FLX] = flx;

  /* Put material */

  WDB[ptr + EVENT_PTR_MAT] = (double)mat;

  /***************************************************************************/
}

/*****************************************************************************/
