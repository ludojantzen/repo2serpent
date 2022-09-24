/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : detidx.c                                       */
/*                                                                           */
/* Created:       2011/07/11 (JLe)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns stat index for detectors                             */
/*                                                                           */
/* Comments: 2017/02/22: Moved FET indexing to separate routine              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DetIdx:"

/*****************************************************************************/

long DetIdx(long det, long ebin, long ubin, long cbin, long mbin, long lbin,
            long zbin, long ybin, long xbin, long tbin)
{
  long ne, nu, nc, nm, nl, nz, ny, nx, nt, nmax, idx, ptr;

  /* Get number of bins */

  ne = (long)RDB[det + DET_N_EBINS];
  nu = (long)RDB[det + DET_N_UBINS];
  nc = (long)RDB[det + DET_N_CBINS];
  nm = (long)RDB[det + DET_N_MBINS];
  nl = (long)RDB[det + DET_N_LBINS];
  nt = (long)RDB[det + DET_N_TBINS];

  /* Mesh bins */

  if ((ptr = (long)RDB[det + DET_PTR_MESH]) > VALID_PTR)
    {
      nx = (long)RDB[ptr + MESH_N0];
      ny = (long)RDB[ptr + MESH_N1];
      nz = (long)RDB[ptr + MESH_N2];
    }
  else
    {
      nx = 1;
      ny = 1;
      nz = 1;
    }

  /* Check FET */

  if ((ptr = (long)RDB[det + DET_FET_PTR_PARAMS]) > VALID_PTR)
    Die(FUNCTION_NAME, "Pointer to FET parameters");

  /* Check bins */

  if ((ebin < 0) || (ebin > ne))
    Die(FUNCTION_NAME, "Invalid ebin");
  if ((ubin < 0) || (ubin > nu))
    Die(FUNCTION_NAME, "Invalid ubin");
  if ((cbin < 0) || (cbin > nc))
    Die(FUNCTION_NAME, "Invalid cbin");
  if ((mbin < 0) || (mbin > nm))
    Die(FUNCTION_NAME, "Invalid mbin");
  if ((lbin < 0) || (lbin > nl))
    Die(FUNCTION_NAME, "Invalid lbin");
  if ((zbin < 0) || (zbin > nz))
    Die(FUNCTION_NAME, "Invalid zbin");
  if ((ybin < 0) || (ybin > ny))
    Die(FUNCTION_NAME, "Invalid ybin");
  if ((xbin < 0) || (xbin > nx))
    Die(FUNCTION_NAME, "Invalid xbin");
  if ((tbin < 0) || (tbin > nt))
    Die(FUNCTION_NAME, "Invalid tbin");

  /* Calculate index */

  nmax = 1;
  idx = 0;

  idx = idx + ebin*nmax;
  nmax = nmax*ne;

  idx = idx + ubin*nmax;
  nmax = nmax*nu;

  idx = idx + cbin*nmax;
  nmax = nmax*nc;

  idx = idx + mbin*nmax;
  nmax = nmax*nm;

  idx = idx + lbin*nmax;
  nmax = nmax*nl;

  idx = idx + xbin*nmax;
  nmax = nmax*nx;

  idx = idx + ybin*nmax;
  nmax = nmax*ny;

  idx = idx + zbin*nmax;
  nmax = nmax*nz;

  idx = idx + tbin*nmax;
  nmax = nmax*nt;

  /* Check index */

  if ((idx < 0) || (idx > (long)RDB[det + DET_N_TOT_BINS]))
    Die(FUNCTION_NAME, "Error in idx");

  /* Return index */

  return idx;
}

/*****************************************************************************/
