/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : SampleMeshDelnu.c                              */
/*                                                                           */
/* Created:       2015/09/29 (VVa)                                           */
/* Last modified: 2016/12/19 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Samples delayed neutrons from precursor concentrations at    */
/*              the beginning of an interval and adds them to que            */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SampleMeshDelnu:"

#define MAX_DELNU_RESAMPLE 10000

/*****************************************************************************/

void SampleMeshDelnu(long id, long np, long idx)
{
  long ptr, pos, seed;
  long loc0, loc1, ng, n0, n1, n2;
  long gbin, xbin, ybin, zbin;
  long prec, precptr, reaptr, lamptr, vplist;
  long cell, mat, rea, erg, new, idx0, i, j;
  double E, mu, u, v, w, wgt2, val;
  double dt, lambda, Wemit, Wall;
  double t0, t1, temit, t, x, y, z;
  double min0, max0, min1, max1, min2, max2;
#ifdef OLDSAMPLING
  double Wmax;
  long part, ibin, tbin, stp;
#endif

  /***************************************************************************/

  /* Get pointer to precursor detector */

  loc0 = (long)RDB[DATA_PTR_PREC_DET];

  /* Close buffer for reading */

  WDB[DATA_BUF_REDUCED] = (double)YES;

  /* Get pointer to source */

  pos = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(pos)", DATA_ARRAY, pos);

  /* Check that source is empty */

  if (ListSize(pos) != 1)
    Die(FUNCTION_NAME, "Source is not empty");

  /* Get pointer to precursor list */

  precptr = (long)RDB[loc0 + PRECDET_PTR_PREC_ARRAY];
  CheckPointer(FUNCTION_NAME, "(precptr)", DATA_ARRAY, precptr);

  /* Get pointer to reaction list */

  reaptr = (long)RDB[loc0 + PRECDET_PTR_REA_ARRAY];
  CheckPointer(FUNCTION_NAME, "(reaptr)", DATA_ARRAY, reaptr);

#ifdef OLDSAMPLING

  /* Get current time bin */

  tbin = (long)RDB[DATA_DYN_TB];

#endif

  /* Get time interval limits */

  t0 = RDB[DATA_TIME_CUT_TMIN];
  t1 = RDB[DATA_TIME_CUT_TMAX];

  /* Calculate time interval length */

  dt = t1 - t0;

  /* Calculate weight of neutrons to emit */

  wgt2 = RDB[loc0 + PRECDET_W_EMIT]/RDB[loc0 + PRECDET_N_EMIT];

  Wall = RDB[loc0 + PRECDET_W_EMIT]*RDB[DATA_NORM_COEF_N];

  /* Get number of group bins */

  ng = (long)RDB[loc0 + PRECDET_NG];

  /* Get pointer to mesh */

  ptr = (long)RDB[loc0 + PRECDET_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get number of mesh bins */

  n0 = (long)RDB[ptr + MESH_N0];
  n1 = (long)RDB[ptr + MESH_N1];
  n2 = (long)RDB[ptr + MESH_N2];

  /* Get mesh limits */

  min0 = RDB[ptr + MESH_MIN0];
  max0 = RDB[ptr + MESH_MAX0];
  min1 = RDB[ptr + MESH_MIN1];
  max1 = RDB[ptr + MESH_MAX1];
  min2 = RDB[ptr + MESH_MIN2];
  max2 = RDB[ptr + MESH_MAX2];

  /* Get pointer to lambda list */

  lamptr = (long)RDB[loc0 + PRECDET_PTR_LAM_ARRAY];
  CheckPointer(FUNCTION_NAME, "(lamptr)", DATA_ARRAY, lamptr);

#ifdef OLDSAMPLING

  /* Get maximum bin value to emit */
  /* Used as majorant for rejection sampling of bins */

  Wmax = RDB[loc0 + PRECDET_MAX_WGT];

  /* Get pointer to statistics */

  stp = (long)RDB[loc0 + PRECDET_PTR_STAT];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);

#endif

  /* Init random number sequence */

  seed = ReInitRNG(idx);
  SEED[id*RNG_SZ] = seed;

  /* Loop until we have sampled toemit neutrons */

  for (i = 0; i < MAX_DELNU_RESAMPLE ; i++)
    {

      /* NB: It might be better to first sample coordinates */
      /* check for fissile material and then get the bin and reject based on mesh bin */
      /* NB: It is also possible to sample the group bin and mesh bin from the unit   */
      /* interval and then sample coordinates from that mesh bin */
#ifdef OLDSAMPLING
      /* Sample group bin to emit from */

      gbin = (long)(drand48()*ng);

      /* Sample mesh bin to emit from */

      xbin = (long)(drand48()*n0);
      ybin = (long)(drand48()*n1);
      zbin = (long)(drand48()*n2);

      ibin = xbin + ybin*n0 + zbin*n0*n1;

      /*
        fprintf(output, "Sampled bin: group %ld mesh %ld\n", gbin, ibin);
      */

      /* Get lambda for group */

      lambda = RDB[lamptr + gbin];

      /* Get pointer for group precursor block */

      prec = (long)RDB[precptr + gbin];

      /* Get number of precursors */

      val = BufVal(stp, tbin, gbin, ibin);

      /* Calculate physical number to emit */

      Wemit = (1-exp(-lambda*dt))*val;

      /* Divide by normalization coefficient to get weight to emit */

      Wemit = Wemit/RDB[DATA_NORM_COEF_N];

      /*
        fprintf(output, "Number of precursors %E, to emit %E, max %E, P = %f\n", val, Wemit, Wmax, Wemit/Wmax);
      */

      /* Reject bin with probability 1 - Wemit/Wmax */

      if (RandF(id) < 1 - Wemit/Wmax)
        continue;
#else

      /* Sample bin to emit in the same way e.g. reaction sampling is done   */
      /* By sampling a fraction of the total emission and then looping over  */
      /* bins subtracting the emission from that bin until we get a negative */
      /* value */

      /* Get pointer to value pair list */

      if ((vplist = (long)RDB[loc0 + PRECDET_PTR_MESH_LIST]) < VALID_PTR)
        Die(FUNCTION_NAME, "Value pair list does not exist");

      /* Sample random number to sample bin with */

      Wemit = RandF(id);

      /* Loop over valuepair list until correct bin is found */

      while (vplist > VALID_PTR)
        {
          /* Get value to emit from this bin */

          val = RDB[vplist + VALUE_PAIR_VAL2];

          /* Subtract this bins fraction from sampled fraction */

          Wemit = Wemit - val/Wall;

          /* If we went negative, this is the bin to emit from */

          if (Wemit < 0)
            break;

          /* Get next bin */

          vplist = NextItem(vplist);
        }

      /* Check if we looped over the whole list without sampling a bin */

      if (vplist < VALID_PTR)
        {
          /* Print warning */

          Warn(FUNCTION_NAME, "Failed sampling mesh bin");

          /* Cycle main sampling loop */

          continue;
        }

      /* Get index of current bin */
      /* idx = xbin + ybin*n0 + zbin*n0*n1 + gbin*n0*n1*n2  */

      idx = (long)RDB[vplist + VALUE_PAIR_VAL1];

      /* We'll still need to figure out the group (mainly for lambda) */
      /* and spatial bin (mainly for sampling the coordinates) */

      /****************************** !!!!!!! *************************/
      /* This will take some time but save some memory since we only  */
      /* need to store idx (and not gbin, xbin, ybin, zbin) in vplist */

      /* Calculate group bin */

      for (gbin = 0; gbin < ng; gbin++)
        {
          if ((gbin + 1)*n0*n1*n2 > idx)
            break;
        }

      /* Check group bin */

      if (gbin >= ng)
        Die(FUNCTION_NAME, "Sampled invalid group bin: %ld", gbin);

      /* Subtract the group specific part from idx and check */

      if ((idx = idx - gbin*n0*n1*n2) < 0)
        Die(FUNCTION_NAME, "Sampling wrong");

      /* Calculate zbin */

      for (zbin = 0; zbin < n2; zbin++)
        {
          if ((zbin + 1)*n0*n1 > idx)
            break;
        }

      /* Check zbin */

      if (zbin >= n2)
        Die(FUNCTION_NAME, "Sampled invalid zbin: %ld", zbin);

      /* Subtract the z-specific part from idx and check */

      if ((idx = idx - zbin*n0*n1) < 0)
        Die(FUNCTION_NAME, "Sampling wrong");

      /* Calculate ybin */

      for (ybin = 0; ybin < n1; ybin++)
        {
          if ((ybin + 1)*n0 > idx)
            break;
        }

      /* Check ybin */

      if (ybin >= n1)
        Die(FUNCTION_NAME, "Sampled invalid ybin: %ld", ybin);

      /* The remainder should be xbin */

      if ((xbin = idx - ybin*n0) < 0)
        Die(FUNCTION_NAME, "Sampling wrong");

      if (xbin >= n0)
        Die(FUNCTION_NAME, "Sampled invalid xbin: %ld", xbin);

      /* Get lambda for group */

      lambda = RDB[lamptr + gbin];

      /* Get pointer for group precursor block */

      prec = (long)RDB[precptr + gbin];

#endif
      /* Passed rejection sampling, emitting from this bin */

      /* Sample emission time */
      /* (time from this moment to emission) */

      temit = -1/lambda*log(1 - (1 - exp(-lambda*dt))*drand48());

      /* Calculate absolute time of emission */

      t = t0 + temit;

      /* Sample emission coordinates */
      /* This may fail if the spatial domain of the bin is mostly */
      /* non fissionable, in that case a new bin is sampled */

      for (j = 0; j < MAX_DELNU_RESAMPLE ; j++)
        {
          /* sample coordinates */

          x = min0 + (double)xbin*(max0 - min0)/(double)n0 +
            drand48()*(max0 - min0)/(double)n0;
          y = min1 + (double)ybin*(max1 - min1)/(double)n1 +
            drand48()*(max1 - min1)/(double)n1;
          z = min2 + (double)zbin*(max2 - min2)/(double)n2 +
            drand48()*(max2 - min2)/(double)n2;

          /* Get cell (if geometry boundary is something else than a square */
          /* we can be outside the geometry */

          if ((cell = WhereAmI(x, y, z, u, v, w, id)) < VALID_PTR)
            {
              /* If there is no cell at these coordinates, resample x,y,z */
              continue;
            }

          /* Get material pointer */

          if ((mat = (long)RDB[cell + CELL_PTR_MAT]) < VALID_PTR)
            {
              /* If there is no material at these coordinates, resample x,y,z */
              continue;
            }

          /* Check that material is fissile */

          if (RDB[mat + MATERIAL_INI_FMASS] == 0.0)
            {
              /* If the material at these coordinates is not fissile, resample  */
              continue;
            }

          /* Coordinates in fissile material -> move on */
          break;
        }

      if (j >= MAX_DELNU_RESAMPLE)
        {
          Warn(FUNCTION_NAME, "Spatial sampling failed");
          continue;
        }

      /* Sample emission energy */
      /* Get pointer to precursor group*/

      loc1 = (long)RDB[precptr + gbin];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Pointer to reaction */

      rea = (long)RDB[reaptr + gbin];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Pointer to precursor groups energy grid */

      erg = (long)RDB[loc1 + PREC_PTR_ERG];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Reset mu (scattering cosine) */

      mu = 0.1;

      /* Sample energy */
      /* Using 1e-6 MeV as incoming neutron energy */
      /* Sampling is probably independent from incoming neutron */
      /* energy, but it has to be checked at some point */

      SampleENDFLaw(rea, erg, 1e-6, &E, &mu, id);

      /* Check that direction cosine was not changed */
      /* If it was changed, the delayed neutron emission reaction */
      /* actually has some kind of directional data, which would */
      /* be weird */

      if (mu != 0.1)
        Die(FUNCTION_NAME, "mu was changed by SampleENDFLaw");

      /* Sample emission direction cosines isotropically */
      /* The direction for the delayed neutron should not */
      /* depend on the direction of the neutron that caused */
      /* the fission */

      IsotropicDirection(&u, &v, &w, id);

      /*
        fprintf(output, "Sampled some stuff for delayed neutron\n");
        fprintf(output, "t = %E\n", t);
        fprintf(output, "E = %E\n", E);
        fprintf(output, "(x,y,z) = (%E, %E, %E)\n", x, y, z);
        fprintf(output, "(u,v,w) = (%E %E %E)\n", u, v, w);
      */

      /* Set outgoing neutron weight */
      /*
        wgt2 = 1.0;
      */
      /* Create the neutron to be emitted and set the parameters */
      /* Duplicate incident neutron */

      new = FromStack(PARTICLE_TYPE_NEUTRON, id);

      /* Put variables */

      WDB[new + PARTICLE_X] = x;
      WDB[new + PARTICLE_Y] = y;
      WDB[new + PARTICLE_Z] = z;

      WDB[new + PARTICLE_U] = u;
      WDB[new + PARTICLE_V] = v;
      WDB[new + PARTICLE_W] = w;

      WDB[new + PARTICLE_E] = E;
      WDB[new + PARTICLE_WGT] = wgt2;
      WDB[new + PARTICLE_PTR_MAT] = (double)mat;

      /* This is nuclide-wise idx (between 1 and 6 or 1 and 8) */

      WDB[new + PARTICLE_DN_GROUP] = RDB[prec + PREC_IDX];
      WDB[new + PARTICLE_DN_LAMBDA] = lambda;

      /* Store absolute emission time */
      /* PARTICLE_T0 is the birth time of the neutron*/
      /* PARTICLE_T  is the current time of the neutron */

      WDB[new + PARTICLE_T0] = t;
      WDB[new + PARTICLE_T] = t;

      /* Delayed neutron emission time */

      WDB[new + PARTICLE_TD] = temit;

      /* Reset thermalization time */

      WDB[new + PARTICLE_TT] = 0.0;

      /* Get fission matrix index */

      idx0 = FissMtxIndex(mat, id);

      /* Put fission matrix index */

      WDB[new + PARTICLE_FMTX_IDX] = (double)idx0;

      /* Put MPI id */

      WDB[new + PARTICLE_MPI_ID] = (double)mpiid;

      /* Reset generation index */

      WDB[new + PARTICLE_GEN_IDX] = 0.0;

#ifdef ScoretheseAfterBuffersAreOpen

      /* Score mesh plotter */

      ScoreMesh(new, mat, 0.0, -1.0, x, y, z, E, t, wgt, 1.0, id);


      /* Score source rate for fission matrix */

      if (idx0 > -1)
        {
          /* Pointer to matrix */

          ptr = (long)RDB[DATA_PTR_FMTX];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          ptr = (long)RDB[ptr + FMTX_PTR_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Score total */

          AddBuf(wgt, 1.0, ptr, id, -1, 0, idx0);

          /* Score prompt (no delayed in source) */

          AddBuf(wgt, 1.0, ptr, id, -1, 1, idx0);
        }

      /* NOTE: Mesh tallyilla on mahdotonta pysya perassa generaatioista */

      if (RDB[new + PARTICLE_GEN_IDX] == 0.0)
        {
          /* Score initial source weight */

          ptr = (long)RDB[RES_INI_SRC_WGT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, wgt, ptr, id, 0);

          /* Score source rate */

          ptr = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, wgt, ptr, id, 0);
        }

      /* Score source rate in fissile and non-fissile materials */

      if (mat > VALID_PTR)
        {
          if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
            AddBuf1D(1.0, wgt, ptr, id, 1);
          else
            AddBuf1D(1.0, wgt, ptr, id, 2);
        }

#endif

      /* Put particle in que */

      ToQue(new, id);

      /* Break from the resample loop */

      break;

    }

  /***************************************************************************/
}

/*****************************************************************************/
