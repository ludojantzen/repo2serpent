/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : elasticscattering.c                            */
/*                                                                           */
/* Created:       2011/02/28 (JLe)                                           */
/* Last modified: 2018/06/11 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Handles elastic scattering for neutrons                      */
/*                                                                           */
/* Comments: - From Serpent 1.1.13 (10.8.2010)                               */
/*                                                                           */
/*           - Noi nuclide temperaturet vois olla datassa jo kT:nä           */
/*                                                                           */
/*           - Vai pitäisikö kT hakea vasta targetvelocity.c:ssä?            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ElasticScattering:"

/*****************************************************************************/

void ElasticScattering(long mat, long rea, long part, double *E,
                       double *u, double *v, double *w, long id)
{
  double awr, muc, Etot;
  double V, Vx, Vy, Vz, Vt, Vtx, Vty, Vtz, Px, Py, Pz;
  double Vcx, Vcy, Vcz, kT, val;
  double pdf, LMS[7] = {0.0}, PS[7] = {0.0}, SS[7] = {0.0};
  long nuc, ptr, loc0;
  long gidx, iene, imat, izai, imu;
  long i, maxi;


  /***************************************************************************/

  /***** Get initial values and check ****************************************/

  /* Check reaction and material pointers */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Pointer to nuclide */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Check initial energy */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);

  /* Get temperature */

  if(RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE)
    kT = GetTemp(mat, id)*KELVIN;
  else
    kT = RDB[nuc + NUCLIDE_TEMP]*KELVIN;

  /* Get awr */

  awr = RDB[nuc + NUCLIDE_AWR];

  /* Initial neutron velocity */

  V = sqrt(*E);

  Vx = *u*V;
  Vy = *v*V;
  Vz = *w*V;

  /* Sample initial target velocity */

  TargetVelocity(rea, *E, &Vtx, &Vty, &Vtz, *u, *v, *w, kT, id);
  Vt = sqrt(Vtx*Vtx + Vty*Vty + Vtz*Vtz);

  /***************************************************************************/

  /***** Remember some values before the collision ***************************/

  /* Total energy */

  if ((Etot = *E + awr*Vt*Vt) < ZERO)
    Die(FUNCTION_NAME, "Error in energy");

  /* Total momentum */

  Px = Vx + awr*Vtx;
  Py = Vy + awr*Vty;
  Pz = Vz + awr*Vtz;

  /* Check */

  if (Px*Px + Py*Py + Pz*Pz < ZERO)
    Die(FUNCTION_NAME, "Error in momentum");

  /***************************************************************************/

  /***** Transformation to C-frame *******************************************/

  /* Calculate velocity of centre-of-mass */

  Vcx = (Vx + awr*Vtx)/(awr + 1.0);
  Vcy = (Vy + awr*Vty)/(awr + 1.0);
  Vcz = (Vz + awr*Vtz)/(awr + 1.0);

  /* Neutron velocities in C-frame */

  Vx = Vx - Vcx;
  Vy = Vy - Vcy;
  Vz = Vz - Vcz;

  V = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);

  /* Target velocities in C-frame */

  Vtx = Vtx - Vcx;
  Vty = Vty - Vcy;
  Vtz = Vtz - Vcz;

  Vt = sqrt(Vtx*Vtx + Vty*Vty + Vtz*Vtz);

  /***************************************************************************/

  /***** Scattering in C-frame ***********************************************/

  /***** Sens *****************************************************************/

  /* Sample scattering cosine in C-frame. */

  muc = SampleMu(rea, -1, *E, LMS, &pdf, id);

  /* Calculate contribution to Legendre moments for sensitivity */

  if (((loc0 = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR) &&
      ((long)RDB[loc0 + SENS_PERT_FLAGS] & SENS_PERT_FLAG_SCATT_MOM))
    {
      /* Calculate contribution of sampled angle to Legendre moments */

      PS[0] = muc;
      PS[1] = 0.5*(3.0*muc*muc - 1.0);
      PS[2] = 0.5*(5.0*muc*muc*muc - 3.0*muc);
      PS[3] = (35.0*muc*muc*muc*muc - 30.0*muc*muc + 3.0)/8.0;
      PS[4] = (63.0*muc*muc*muc*muc*muc - 70.0*muc*muc*muc + 15*muc)/8.0;
      PS[5] = (231.0*muc*muc*muc*muc*muc*muc - 315.0*muc*muc*muc*muc + 105*muc*muc - 5.0)/16.0;
      PS[6] = (429.0*muc*muc*muc*muc*muc*muc*muc - 693.0*muc*muc*muc*muc*muc + 315*muc*muc*muc - 35.0*muc)/16.0;

      /* Get maximum scattering moment */

      maxi = (long)RDB[loc0 + SENS_MAX_SCATT_MOM];

      /* Calculate product between sampled angle Legendre moments */
      /* and those of the full scattering distribution at energy E */

      for (i = 0; i < maxi; i++)
        SS[i] = (0.5*(2.0*((double)i + 1.0) + 1.0)*PS[i]*LMS[i])/pdf;

      /* Find sub-indices for sensitivity label */

      FindSensIndices(mat, rea, *E, &imat, &izai, &iene);

      /* Create events */

      if ((izai >= 0) && (imat >= 0))
        {
          /* Get pointer to index list */

          ptr = (long)RDB[loc0 + SENS_PTR_PERT_INDICES];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Store events for each scattering moment */

          for (i = 0; i < maxi; i++)
            {
              gidx = CompressSensLabel(imat, izai,
                                       (long)RDB[ptr + ELA_P1_IDX + i],
                                       iene, COLL_HIT);
              StoreSensEvent(part, gidx, *E, SS[i], id);

              /* Score custom perturbations for Legendre moments */

              ScorePerturbations(part, mat, rea, (long)RDB[ptr + ELA_P1_IDX + i], *E, (double)COLL_HIT, SS[i], id);
            }
        }
    }

  /* Perturbation of scattering cosine */

  if (((loc0 = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR) &&
      ((long)RDB[loc0 + SENS_PERT_FLAGS] & SENS_PERT_FLAG_ELA_MU))
    {
      /* Perturbation of elastic scattering cosine */

      /* Get pointer to index list */

      ptr = (long)RDB[loc0 + SENS_PTR_PERT_INDICES];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Find indices */

      FindSensIndices(mat, rea, *E, &imat, &izai, &iene);

      /* Only score these for included materials and zais */

      if ((imat >= 0) && (izai >= 0))
        {
          /* Find index for scattering cosine/angle */

          imu = FindSensMuIndex(muc);

          /* Compress the label */

          gidx = CompressSensLabel(imat, izai, (long)RDB[ptr + ELA_MU_IDX], iene, COLL_HIT);

          /* Create event for accepted scattering cosine */

          StoreSensEvent(part, gidx, *E, muc, id);

          /* Sample rejected scattering cosine */

          val = SampleMu(rea, -1, *E, NULL, NULL, id);

          /* Find index for scattering cosine/angle */

          imu = FindSensMuIndex(val);

          /* Create event for rejected scattering cosine */

          StoreSensEvent(part, gidx*COLL_MISS, *E, val, id);
        }
    }

  /****************************************************************************/

  /* Calculate direction cosines */

  *u = Vx/V;
  *v = Vy/V;
  *w = Vz/V;

  /* Sanity check for mu and direction vectors (for NAN's etc.) */

  CheckValue(FUNCTION_NAME, "muc", "", muc, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "u", "", *u, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "v", "", *v, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "w", "", *w, -1.01, 1.01);

  /* Rotate */

  AziRot(muc, u, v, w, id);

  /* Velocities after collision. */

  Vx = *u*V;
  Vy = *v*V;
  Vz = *w*V;

  Vtx = -*u*Vt;
  Vty = -*v*Vt;
  Vtz = -*w*Vt;

  /***************************************************************************/

  /***** Transformation back to L-frame **************************************/

  /* Neutron velocities in L-frame */

  Vx = Vx + Vcx;
  Vy = Vy + Vcy;
  Vz = Vz + Vcz;

  V = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);

  /* Target velocities in L-frame */

  Vtx = Vtx + Vcx;
  Vty = Vty + Vcy;
  Vtz = Vtz + Vcz;

  Vt = sqrt(Vtx*Vtx + Vty*Vty + Vtz*Vtz);

  /* Set neutron energy */

  *E = (Vx*Vx + Vy*Vy + Vz*Vz);

  /* Set direction cosines */

  *u = Vx/V;
  *v = Vy/V;
  *w = Vz/V;

  /***************************************************************************/

  /***** Check conservation of energy and momentum ***************************/

#ifdef DEBUG

  /* Check relative change in total energy */

  Etot = 1.0 - (*E + awr*Vt*Vt)/Etot;
  CheckValue(FUNCTION_NAME, "E", "", Etot, -1E-5, 1E-4);

  /* Check relative change in momentum components */

  if (Px > 0.0)
    {
      Px = 1.0 - (Vx + awr*Vtx)/Px;
      CheckValue(FUNCTION_NAME, "Px", "", Px, -1E-3, 1E-3);
    }

  if (Py > 0.0)
    {
      Py = 1.0 - (Vy + awr*Vty)/Py;
      CheckValue(FUNCTION_NAME, "Py", "", Py, -1E-3, 1E-3);
    }

  if (Pz > 0.0)
    {
      Pz = 1.0 - (Vz + awr*Vtz)/Pz;
      CheckValue(FUNCTION_NAME, "Pz", "", Pz, -1E-3, 1E-3);
    }

#endif

  /***************************************************************************/
}

/*****************************************************************************/
