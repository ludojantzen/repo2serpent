/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : atomicrelaxation.c                             */
/*                                                                           */
/* Created:       2014/07/12 (TKa)                                           */
/* Last modified: 2017/04/24 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Handles atomic relaxation process                            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AtomicRelaxation:"

/*****************************************************************************/

void AtomicRelaxation(long mat, long rea, long part, long ss0, double x,
                      double y, double z, double wgt, double t, double *Edep,
                      long id) {

  long ptd, ptr, ss, ssj, ssk, subj, subk, tr, newp, ssqueuefront, ssqueueback,
      ptretr, ptrletr, ptrsubj, ptrsubk, Nssd, nuc, ptralias1, ptralias2,
      ptrcutoff, ptridxrad, nbr, idx, k;
  double ebi0, etr, letr, EdepTTB, Etrtot, u, v, w, R;
  const double *d2imap, *ntrd, *ebid, *elnd;
  long ssqueue[PHOTON_ZMAX];
  long vacancy[PHOTON_NSS_MAX];
  double Emin = RDB[mat + MATERIAL_PHOTON_ECUT];
  long useTTB = (long)RDB[DATA_PHOTON_USE_TTB] == YES;

  /* Check the initial vacancy subshell index ss0 */
  if (ss0 < 0)
    Die(FUNCTION_NAME, "Negative primary subshell index");

  /* Get nuclide pointer */
  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Pointer to relaxation data */
  ptd = (long)RDB[nuc + NUCLIDE_PTR_RELAX];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  /* Number of shells in the data */
  Nssd = (long)RDB[ptd + RELAX_NSS];
  CheckValue(FUNCTION_NAME, "NSS", "", (double)Nssd, 1, PHOTON_NSS_MAX);

  if (ss0 >= Nssd) {
    /* Relaxation data not found */
    Die(FUNCTION_NAME, "Atomic relaxation data not found for the shell index "
                      "%ld of element Z=%ld", (long)RDB[nuc + NUCLIDE_Z], ss0);
  }

  ptr = (long)RDB[ptd + RELAX_D2IMAP];
  CheckPointer(FUNCTION_NAME, "(D2IMAP)", DATA_ARRAY, ptr);
  d2imap = &RDB[ptr];

  ptr = (long)RDB[ptd + RELAX_NTR];
  CheckPointer(FUNCTION_NAME, "(NTR)", DATA_ARRAY, ptr);
  ntrd = &RDB[ptr];

  ptr = (long)RDB[ptd + RELAX_EBI];
  CheckPointer(FUNCTION_NAME, "(EBI)", DATA_ARRAY, ptr);
  ebid = &RDB[ptr];

  ptr = (long)RDB[ptd + RELAX_ELN];
  CheckPointer(FUNCTION_NAME, "(ELN)", DATA_ARRAY, ptr);
  elnd = &RDB[ptr];

  /* Pointers to electron shell data */
  ptretr = (long)RDB[ptd + RELAX_ETR];
  CheckPointer(FUNCTION_NAME, "(ETR)", DATA_ARRAY, ptretr);

  ptrletr = (long)RDB[ptd + RELAX_LETR];
  CheckPointer(FUNCTION_NAME, "(LETR)", DATA_ARRAY, ptrletr);

  ptridxrad = (long)RDB[ptd + RELAX_IDX_RAD];
  CheckPointer(FUNCTION_NAME, "(IDX_RAD)", DATA_ARRAY, ptridxrad);

  ptrsubj = (long)RDB[ptd + RELAX_SUBJ];
  CheckPointer(FUNCTION_NAME, "(SUBJ)", DATA_ARRAY, ptrsubj);

  ptrsubk = (long)RDB[ptd + RELAX_SUBK];
  CheckPointer(FUNCTION_NAME, "(SUBK)", DATA_ARRAY, ptrsubk);

  ptralias1 = (long)RDB[ptd + RELAX_ALIAS1];
  CheckPointer(FUNCTION_NAME, "(ALIAS1)", DATA_ARRAY, ptralias1);

  ptralias2 = (long)RDB[ptd + RELAX_ALIAS2];
  CheckPointer(FUNCTION_NAME, "(ALIAS2)", DATA_ARRAY, ptralias2);

  ptrcutoff = (long)RDB[ptd + RELAX_CUTOFF];
  CheckPointer(FUNCTION_NAME, "(CUTOFF)", DATA_ARRAY, ptrcutoff);

  /* Binding energy of the initial vacancy subshell */
  ebi0 = ebid[ss0];

  /* Check the lower energy limit */
  if (ebi0 < Emin) {
    /* The whole binding energy is deposited locally  */
    *Edep = ebi0;
    return;
  }

  /* Initialize vacancy array */
  memset(vacancy, 0, PHOTON_NSS_MAX*sizeof(double));

  /* Set the initial vacancy */
  vacancy[ss0] = 1;

  /* Intialize subshell queue */
  ssqueuefront = 0;
  ssqueueback = 0;
  ssqueue[ssqueueback] = ss0;

  /* Total transported energy */
  Etrtot = 0.0;

  while (ssqueuefront <= ssqueueback) {

    ss = ssqueue[ssqueuefront++];

    /* Check the subshell index */
    if (ss == -1) {
      /* Data not available for the subshell (ebi below DATA_PHOTON_EMIN) */
      /* NOTE: Vacancy array is not changed */
      continue;
    }

    /* Decrease primary subshell vacancy */
    vacancy[ss]--;

    /* Check the binding energy */
    if (ebid[ss] < Emin)
      continue;

    while (1) {

      /* Sample transition using Walker's alias method */
      R = RandF(id)*ntrd[ss];
      k = (long)R;
      if (R - (double)k < RDB[(long)RDB[ptrcutoff + ss] + k])
        tr = (long)RDB[(long)RDB[ptralias1 + ss] + k];
      else
        tr = (long)RDB[(long)RDB[ptralias2 + ss] + k];

      CheckValue(FUNCTION_NAME, "tr", "", (double)tr, 0, ntrd[ss] - 1.0);

      subj = (long)RDB[(long)RDB[ptrsubj + ss] + tr];
      subk = (long)RDB[(long)RDB[ptrsubk + ss] + tr];

      ssj = (long)d2imap[subj];
      ssk = (long)d2imap[subk];

      /* Check vacancies */
      /* NOTE: Vacancies are checked only from the available subshells */
      if ((ssj != -1) && (vacancy[ssj] == elnd[ssj]))
        continue;
      else if ((ssk != -1) && (vacancy[ssk] == elnd[ssk]))
        continue;
      else
        break;
    }

    /* Transition energy */
    etr = RDB[(long)RDB[ptretr + ss] + tr];

    /* Logarithm of the transition energy */
    letr = RDB[(long)RDB[ptrletr + ss] + tr];

    if (ssj >= 0) {
      /* Increase secondary subshell vacancy */
      vacancy[ssj]++;

      /* Put secondary subshell to queue */
      ssqueue[++ssqueueback] = ssj;

      CheckValue(FUNCTION_NAME, "ssqueueback", "Maximum number of vacancies "
                 "exceeded", (double)ssqueueback, (double)ssqueuefront,
                 PHOTON_ZMAX-1);
    }

    if (subk == 0) {
      /* Radiative transition */

      if (etr > Emin) {

        /* Sample direction isotropically */
        IsotropicDirection(&u, &v, &w, id);

        /* Create a new photon */
        newp = DuplicateParticle(part, id);

        /* Put variables */
        WDB[newp + PARTICLE_X] = x;
        WDB[newp + PARTICLE_Y] = y;
        WDB[newp + PARTICLE_Z] = z;
        WDB[newp + PARTICLE_WGT] = wgt;
        WDB[newp + PARTICLE_T] = t;
        WDB[newp + PARTICLE_E] = etr;
        WDB[newp + PARTICLE_U] = u;
        WDB[newp + PARTICLE_V] = v;
        WDB[newp + PARTICLE_W] = w;
        WDB[newp + PARTICLE_LE] = letr;   /* log energy */
        WDB[newp + PARTICLE_PTR_MAT] = (double)mat;

        /* Set array index used in photoelectric.c */
        WDB[newp + PARTICLE_PE_AR_IDX] = RDB[(long)RDB[ptridxrad + ss] + tr];

        /* Check the index */
        CheckValue(FUNCTION_NAME, "PE_AR_IDX", "",
                   RDB[newp + PARTICLE_PE_AR_IDX], 0, PHOTON_NRADTR_MAX);

        /* Set photon emission type */
        WDB[newp + PARTICLE_PHOTON_TYPE] = PHOTON_TYPE_FLUOR;

        /* Score particle balance (JLe 27.1.2017 / 2.1.28) */

        ptr = (long)RDB[RES_G_BALA_SRC];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
        AddBuf(1.0, 1.0, ptr, id, -1, BALA_G_SRC_FLUOR, 0);
        AddBuf(wgt, 1.0, ptr, id, -1, BALA_G_SRC_FLUOR, 1);
        AddBuf(wgt*etr, 1.0, ptr, id, -1, BALA_G_SRC_FLUOR, 2);

        /* Put photon in queue */
        ToQue(newp, id);

        /* Increase transported energy */
        Etrtot += etr;

      }
    }
    else {
      /* Non-radiative transition */

      if (ssk >= 0) {
        /* Increase tertiary subshell vacancy */
        vacancy[ssk]++;

        /* Put tertiary subshell to queue */
        ssqueue[++ssqueueback] = ssk;

        CheckValue(FUNCTION_NAME, "ssqueueback", ": Maximum number of "
                   "vacancies exceeded", (double)ssqueueback,
                   (double)ssqueuefront, PHOTON_ZMAX-1);
      }

      if ((etr > Emin) && useTTB) {

        /* Sample the number of TTB photons */
        TTByield(mat, etr, 1, &letr, MT_ELECTRON_AUGER, &nbr, &idx, id);

        if (nbr > 0) {
          /* Sample electron direction isotropically */
          IsotropicDirection(&u, &v, &w, id);

          /* Use TTB-approximation for the electron, get deposited energy */
          TTBenergy(mat, part, etr, x, y, z, u, v, w, wgt, t,
                    MT_ELECTRON_AUGER, idx, nbr, letr, &EdepTTB, id);
        }
        else {
          /* The electron energy is deposited locally */
          EdepTTB = etr;
        }
      }
      else {
        /* The electron energy is deposited locally */
        EdepTTB = etr;
      }

      CheckValue(FUNCTION_NAME, "EdepTTB", "", EdepTTB, -INFTY, etr);

      /* Increase transported energy */
      Etrtot += etr - EdepTTB;

    }
  }

  /* Set deposited energy.
   * NOTE: Deposited energy can be negative when the vacancy is created by
   * Compton scattering (see comptonscattering.c). It is also possible (but
   * rare) that numerical accuracy in the ENDF data causes a negative
   * energy. */
  *Edep = ebi0 - Etrtot;

  /* NOTE: Is the lower limit adequate? */
  CheckValue(FUNCTION_NAME, "Edep", "", *Edep, -1.0e-3, ebi0);

}

/*****************************************************************************/
