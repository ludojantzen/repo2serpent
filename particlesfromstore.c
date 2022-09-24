/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : particlesfromstore.c                           */
/*                                                                           */
/* Created:       2016/06/07 (VVa)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Retrieves particles from store and puts them into a bank or  */
/*              to source                                                    */
/* Comments:   -Add some neutron/gamma check?                                */
/*             -Currently works for neutrons                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ParticlesFromStore:"

/*****************************************************************************/

void ParticlesFromStore()
{
  long ptr, pos, id, nt, nm, ng, nbin, j, k, tb;
  long nb, i, new, stp, loc0, npart;
  long tmplong, flag;
  char tmpstr[MAX_STR];
  double val;
  FILE *fp;

  /* Stores are only used in coupled transient calculations */

  if (RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_DYN)
    return;

  /*********************/
  /* File based stores */
  /*********************/

#ifdef DNPRINT
  fprintf(outp, "particlesfromstore.c -->\n");
#endif

  /* Create store file name*/

  flag = (long)RDB[DATA_BOI_STORE_NAME];
  CheckValue(FUNCTION_NAME, "flag", "", flag, 0, 1);

  if (flag == 0)
    sprintf(tmpstr, "%s.storeA%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  else
    sprintf(tmpstr, "%s.storeB%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);

  /* Open file for reading */

  fp = fopen(tmpstr, "r");

  /* Get current batch number */

  nb = (long)RDB[DATA_CYCLE_IDX];

  /* Seek to correct position */

  if (nb > 0)
    {
      /* Get correct position (offset) */

      pos = (long)RDB[DATA_NEUTRON_STORE_POS];
      CheckValue(FUNCTION_NAME, "(pos)", "", pos, 0, INFTY);

      /* Seek to correct position */

      if (fseek(fp, pos, SEEK_SET))
        Die(FUNCTION_NAME, "Could not seek in precursor source file");

    }

  /* Read batch number */

  if (fread(&tmplong, sizeof(long), 1, fp) == 0)
    Die(FUNCTION_NAME, "Could not read batch number from file %s", tmpstr);

  /* Check batch number */

  if (tmplong != nb)
    Die(FUNCTION_NAME, "Wrong batch number read from file: %ld (should be %ld)",
        tmplong, nb);

  /* Read number of particles */

  if (fread(&npart, sizeof(long), 1, fp) == 0)
    Die(FUNCTION_NAME, "Could not read batch number from file %s", tmpstr);

#ifdef DNPRINT
  fprintf(outp, "Reading %ld neutrons from store for batch %ld\n", npart, tmplong);
#endif

  /* Reset id */

  id = 0;

  /* Read npart particles */

  for (i = 0; i < npart; i++)
    {

      /* Get template for new neutron */
      /* TODO: Gammas?                */

      new = FromStack(PARTICLE_TYPE_NEUTRON, id);

      /* Read particle from file */

      if (fread(&WDB[new + LIST_DATA_SIZE], sizeof(double),
                (PARTICLE_BLOCK_SIZE - LIST_DATA_SIZE), fp) !=
          (PARTICLE_BLOCK_SIZE - LIST_DATA_SIZE))
        Die(FUNCTION_NAME, "Could not read enough bytes from file %s", tmpstr);

      /* Check type */

      if (RDB[new + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
        Die(FUNCTION_NAME, "Stores not implemented for gammas");

      /* Put particle to bank */

      ToBank(new, id++);

      /* Check id */

      if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
        id = 0;

    }

  /* Get current position in neutron/gamma file */

  pos = ftell(fp);

  /* Store current position in neutron/gamma file */

  WDB[DATA_NEUTRON_STORE_POS] = (double)pos;

  /* Close input file */

  fclose(fp);

  /*******************************/
  /* File based precursor stores */
  /*******************************/

  /* Get pointer to precursor detector or return */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  /* Get current time bin index */
  /* No plus one since reading the BOI values */

  tb = (long)RDB[DATA_DYN_TB];

  /* Get number of time bins */

  nt = (long)RDB[loc0 + PRECDET_NT];

  /* Check current time bin index */

  if (tb >= nt)
    Die(FUNCTION_NAME, "Current time interval %ld, number of time bins %ld", tb, nt);

  /* Get pointer to mesh */

  ptr = (long)RDB[loc0 + PRECDET_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Calculate number of mesh bins */

  nm = 1;
  nm = nm*(long)RDB[ptr + MESH_N0];
  nm = nm*(long)RDB[ptr + MESH_N1];
  nm = nm*(long)RDB[ptr + MESH_N2];

  /* Group bins */

  ng = (long)RDB[loc0 + PRECDET_NG];

  /* Calculate number of bins to save */

  nbin = nm*ng;

  /* Create store file name*/

  flag = (long)RDB[DATA_BOI_STORE_NAME];
  CheckValue(FUNCTION_NAME, "flag", "", flag, 0, 1);

  if (flag == 0)
    sprintf(tmpstr, "%s.storepmeshA%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  else
    sprintf(tmpstr, "%s.storepmeshB%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);

  /* Get index of batch */

  nb = (long)RDB[DATA_CYCLE_IDX];

  /* Open file for reading */

  fp = fopen(tmpstr, "r");

  /* Seek to correct position */

  if (nb > 0)
    {
      /* Get correct position (offset) */

      pos = (long)RDB[DATA_PRECURSOR_MESH_STORE_POS];
      CheckValue(FUNCTION_NAME, "(pos)", "", pos, 0, INFTY);

      /* Seek to correct position */

      if (fseek(fp, pos, SEEK_SET))
        Die(FUNCTION_NAME, "Could not seek in precursor source file");

    }

  /* Read batch number */

  if (fread(&tmplong, sizeof(long), 1, fp) != 1)
    Die(FUNCTION_NAME, "Could not read batch number from file %s", tmpstr);

  /* Check batch number */

  if (nb != tmplong)
    Die(FUNCTION_NAME, "Wrong batch number read from mesh file %ld (expected %ld)",
        tmplong, nb);

  /* Read number of bins */

  if (fread(&tmplong, sizeof(long), 1, fp) != 1)
    Die(FUNCTION_NAME, "Could not read number of bins from file %s", tmpstr);

  /* Check number of bins */

  if (nbin != tmplong)
    Die(FUNCTION_NAME, "Wrong number of bins read from mesh file %ld (expected %ld)",
        tmplong, nbin);

  /* Get pointer to statistics */

  stp = (long)RDB[loc0 + PRECDET_PTR_STAT];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);

  /* Loop over bins to read bin values */

  for (j = 0; j < ng; j++)
    for (k = 0; k < nm; k++)
      {
        /* Read bin value from file */

        if (fread(&val, sizeof(double), 1, fp) != 1)
          Die(FUNCTION_NAME, "Could not read bin value from %s", tmpstr);

        /* Put bin value to the BOI bin of current time interval */

        AddBuf(val, 1.0, stp, 0, -1, tb, j, k);

      }

  /* Get current position in precursor file */

  pos = ftell(fp);

  /* Store current position in precursor file */

  WDB[DATA_PRECURSOR_MESH_STORE_POS] = (double)pos;

  /* Close mesh file */

  fclose(fp);

  /* Return if point-wise precursors don't have to be read */

  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] != PREC_MODE_POINT)
    return;

  /*********************************************************/
  /* Point-wise precursors                                 */

  /* Get pointer to precursor source */

  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check that precursor source is empty */

  if (ListSize(ptr) != 1)
    Die(FUNCTION_NAME, "Source is not empty");

  /* Create store file name*/

  flag = (long)RDB[DATA_BOI_STORE_NAME];
  CheckValue(FUNCTION_NAME, "flag", "", flag, 0, 1);

  if (flag == 0)
    sprintf(tmpstr, "%s.storeprecA%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);
  else
    sprintf(tmpstr, "%s.storeprecB%d", GetText(DATA_PTR_INPUT_FNAME), mpiid);

  /* Open file for reading */

  fp = fopen(tmpstr, "r");

  /* Seek to correct position */

  if (nb > 0)
    {
      /* Get correct position (offset) */

      pos = (long)RDB[DATA_PRECURSOR_STORE_POS];
      CheckValue(FUNCTION_NAME, "(pos)", "", pos, 0, INFTY);

      /* Seek to correct position */

      if (fseek(fp, pos, SEEK_SET))
        Die(FUNCTION_NAME, "Could not seek in precursor source file");

    }

  /* Read batch number */

  if (fread(&tmplong, sizeof(long), 1, fp) == 0)
    Die(FUNCTION_NAME, "Could not read batch number from file %s", tmpstr);

  /* Check batch number */

  if (tmplong != nb)
    Die(FUNCTION_NAME, "Wrong batch number read from file: %ld (should be %ld)", tmplong, nb);

  /* Read number of particles */

  if (fread(&npart, sizeof(long), 1, fp) == 0)
    Die(FUNCTION_NAME, "Could not read batch number from file %s", tmpstr);

#ifdef DNPRINT
  fprintf(outp, "Reading %ld precursors from store for batch %ld\n", npart, tmplong);
#endif

  /* Reset id */

  id = 0;

  /* Read npart particles */

  for (i = 0; i < npart; i++)
    {

      /* Get template for new precursor */

      new = FromStack(PARTICLE_TYPE_PRECURSOR, id++);

      /* Check id */

      if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
        id = 0;

      /* Read particle from file */

      if (fread(&WDB[new + LIST_DATA_SIZE], sizeof(double), (PARTICLE_BLOCK_SIZE - LIST_DATA_SIZE), fp) !=
          (PARTICLE_BLOCK_SIZE - LIST_DATA_SIZE))
        Die(FUNCTION_NAME, "Could not read enough bytes from file %s", tmpstr);


      /* Put precursor to precursor source (will be sorted later) */

      AddItem(DATA_PART_PTR_PSOURCE, new);

    }

  /* Get current position in neutron/gamma file */

  pos = ftell(fp);

  /* Store current position in neutron/gamma file */

  WDB[DATA_PRECURSOR_STORE_POS] = (double)pos;

  /* Close input file */

  fclose(fp);

#ifdef DNPRINT
  fprintf(outp, "<-- particlesfromstore.c\n");
#endif
  /***************************************************************************/
}

/*****************************************************************************/
