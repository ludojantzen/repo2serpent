/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : updatefinixpower.c                             */
/*                                                                           */
/* Created:       2013/03/27 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Updates the power density for Finix nodes                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"

#define FUNCTION_NAME "UpdateFinixPower:"

/*****************************************************************************/
/* Calculates power density at axial region zi, at radius r*/
double GetPden(double r, long buf, long zi, double h, long nr, double rf0, double rf1, long ptr)
{
  long ri=0;
  long nz, idx;
  double val;

  /* Get correct radial bin index */

  while(r*r >= ri/(double)nr*(rf1*rf1-rf0*rf0)+rf0*rf0)
    {
      ri++;
    }

  ri--;

  if(ri == nr)
    ri--;

  /* Get number of axial, angular and time bins */

  nz = (long)RDB[ptr + FUEP_NZ];

  /* Calculate global index */

  idx = zi + nz*0 + nz*1*ri;

  /* Convert power to power density */

  val = RDB[buf + idx]/(h*PI*((rf1*rf1-rf0*rf0)/(double)nr));

  /* Return power density */

  return val;
}

void UpdateFinixPower(long fib, long loc0)
{
  long ptr, i, k, n, m, limptr;
  long nr, nz, pbuf2;
  double r0, r1, sum, dt;
  double rb, rf0, rf1, A, B, C, r_;
  double zmin, zmax, h, tming, tmaxg;
  Rod *rod;
  Boundary_conditions *bc;
  Results *results;
  Options *options;

  /* Check if found */

  if(loc0 < VALID_PTR)
    Die(FUNCTION_NAME, "Could not find universe");

  /* Pointer to power buffer */

  /* These are private to MPI tasks (haven't been broadcast anywhere) */
  /* Should be identical though */

  /* Created in AllocInterfaceStat() */

  pbuf2 = (long)RDB[loc0 + IFC_FUEP_PTR_POWER_REL];
  CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, pbuf2);

  /* Created in ReadIFCFB */

  ptr = (long)RDB[loc0 + IFC_FUEP_OUT_PTR_LIM];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  limptr = ptr;

  /* Get number of bins in different dimensions */

  nr = (long)RDB[ptr + FUEP_NR];
  nz = (long)RDB[ptr + FUEP_NZ];

  /* Get transport time bin limits */

  tmaxg = RDB[DATA_TIME_CUT_TMAX];

  tming = RDB[DATA_TIME_CUT_TMIN];

  /* Calculate time interval length (needed for converting energy to power) */

  if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    dt = 1.0;
  else
    dt = tmaxg - tming;

  /* Get FINIX pointers */

  rod = (Rod *)((long)RDB[fib + FINIX_PTR_ROD]);
  options = (Options *)((long)RDB[fib + FINIX_PTR_OPTIONS]);
  results = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);
  bc = (Boundary_conditions *)((long)RDB[fib + FINIX_PTR_BC]);

  /* Loop over axial regions */
  for (i = 0 ; i < options->axial_nodes ; i++)
    {
      /* Get rmin and rmax (not squared)*/

      rf0 = rod->pellet_inner_radius*100.0;
      rf1 = rod->pellet_outer_radius*100.0;

      /* Get number of radial power nodes */

      nr = (long)RDB[ptr + FUEP_NR];

      /* This is the radial index of the current finix node */
      n = 0;

      /* Get the z-limits of the current tally */

      zmin = RDB[ptr + FUEP_ZMIN];
      zmax = RDB[ptr + FUEP_ZMAX];

      /* Check 2D calculation */

      if ((zmin < -INFTY / 100.0) || (zmax > INFTY / 100.0))
        h = 1.0;
      else
        h = (zmax-zmin)/RDB[ptr + FUEP_NZ];

      /* Calculate total power in this axial zone for LHR */

      /* Reset sum*/

      sum = 0.0;

      /* Loop over radial bins and add to sum */

      for (k = 0 ; k < nr; k++)
        sum += RDB[pbuf2 + i + 0*nz + k*1*nz];

      /* Zero LHR can lead to problems, change it to something small */

      if (sum == 0)
        sum = h*1E-6;

      /* Calculate and store linear heat rate */

      /* Relaxation is handled on the power tally in */
      /* relaxinterfacepower.c */

      bc->linear_power[i] = sum/h*100/dt;

      fprintf(outp, "LHR = %E\n", sum/h/dt);

      while (n < options->pellet_radial_nodes)
        {

          /* First finix node is a special case */

          if (n == 0)
            {

              /* Get inner and outer radius */

              r0 = results->radial_node_position_cold[i][n]*100;
              r1 = results->radial_node_position_cold[i][n+1]*100;

              /* Get power densities at inner and outer radii */

              A = GetPden(r0, pbuf2, i, h, nr, rf0, rf1, limptr);
              B = GetPden(r1, pbuf2, i, h, nr, rf0, rf1, limptr);

              /* This is the pden in the middle of the nodes (node boundary) */
              C = GetPden((r0+r1)/2.0, pbuf2, i, h, nr, rf0, rf1, limptr);

              /* If there is a power tally boundary in the first "control volume"*/
              if (A != C)
                {
                  /* Find the power tally boundary */
                  m = 0;
                  while (r0 >= (rb = sqrt(m/(double)nr*(rf1*rf1-rf0*rf0)+rf0*rf0)))
                    m++;

                  /* If there are two or more boundaries between these nodes things get hard */

                  if (C != B)
                    Die(FUNCTION_NAME, "Please increase the number of Finix radial nodes");

                  /* Calculate power density by a very special average */
                  /* To conserve the power both in axial slice and in this node volume */

                  C = (A*(rb*rb-r0*r0) +
                       B*(-2.0/3.0*(r1*r1+r1*r0+r0*r0)+r0*(r1+r0)+r1*r1-rb*rb))
                    /(-2.0/3.0*(r1*r1+r1*r0+r0*r0)+r0*(r1+r0)+r1*r1-r0*r0);
                }

            }
          else if (n == options->pellet_radial_nodes - 1)
            {

              /* The last finix node is also a special case */

              r0 = results->radial_node_position_cold[i][n]*100;
              r_ = results->radial_node_position_cold[i][n-1]*100;

              /* Get power densities at inner and outer radii */

              A = GetPden(r_, pbuf2, i, h, nr, rf0, rf1, limptr);
              B = GetPden(r0, pbuf2, i, h, nr, rf0, rf1, limptr);

              /* This is the pden in the middle of the nodes */
              /* +1E-6 because the intervals are (r_,r1]*/
              C = GetPden((r_+r0)/2.0+1E-6, pbuf2, i, h, nr, rf0, rf1, limptr);

              /* If there is a power tally boundary in the last "control volume"*/
              if (B != C)
                {
                  /* Find the power tally boundary */
                  m = 0;
                  while (r_ >=  (rb = sqrt(m/(double)nr*(rf1*rf1-rf0*rf0)+rf0*rf0)))
                    m++;


                  /* If there are two or more boundaries between these nodes things get hard */
                  if(C != A)
                    Die(FUNCTION_NAME, "Please increase the number of Finix radial nodes");

                  /* Calculate power density by a very special average */
                  /* To conserve the power both in axial slice and in this node volume */
                  C = (A*(2.0/3.0*(r0*r0+r0*r_+r_*r_)-r_*(r0+r_)+rb*rb-r0*r0)+B*(r0*r0-rb*rb))/
                    (2.0/3.0*(r0*r0+r0*r_+r_*r_)-r_*(r0+r_));
                }

            }
          else
            {

              /* Get radii of adjacent nodes */

              r_ = results->radial_node_position_cold[i][n-1]*100;
              r0 = results->radial_node_position_cold[i][n]*100;
              r1 = results->radial_node_position_cold[i][n+1]*100;

              /* Get power densities at inner and outer radii of node volume */
              /* and at node point */
              A = GetPden((r_+r0)*0.5+1E-6, pbuf2, i, h, nr, rf0, rf1, limptr);
              B = GetPden((r0+r1)*0.5, pbuf2, i, h, nr, rf0, rf1, limptr);
              C = GetPden(r0, pbuf2, i, h, nr, rf0, rf1, limptr);

              /* If there is a power tally boundary in the "control volume" around this node */

              if (A != B)
                {

                  /* Find boundary*/
                  m = 0;
                  while ((r_+r0)*0.5+1E-6 >=  (rb = sqrt(m/(double)nr*(rf1*rf1-rf0*rf0)+rf0*rf0)))
                    m++;


                  /* If there are two or more boundaries between these nodes things get hard */
                  if ((r0+r1)*0.5 > sqrt((m+1)/(double)nr*(rf1*rf1-rf0*rf0)+rf0*rf0))
                    Die(FUNCTION_NAME, "Please increase the number of Finix radial nodes");

                  /* Get power densities at the adjacent node points */

                  A = GetPden(r_, pbuf2, i, h, nr, rf0, rf1, limptr);
                  B = GetPden(r1, pbuf2, i, h, nr, rf0, rf1, limptr);

                  /* Calculate power density by a very special average */
                  /* To conserve the power both in axial slice and in this node volume */

                  C = (A*(2.0/3.0*(r0*r0+r0*r_+r_*r_)-r_*(r0+r_)+rb*rb-r0*r0)+
                     B*(-2.0/3.0*(r1*r1+r1*r0+r0*r0)+r0*(r1+r0)+r1*r1-rb*rb))/
                    (2.0/3.0*(r0*r0+r0*r_+r_*r_)-r_*(r0+r_)
                     -2.0/3.0*(r1*r1+r1*r0+r0*r0)+r0*(r1+r0)
                     +r1*r1-r0*r0);
                }

            }

          /* Store new power density */

          /* The thing stored is not precisely power density but some kind of */
          /* power factor (see FINIX source for more information) */

          bc->power_distr[i][n] = C*1E6*PI/bc->linear_power[i]*
            results->radial_node_position[i][options->pellet_radial_nodes]*
            results->radial_node_position[i][options->pellet_radial_nodes]/dt;

          /* Move to the next node */

          n++;
        }
    }

  /*printf("\n");*/

}

#endif

/*****************************************************************************/
