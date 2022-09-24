/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : rand64.c                                       */
/*                                                                           */
/* Created:       2011/03/03 (TVi)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Main RNG routine                                             */
/*                                                                           */
/* Comments: This pseudo-random integer generator was introduced in the same */
/*           source [1] as the MCNP5 random number generators 2-4 (rng).     */
/*           The most significant 54 bits of the sequence pass the DIEHARD   */
/*           tests [2]. Lowermost 12 bits fail the tests, but, fortunately,  */
/*           these bit can be easily omitted in the floating point number    */
/*           generation.                                                     */
/*                                                                           */
/*           [1] L'Ecuyer, Tables of Linear Congruantial Generators of       */
/*               Different Sizes and Good Lattice Structure, Mathematics of  */
/*               Computation, vol. 68, 225, pp. 249--260 (1999).             */
/*                                                                           */
/*           [2] G.S. Marsaglia, The DIEHARD Battery of Tests of Randomness, */
/*               http://stat.fsu.edu/pub/diehard                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "Rand64:"

/*****************************************************************************/

double Rand64(unsigned long *seed)
{
  /* The actual rng */
  
  *seed = *seed*2862933555777941757 + 12345;

  /* Conversion to floating point number in interval [0,1) */
  
  return (double)(*seed >> 12)/0x0010000000000000;
}

/*****************************************************************************/

