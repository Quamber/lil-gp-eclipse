/*  lil-gp Genetic Programming System, version 1.0, 11 July 1995
 *  Copyright (C) 1995  Michigan State University
 * 
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of version 2 of the GNU General Public License as
 *  published by the Free Software Foundation.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *  
 *  Douglas Zongker       (zongker@isl.cps.msu.edu)
 *  Dr. Bill Punch        (punch@isl.cps.msu.edu)
 *
 *  Computer Science Department
 *  A-714 Wells Hall
 *  Michigan State University
 *  East Lansing, Michigan  48824
 *  USA
 *  
 */

#include "lilgp.h"
/*
 * adapted from the RAN3 routine (in Fortran - ugh) in "Numerical Recipes:
 * The Art of Scientific Computing", page 199.  They adapted it from
 * (surprise) Knuth, Seminumerical Algorithms.
 */

static double mbig = 10000000.0; /* this can be changed to any "big" number. */ 
static double mseed = 1618033.0; /* this can be changed to any smaller but still
				    "big" number. */
static double mz = 0.0;
static double ma[55];       /* the number 55 is special -- see Knuth. */
static int inext, inextp;

/* random_seed()
 *
 * seeds the random number generator using the given int.
 */

void random_seed ( long seed )
{
     int i, i1, k;
     double mj, mk;
     
     mj = mseed - seed;
     mj = mj-((int)(mj/mbig)*mbig);
     ma[54] = mj;
     mk = 1;
     
     for ( i = 1; i < 55; ++i )
     {
          i1 = ((21*i) % 55)-1;
          ma[i1] = mk;
          mk = mj - mk;
          if ( mk < mz )
               mk += mbig;
          mj = ma[i1];
     }
     
     for ( k = 0; k < 4; ++k )
          for ( i = 0; i < 55; ++i )
          {
               ma[i] = ma[i] - ma[(i+30)%55];
               if ( ma[i] < mz )
                    ma[i] += mbig;
          }
     
     inext = 0;
     inextp = 31;     /* the number 31 is special -- see Knuth. */
     
}
     
/* random_int()
 *
 * returns an integer randomly selected from the uniform distribution
 * over the interval [0,max).
 */

int random_int ( int max )
{
     double v = random_double();
     return (int)(v*(double)max);
}

/* random_double()
 *
 * returns a double randomly selected from the uniform distribution
 * over the interval [0,1).
 */

double random_double ( void )
{

     double mj;

     inext = (inext+1)%55;
     inextp = (inextp+1)%55;
     
     mj = ma[inext] - ma[inextp];
     if ( mj < mz )
          mj = mj + mbig;
     ma[inext] = mj;

     return mj/mbig;
}

/* random_get_state()
 *
 * allocates a memory block, saves the state of the random number
 * generator in it, and returns the address.  puts the number of
 * bytes in the block into *size.
 */

void *random_get_state ( int *size )
{
     unsigned char *buffer;
     double *db;
     int i;

     *size = sizeof(double)*58+2*sizeof(int);
     
     buffer = (unsigned char *)MALLOC ( *size );
     db = (double *)buffer;

     db[0] = mbig;
     db[1] = mseed;
     db[2] = mz;
     for ( i = 0; i < 55; ++i )
          db[i+3] = ma[i];
     ((int *)(buffer+58*sizeof(double)))[0] = inext;
     ((int *)(buffer+58*sizeof(double)))[1] = inextp;

#ifdef DEBUG
     fprintf ( stderr, "writing random state: %lf, %lf, %lf, %d, %d\n",
             mbig, mseed, mz, inext, inextp );
#endif

     return buffer;
     
}

/* random_set_state()
 *
 * restores the random number generator state using a block of
 * data previously returned by random_get_state().
 */

void random_set_state ( void *buffer )
{
     unsigned char *cb;
     double *db;
     int i;

     cb = (unsigned char *)buffer;
     db = (double *)buffer;

     mbig = db[0];
     mseed = db[1];
     mz = db[2];
     for ( i = 0; i < 55; ++i )
          ma[i] = db[i+3];
     inext = ((int *)(cb+58*sizeof(double)))[0];
     inextp = ((int *)(cb+58*sizeof(double)))[1];

#ifdef DEBUG
     fprintf ( stderr, "reading random state: %lf, %lf, %lf, %d, %d\n",
             mbig, mseed, mz, inext, inextp );
#endif
     
}



