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

#include <math.h>
#include <stdio.h>

#include "kernel/lilgp.h"

DATATYPE f_multiply ( int tree, farg *args )
{
     return args[0].d * args[1].d;
}

#ifdef TOLERANCE_ZERO

/* fixes the DEC Alpha "high performance arithmetic trap" error.
 * provided by Glen Ropella.
 */

DATATYPE f_protdivide ( int tree, farg *args )
{
     if (args[1].d >= -0.00000000001 || args[1].d <= 0.00000000001)
          return 1.0;
     else
          return args[0].d / args[1].d;
}

#else

DATATYPE f_protdivide ( int tree, farg *args )
{
     if ( args[1].d == 0.0 )
          return 1.0;
     else
          return args[0].d / args[1].d;
}

#endif

DATATYPE f_add ( int tree, farg *args )
{
     return args[0].d + args[1].d;
}

DATATYPE f_subtract ( int tree, farg *args )
{
     return args[0].d - args[1].d;
}

DATATYPE f_sin ( int tree, farg *args )
{
     return sin(args[0].d);
}

DATATYPE f_cos ( int tree, farg *args )
{
     return cos(args[0].d);
}

DATATYPE f_exp ( int tree, farg *args )
{
     return exp(args[0].d);
}

DATATYPE f_rlog ( int tree, farg *args )
{
     if ( args[0].d == 0.0 )
          return 0.0;
     else
          return log ( fabs ( args[0].d ) );
}

DATATYPE f_indepvar ( int tree, farg *args )
{
     return g.x;
}

void f_erc_gen ( DATATYPE *r )
{
     *r = (random_double()*2.0) - 1.0;
}

char *f_erc_print ( DATATYPE d )
{
     static char buffer[20];

     sprintf ( buffer, "%.5f", d );
     return buffer;
}
