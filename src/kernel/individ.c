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

/* print_individual_stdout()
 *
 * prints the given individual to stdout.  used with the "call" command
 * in gdb (stdout is a #defined symbol and so is not available to the
 * debugger, at least under Solaris...)
 */

void print_individual_stdout ( individual *ind )
{
     print_individual ( ind, stdout );
}

/* print_individual()
 *
 * prints the given individual to the given FILE *.
 */

void print_individual ( individual *ind, FILE *f )
{
     int j;

     for ( j = 0; j < tree_count; ++j )
     {
          fprintf ( f, "   %s: ", tree_map[j].name );
          print_tree ( ind->tr[j].data, f );
     }
}

/* pretty_print_individual()
 *
 * pretty-prints the individual to the given FILE *.  shows expression
 * structure via indentation.
 */

void pretty_print_individual ( individual *ind, FILE *f )
{
     int j;

     for ( j = 0; j < tree_count; ++j )
     {
          fprintf ( f, "%s:", tree_map[j].name );
          pretty_print_tree ( ind->tr[j].data, f );
     }
}

/* individual_size()
 *
 * returns the total number of nodes in an individual.
 */

int individual_size ( individual *ind )
{
     int j, k = 0;
     
     for ( j = 0; j < tree_count; ++j )
          k += ind->tr[j].nodes;

     return k;
}

/* individual_depth()
 *
 * returns the depth of an individual (maximum of the depths
 * of its trees).
 */

int individual_depth ( individual *ind )
{
     int i, j, k = 0;

     for ( j = 0; j < tree_count; ++j )
     {
          if ( ( i = tree_depth ( ind->tr[j].data ) ) > k )
               k = i;
     }
     return k;
}
     
/* duplicate_individual()
 *
 * duplicates an individual.
 */

void duplicate_individual ( individual *to, individual *from )
{
     int j;
     for ( j = 0; j < tree_count; ++j )
          copy_tree ( to->tr+j, from->tr+j );
     to->r_fitness = from->r_fitness;
     to->s_fitness = from->s_fitness;
     to->a_fitness = from->a_fitness;
     to->hits = from->hits;
     to->evald = from->evald;
     to->flags = from->flags;
}

