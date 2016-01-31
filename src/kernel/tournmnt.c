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
typedef struct
{
     int count;
} tournament_data;

/* select_tournament_context()
 *
 * returns a selection context for the tournament selection method.
 */

sel_context *select_tournament_context ( int op, sel_context *sc,
                                        population *p, char *string )
{
     char **argv;
     int i, j;
     tournament_data *td;
     
     switch ( op )
     {
        case SELECT_INIT:

          sc = (sel_context *)MALLOC ( sizeof ( sel_context ) );
	  /* fill in fields of the selection context. */
          sc->p = p;
          sc->select_method = select_tournament;
          sc->context_method = select_tournament_context;

	  /* store the tournament data record. */
          td = (tournament_data *)MALLOC ( sizeof ( tournament_data ) );

	  /* parse the options. */
          td->count = 2;
          j = parse_o_rama ( string, &argv );
          for ( i = 1; i < j; ++i )
          {
	       /* "size" is the only valid option. */
               if ( strcmp ( argv[i], "size" ) == 0 )
                    td->count = atoi ( argv[++i] );
               else
                    error ( E_FATAL_ERROR, "unknown tournament option \"%s\".",
                           argv[i] );
          }

	  free_o_rama ( j, &argv );
          
          if ( td->count <= 0 )
               error ( E_FATAL_ERROR,
                      "tournament size must be at least 1.  (%s)", string );

          sc->data = (void *)td;
          return sc;
          break;
          
        case SELECT_CLEAN:

          td = (tournament_data *)(sc->data);
          FREE ( sc->data );
          FREE ( sc );
          return NULL;
          break;
     }

     return NULL;
}

/* select_tournament()
 *
 * does a tournament selection.  randomly picks (uniformly) a number (size)
 * of individuals, then selects the best one among those.
 */

int select_tournament ( sel_context *sc )
{
     int i, j, k;
     tournament_data *td;
     population *p;

     td = (tournament_data *)(sc->data);
     p = sc->p;

     j = -1;
     for ( i = 0; i < td->count; ++i )
     {
	  /* pick another individual. */
          k = random_int ( p->size );
	  /* save it if it is better than the current best. */
          if ( j == -1 || p->ind[k].a_fitness > p->ind[j].a_fitness )
               j = k;
     }
     
     return j;
}


