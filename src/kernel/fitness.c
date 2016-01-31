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
/* select_afit_context()
 *
 * creates context for the fitness selection method.
 */

sel_context *select_afit_context ( int op, sel_context *sc, population *p,
                                  char *string )
{
     interval_data *id;
     int i, j = 0;

     switch ( op )
     {
        case SELECT_INIT:
          sc = (sel_context *)MALLOC ( sizeof ( sel_context ) );
          sc->p = p;
          sc->select_method = select_interval;
          sc->context_method = select_afit_context;

	  /* the interval_data structure (used with select_interval()) is
	     essentially a list of interval widths and indices for each
	     individual. */
	  
          id = (interval_data *)MALLOC ( sizeof ( interval_data ) );

          id->ri = (reverse_index *)MALLOC ( (p->size+1) *
                                           sizeof ( reverse_index ) );
          id->total = 0.0;
          id->count = p->size;
          id->ri[j].fitness = 0.0;
          id->ri[j].index = -1;
          ++j;
          
          for ( i = 0; i < p->size; ++i )
          {
	       /* interval width is the adjusted fitness. */
               id->total += p->ind[i].a_fitness;
               id->ri[j].fitness = id->total;
               id->ri[j].index = i;
               ++j;
          }

          sc->data = (void *)id;
          return sc;
          break;

        case SELECT_CLEAN:
          id = (interval_data *)(sc->data);
          FREE ( id->ri );
          
          FREE ( sc->data );
          FREE ( sc );
          return NULL;
          break;
     }

     return NULL;
}

/* select_inverse_afit_context()
 *
 * creates context for inverse_fitness selection method.
 */

sel_context *select_inverse_afit_context ( int op, sel_context *sc,
                                          population *p, char *string )
{
     interval_data *id;
     int i, j = 0;

     switch ( op )
     {
        case SELECT_INIT:
          sc = (sel_context *)MALLOC ( sizeof ( sel_context ) );
          sc->p = p;
          sc->select_method = select_interval;
          sc->context_method = select_inverse_afit_context;

	  /** use select_interval() to do the selection. **/
	  
          id = (interval_data *)MALLOC ( sizeof ( interval_data ) );

          id->ri = (reverse_index *)MALLOC ( (p->size+1) *
                                           sizeof ( reverse_index ) );
          id->total = 0.0;
          id->count = p->size;
          id->ri[j].fitness = 0.0;
          id->ri[j].index = -1;
          ++j;
          
          for ( i = 0; i < p->size; ++i )
          {
	       /* interval width is inverse of adjusted fitness. */
               id->total += 1.0/p->ind[i].a_fitness;
               id->ri[j].fitness = id->total;
               id->ri[j].index = i;
               ++j;
          }
          
          sc->data = (void *)id;
          return sc;
          break;

        case SELECT_CLEAN:
          id = (interval_data *)(sc->data);
          FREE ( id->ri );
          
          FREE ( sc->data );
          FREE ( sc );
          return NULL;
          break;
     }

     return NULL;
}

/* select_afit_overselect_context()
 *
 * creates context for fitness_overselect selection method.
 */

sel_context *select_afit_overselect_context ( int op, sel_context *sc,
                                             population *p, char *string )
{
     interval_data *id;
     int i, j;
     double total;
     double group1_cutoff = 0.32;
     int cutoffset;
     double group1_selection = 0.8;
     double cutoff;
     double temp;
     char **argv;

     switch ( op )
     {
        case SELECT_INIT:

          sc = (sel_context *)MALLOC ( sizeof ( sel_context ) );
          sc->p = p;
          sc->select_method = select_interval;
          sc->context_method = select_afit_overselect_context;

	  /** parse the options string. **/
	  
          j = parse_o_rama ( string, &argv );

          cutoffset = 0;
          for ( i = 1; i < j; ++i )
          {
               if ( strcmp ( argv[i], "cutoff" ) == 0 )
               {
                    group1_cutoff = strtod ( argv[++i], NULL );
                    cutoffset = 1;
               }
               else if ( strcmp ( argv[i], "proportion" ) == 0 )
                    group1_selection = strtod ( argv[++i], NULL );
               else
                    error ( E_FATAL_ERROR, "unknown fitness_overselect option \"%s\".",
                           argv[i] );
          }

	  /* if the cutoff was not set manually, then set it based on
	     the population size. */
          if ( !cutoffset )
          {
               group1_cutoff = 320.0/p->size;
               if ( group1_cutoff < 0.0 || group1_cutoff > 1.0 )
                    group1_cutoff = 0.32;
          }

	  free_o_rama ( j, &argv );
	  
          if ( group1_cutoff < 0.0 || group1_cutoff > 1.0 )
               error ( E_FATAL_ERROR, "Overselected fitness cutoff out of range.  (%s)", string );

          if ( group1_selection < 0.0 || group1_selection > 1.0 )
               error ( E_FATAL_ERROR, "Overselected fitness proportion out of range.  (%s)", string );
          
          id = (interval_data *)MALLOC ( sizeof ( interval_data ) );

          id->ri = (reverse_index *)MALLOC ( (p->size+1) *
                                            sizeof ( reverse_index ) );
          id->total = 0.0;
          id->count = p->size;

          id->ri[0].fitness = 0.0;
          id->ri[0].index = -1;
          j = 1;
          
          /* store the fitness values in the reverse_index */
          total = 0.0;
          for ( i = 0; i < p->size; ++i )
          {
               total += p->ind[i].a_fitness;
               id->ri[j].fitness = p->ind[i].a_fitness;
               id->ri[j].index = i;
               ++j;
          }

          /* (sort lowest first) */
          qsort ( (id->ri)+1, p->size, sizeof ( reverse_index ),
                 rev_ind_compare );

	  /* find the top individuals accounting for (cutoff) of the fitness,
	     and multiply their interval width by the selection.  multiply
	     all the others by (1-selection). */
          cutoff = total * (1.0-group1_cutoff);
          total = 0.0;
          for ( i = 1; i < p->size+1; ++i )
          {
               if ( total >= cutoff )
                    temp = id->ri[i].fitness * group1_selection;
               else
                    temp = id->ri[i].fitness * (1.0-group1_selection);
               total += id->ri[i].fitness;
               id->ri[i].fitness = temp;
          }

          /* now convert to cumulative. */
          total = 0.0;
          for ( i = 1; i < p->size+1; ++i )
          {
               total += id->ri[i].fitness;
               id->ri[i].fitness = total;
          }

          id->total = total;

          sc->data = (void *)id;
          return sc;
          break;
          
        case SELECT_CLEAN:
          id = (interval_data *)(sc->data);
          FREE ( id->ri );
          
          FREE ( sc->data );
          FREE ( sc );
          return NULL;
          break;
     }

     return NULL;
}

