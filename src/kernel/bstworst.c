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
     int next;
     int *list;
} bestworst_data;

static population *selectpop = NULL;

/* select_bestworst()
 *
 * do the actual selection for both the best and worst methods.  both these
 * just create a sorted list of individuals to return, so this function just
 * returns the next in the list.
 */

int select_bestworst ( sel_context *sc )
{
     bestworst_data *bwd;
     bwd = (bestworst_data *)(sc->data);
     return bwd->list[bwd->next++];
}

/* select_best_context()
 *
 * Sets up the best selection method.
 */

sel_context *select_best_context ( int op, sel_context *sc,
                                  population *p, char *string )
{
     int i, j;
     bestworst_data *bwd;
     
     switch ( op )
     {
        case SELECT_INIT:
          sc = (sel_context *)MALLOC ( sizeof ( sel_context ) );
          sc->p = p;
          sc->select_method = select_bestworst;
          sc->context_method = select_best_context;

	  /** the method-specific part is a sorted list of individuals
	    (from best to worst).  the selection function just returns
	    the next element off the list.  the list is built here. **/
	  
          bwd = (bestworst_data *)MALLOC ( sizeof ( bestworst_data ) );
          bwd->list = (int *)MALLOC ( (p->size+1) * sizeof ( int ) );
          bwd->next = 0;
          
          j = 0;
          for ( i = 0; i < p->size; ++i )
               bwd->list[j++] = i;

          selectpop = p;
          qsort ( bwd->list, p->size, sizeof(int), select_best_compare );
          selectpop = NULL;

          sc->data = (void *)bwd;
          return sc;
          break;

        case SELECT_CLEAN:
          bwd = (bestworst_data *)(sc->data);
          FREE ( bwd->list );
          
          FREE ( sc->data );
          FREE ( sc );
          return NULL;
          break;
     }

     return NULL;
}
 
/* select_best_compare()
 *
 * comparison function for qsort() to sort a list of indices in order
 * of fitness (best first).
 */

int select_best_compare ( const void *a, const void *b )
{
     if ( selectpop->ind[*(int *)a].a_fitness ==
          selectpop->ind[*(int *)b].a_fitness )
          return 0;
     else if ( selectpop->ind[*(int *)a].a_fitness <
          selectpop->ind[*(int *)b].a_fitness )
          return 1;
     else 
          return -1;
}

/* select_worst_context()
 *
 * sets up the worst selection method.  identical to select_best_context(),
 * except for the call to qsort().
 */

sel_context *select_worst_context ( int op, sel_context *sc, population *p,
                                   char *string )
{
     int i, j;
     bestworst_data *bwd;
     
     switch ( op )
     {
        case SELECT_INIT:
          sc = (sel_context *)MALLOC ( sizeof ( sel_context ) );
          sc->p = p;
          sc->select_method = select_bestworst;
          sc->context_method = select_worst_context;
          
          bwd = (bestworst_data *)MALLOC ( sizeof ( bestworst_data ) );
          bwd->list = (int *)MALLOC ( (p->size+1) * sizeof ( int ) );
          bwd->next = 0;
          
          j = 0;
          for ( i = 0; i < p->size; ++i )
               bwd->list[j++] = i;
          selectpop = p;
          qsort ( bwd->list, p->size, sizeof(int), select_worst_compare );
          selectpop = NULL;

          sc->data = (void *)bwd;
          return sc;
          break;

        case SELECT_CLEAN:
          bwd = (bestworst_data *)(sc->data);
          FREE ( bwd->list );
          
          FREE ( sc->data );
          FREE ( sc );
          return NULL;
          break;
     }

     return NULL;
}

/* select_worst_compare()
 *
 * comparison function for qsort() to sort a list of indices by fitness,
 * worst first.
 */

int select_worst_compare ( const void *a, const void *b )
{
     if ( selectpop->ind[*(int *)a].a_fitness ==
          selectpop->ind[*(int *)b].a_fitness )
          return 0;
     else if ( selectpop->ind[*(int *)a].a_fitness <
          selectpop->ind[*(int *)b].a_fitness )
          return -1;
     else
          return 1;
}
          
/* select_random_context()
 *
 * sets up the random selection method.
 */

sel_context *select_random_context ( int op, sel_context *sc,
                                    population *p, char *string )
{
     switch ( op )
     {
        case SELECT_INIT:
          sc = (sel_context *)MALLOC ( sizeof ( sel_context ) );
          sc->p = p;
          sc->data = NULL;
          sc->select_method = select_random;
          sc->context_method = select_random_context;
          return sc;
          break;

        case SELECT_CLEAN:
          FREE ( sc );
          return NULL;
          break;
     }

     return NULL;
}

/* select_random()
 *
 * picks an individual at random a returns it.
 */

int select_random ( sel_context *sc )
{
     return random_int ( sc->p->size );
}

