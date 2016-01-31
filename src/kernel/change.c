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

/* the table of operator names and initialization functions.  extend this
 * table whenever you add a new operator.  {NULL,NULL} marks the end of
 * the table.
 */

operator operator_table[] =
{ { "crossover",      operator_crossover_init },
  { "reproduction",   operator_reproduce_init },
  { "mutation",       operator_mutate_init },
  { NULL, NULL } };

/* change_population()
 *
 * breed the new population.
 */

population *change_population ( population *oldpop, breedphase *bp )
{
     population *newpop;
     int i, j;
     int numphases;
     double totalrate = 0.0;
     double r, r2;
     int prob_oper = atoi ( get_parameter ( "probabilistic_operators" ) );

     /* allocate the new population. */
     newpop = allocate_population ( oldpop->size );

     /* the first element of the breedphase table is a dummy -- its
	operator field stores the number of phases. */
     numphases = bp[0].operator;

     /* call the start method for each phase. */
     for ( i = 1; i <= numphases; ++i )
     {
          totalrate += bp[i].rate;
          if ( bp[i].operator_start )
               bp[i].operator_start ( oldpop, bp[i].data );
     }

     /* now fill the new population. */
     while ( newpop->next < newpop->size )
     {

	  /** select an operator, either stochastically or not depending on
	    the probabilistic_operators parameter. **/
          if ( prob_oper )
               r = totalrate * random_double();
          else
               r = totalrate * ((double)newpop->next/(double)newpop->size);

          r2 = bp[1].rate;
          for ( i = 1; r2 < r; )
               r2 += bp[++i].rate;
#ifdef DEBUG
          fprintf ( stderr, "picked %10.3lf; operator %d\n", r, i );
#endif

	  /* call the phase's method to do the operation. */
          if ( bp[i].operator_operate )
               bp[i].operator_operate ( oldpop, newpop, bp[i].data );
     }

     /* call each phase's method to do cleanup. */
     for ( i = 1; i <= numphases; ++i )
     {
          if ( bp[i].operator_end )
               bp[i].operator_end ( bp[i].data );
     }

     /* mark all the ERCs referenced in the new population. */
     for ( i = 0; i < newpop->size; ++i )
          for ( j = 0; j < tree_count; ++j )
               reference_ephem_constants ( newpop->ind[i].tr[j].data, 1 );

     /* free the old population. */
     free_population ( oldpop );

     return ( newpop );
     
}

/* free_breeding()
 *
 * this frees the breedphase table for each subpopulation.
 */

void free_breeding ( multipop *mpop )
{
     int i;

     for ( i = 0; i < mpop->size; ++i )
     {
          free_one_breeding ( mpop->bpt[i] );
          FREE ( mpop->bpt[i] );
     }
     FREE ( mpop->bpt );
     mpop->bpt = NULL;
}

/* free_one_breeding()
 *
 * this frees the breedphase table for a single subpopulation.
 */

void free_one_breeding ( breedphase *bp )
{
     int i;

     for ( i = 1; i <= bp[0].operator; ++i )
     {
          if ( bp[i].operator_free )
               bp[i].operator_free ( bp[i].data );
     }
}

/* initialize_breeding()
 *
 * this builds the breedphase table for each subpopulation.
 */

void initialize_breeding ( multipop *mpop )
{
     char pnamebuf[100];
     int i;

     mpop->bpt = (breedphase **)MALLOC ( mpop->size * sizeof ( breedphase * ) );
     
     for ( i = 0; i < mpop->size; ++i )
     {
          sprintf ( pnamebuf, "subpop[%d].", i+1 );
          mpop->bpt[i] = initialize_one_breeding ( pnamebuf );
     }

}

/* initialize_one_breeding()
 *
 * this builds the breedphase table for a single subpopulation. */

breedphase * initialize_one_breeding ( char *prefix )
{
     char pnamebuf[100];
     char *param, *param2;
     int i, j;
     double rate;
     int errors = 0;
     breedphase *bp;
     operator *op;
     char *name, *namep;

     /* get the number of phases. */
     param = get_breed_parameter ( prefix, "breed_phases" );
     if ( param == NULL )
          error ( E_FATAL_ERROR, "no value specified for \"%sbreed_phases\".",
                 prefix );
     j = atoi ( param );
     if ( j <= 0 )
          error ( E_FATAL_ERROR,
                 "\"%sbreed_phases\" must be greater than zero.", prefix );

     /* the first record of the table is a dummy -- its operator field
	contains the number of phases. */
     bp = (breedphase *)MALLOC ( (j+1) * sizeof ( breedphase ) );
     bp[0].operator = j;
     bp[0].rate = 0.0;
     bp[0].operator_start = NULL;
     bp[0].operator_end = NULL;
     bp[0].operator_free = NULL;
     bp[0].operator_operate = NULL;

     /* for each phase... */
     for ( i = 0; i < j; ++i )
     {
          bp[i+1].operator_start = NULL;
          bp[i+1].operator_end = NULL;
          bp[i+1].operator_free = NULL;
          bp[i+1].operator_operate = NULL;

	  /* get the operator string (name and options) */
          param = get_breed_parameter ( prefix, "breed[%d].operator", i+1 );
          if ( param == NULL )
          {
               ++errors;
               error ( E_ERROR,
                      "no value specifed for \"%sbreed[%d].operator\".",
                      prefix, i+1 );
               continue;
          }

	  /* isolate the name portion of the string. */
          name = (char *)MALLOC ( ( strlen ( param ) + 1 ) * sizeof ( char ) );
          namep = name;
          for ( param2 = param, namep = name;
               *param2 && *param2 != ',' && *param2 != '\n';
               ++param2 )
               if ( !isspace(*param2) )
                    *(namep++) = *param2;
          *namep = 0;
          if ( ! *param2 )
               --param2;

          /* look up the name in the table of operators. */
          op = operator_table;
          while ( op->name )
          {
               if ( strcmp ( op->name, name ) == 0 )
                    break;
               ++op;
          }

          FREE ( name );

          if ( op->name )
	       /* call the operator's initialization function to fill in the fields
		  of the table for this phase. */
               errors += op->func ( param2+1, bp+i+1 );
          else
	       /* the specified operator is not in the table. */
               error ( E_FATAL_ERROR,
                      "%s: \"%s\" is not a known operator.", pnamebuf, param );

	  /* get the rate for this phase. */
          param = get_breed_parameter ( prefix, "breed[%d].rate", i+1 );
          if ( param == NULL )
          {
               ++errors;
               error ( E_ERROR,
                      "no value specified for \"%sbreed[%d].rate\".",
                      prefix, i+1 );
          }
          else
          {
               rate = strtod ( param, NULL );
               if ( rate < 0.0 )
               {
                    ++errors;
                    error ( E_FATAL_ERROR,
                           "\"%sbreed[%d].rate\" must be nonnegative.",
                           prefix, i+1 );
               }
               else if ( rate == 0.0 )
               {
                    error ( E_WARNING,
                           "\"%sbreed[%d].rate\" is zero; is this correct?",
                           prefix, i+1 );
               }
               bp[i+1].rate = rate;
          }
     }

     /* if any errors occurred, stop now. */
     if ( errors )
          error ( E_FATAL_ERROR, "Errors have occurred, aborting." );
     
     return bp;
}

/* rebuild_breeding()
 *
 * rebuilds the breeding table from the parameter database.  called from
 * user code when the breeding parameters change mid-run.
 */

void rebuild_breeding ( multipop *mpop )
{
     free_breeding ( mpop );
     initialize_breeding ( mpop );
}

/* get_breed_parameter()
 *
 * format and following arguments are passed to sprintf to form a string.
 * looks for a parameter called "<prefix><string>", and returns its value.
 * if it does not exist, returns the value of a parameter called "<string>".
 */

char *get_breed_parameter ( char *prefix, char *format, ... )
{
     char *param;
     static char pnamebuf[200];
     int len = strlen(prefix);
     va_list ap;

     strcpy ( pnamebuf, prefix );
     va_start ( ap, format );
     vsprintf ( pnamebuf+len, format, ap );
     va_end ( ap );
     
     param = get_parameter ( pnamebuf );
     if ( param == NULL )
          return get_parameter ( pnamebuf+len );
     else
          return param;
}
