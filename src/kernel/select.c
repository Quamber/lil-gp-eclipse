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
/* the table listing selection method names and the functions which
   create selection contexts.  extend this table whenever you add a
   new selection method.  {NULL,NULL} marks the end of the table. */

select_method select_method_table[] =
{ { "fitness",            select_afit_context },
  { "fitness_overselect", select_afit_overselect_context },
  { "tournament",         select_tournament_context },
  { "inverse_fitness",    select_inverse_afit_context },
  { "best",               select_best_context },
  { "worst",              select_worst_context },
  { "random",             select_random_context },
  { NULL, NULL } };

/* context_func_ptr()
 *
 * looks for the named string in the selection method table
 * and returns the corresponding context method.
 */

select_context_func_ptr get_select_context ( char *string )
{
     int i, j, k;
     char *name;
     select_method *s = select_method_table;

     /* pull off the name section of the string. */
     for ( i = 0; string[i] != 0 && string[i] != ','
          && string[i] != '\n'; ++i );
     name = (char *)MALLOC ( i+1 );
     k = 0;
     for ( j = 0; j < i; ++j )
          if ( !isspace(string[j]) )
               name[k++] = string[j];
     name[k] = 0;

     /* search the table. */
     while ( s->name != NULL )
     {
          if ( strcmp ( s->name, name ) == 0 )
               break;
          ++s;
     }
     FREE ( name );

     /* return the function (NULL if it wasn't found). */
     return s->func;
}

/* exists_select_method()
 *
 * returns 1 if the named selection method exists, 0 otherwise.
 */

int exists_select_method ( char *string )
{
     return ( get_select_context ( string ) != NULL );
}

/* free_o_rama()
 *
 * frees the argv-style array produced by parse_o_rama().
 */

void free_o_rama ( int j, char ***argv )
{
     int i;

     for ( i = 0; i < j; ++i )
          FREE ( (*argv)[i] );
     FREE ( *argv );
     *argv = NULL;
}

/* parse_o_rama()
 *
 * breaks a string into an argv-style array.  field delimiters
 * are newlines, commas, equal-signs, NULLS, and open-parentheses.
 * no field breaking occurs within a pair of nested parentheses.
 * any whitespace is removed from the string.  any zero-length
 * fields are ignored.
 *
 * currently seg faults parentheses are mismatched.  this should
 * be fixed.
 */

int parse_o_rama ( char *string, char ***argv )
{
     int i, j;
     char *p, *m, *o;
     int parendepth = 0;
     char **fargv;
     int nonblank = 0;

     /** this is not commented because it sorely needs to be
       rewritten. **/
     
     /** pass through the string counting fields **/
     for ( j = 1, p = string; *p; ++p )
     {
          j += (parendepth==0)&&(*p=='\n'||*p=='='||*p==','||*p==0||*p=='(');
          if ( *p == '(' )
               ++parendepth;
          else if ( *p == ')' )
          {
               --parendepth;
               if ( parendepth == 0 )
                    ++j;
               else if ( parendepth < 0 )
                    return -1;
          }
     }

     if ( *string == 0 )
          j = 0;

     /* allocate space for the fields. */
     fargv = (char **)MALLOC ( j * sizeof ( char * ) );
     p = m = string;
     parendepth = 0;
     for ( i = 0; i < j; ++i )
     {
          while ( parendepth || ( *p && *p != '=' && *p != ',' && *p != '\n'
                                 && *p != '(' && *p != ')' ) )
          {
               if ( *p == '(' )
                    ++parendepth;
               else if ( *p == ')' )
               {
                    --parendepth;
                    if ( !parendepth )
                         --p;
               }
               ++p;
          }
          if ( *p == '(' )
               ++parendepth;
          fargv[i] = (char *)MALLOC ( (p-m+1) * sizeof ( char ) );
          o = fargv[i];
          while ( m < p )
          {
               if ( !isspace(*m) )
                    *(o++) = *m;
               ++m;
          }
          *o = 0;
          if ( fargv[i][0] != 0 )
               ++nonblank;
          ++p;
          ++m;
     }

     *argv = (char **)MALLOC ( nonblank * sizeof ( char * ) );
     nonblank = 0;
#ifdef DEBUG
     printf ( "parse_o_rama:  %d fields\n", nonblank );
#endif
     for ( i = 0; i < j; ++i )
     {
          if ( fargv[i][0] )
          {
#ifdef DEBUG
               printf ( "   [%s]\n", fargv[i] );
#endif
               (*argv)[nonblank++] = fargv[i];
          }
          else
               FREE ( fargv[i] );
     }
     FREE ( fargv );
     
     return nonblank;
}

/* rev_ind_compare()
 *
 * comparison function for sorting a reverse_index table by increasing
 * fitness.
 */

int rev_ind_compare ( const void *a, const void *b )
{
     if ( ((reverse_index *)a)->fitness > ((reverse_index *)b)->fitness )
          return 1;
     else if ( ((reverse_index *)a)->fitness < ((reverse_index *)b)->fitness )
          return -1;
     else
          return 0;
}

/* select_interval()
 *
 * for selection methods which can be expressed as randomly selecting an
 * individual, where each individual has some fixed probability of
 * being selected, this efficiently does the selection.
 *
 * the selection_context's data field must point to an interval_data
 * structure, which contains (essentially) a list of consecutive intervals
 * and which individuals they correspond to.  this function chooses a random
 * number in the whole range of the intervals and uses binary search to
 * locate which interval that falls in, returning the corresponding
 * index.
 */

int select_interval ( sel_context *sc )
{
     double rval;
     int middle;
     interval_data *id = sc->data;
     int low = 0, high = id->count;
     
     rval = random_double() * id->total;

#ifdef DEBUG_INTERVAL
     printf ( "random value is %.6f (%.6f)\n", rval, id->total );
#endif
     
     while ( low < high-1 )
     {
#ifdef DEBUG_INTERVAL
          printf ( "current range is %d (%.6f) to %d (%.6f)\n",
                  low, id->ri[low].fitness,
                  high, id->ri[high].fitness );
#endif          
          
          middle = (low+high)/2;
          if ( rval >= id->ri[middle].fitness )
               low = middle;
          else
               high = middle;
     }

#ifdef DEBUG_INTERVAL
     printf ( "selected value %d (%.6f)\n", high,
             id->ri[high].fitness );
#endif

     if ( id->ri[high].index == -1 )
     {
	  /* this shouldn't ever happen either, but I'm nervous about
	     off-by-one errors on binary searches. */
          fprintf ( stderr, "afitness select misfired.\n" );
     }
     
     return id->ri[high].index;
     
}

