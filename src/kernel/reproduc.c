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
     char *sname;
     sel_context *sc;
} reproduce_data;

/* operator_reproduce_init()
 *
 * called to parse reproduction options and initialize one record of a
 * breedphase table.
 */

int operator_reproduce_init ( char *options, breedphase *bp )
{
     int errors = 0;
     reproduce_data *rd;
     int i, j;
     char **argv;

     rd = (reproduce_data *)MALLOC ( sizeof ( reproduce_data ) );

     /* fill in the breedphase record. */
     bp->operator = OPERATOR_REPRODUCE;
     bp->data = (void *)rd;
     bp->operator_free = operator_reproduce_free;
     bp->operator_start = operator_reproduce_start;
     bp->operator_end = operator_reproduce_end;
     bp->operator_operate = operator_reproduce;

     rd->sname = NULL;

     j = parse_o_rama ( options, &argv );

     for ( i = 0; i < j; ++i )
     {
	  /* parse "select" option. */
          if ( strcmp ( "select", argv[i] ) == 0 )
          {
               if ( !exists_select_method ( argv[++i] ) )
               {
                    ++errors;
                    error ( E_ERROR, "reproduction: \"%s\" is not a known selection method.",
                           argv[i] );
               }
               FREE ( rd->sname );
               rd->sname = (char *)MALLOC ( (strlen(argv[i])+1) * sizeof ( char ) );
               strcpy ( rd->sname, argv[i] );
          }
          else
          {
               ++errors;
               error ( E_ERROR, "reproduction: unknown option \"%s\".",
                      argv[i] );
          }
     }

     free_o_rama ( j, &argv );
     
     if ( rd->sname == NULL )
     {
          ++errors;
          error ( E_ERROR, "reproduction: no selection method specified." );
     }

#ifdef DEBUG
     if ( !errors )
     {
          printf ( "reproduction options:\n" );
          printf ( "   selection: %s\n", rd->sname==NULL?"NULL":rd->sname );
     }
#endif
     
     return errors;
}         

/* operator_reproduce_free()
 *
 * frees the reproduction-specific data of a breedphase record.
 */

void operator_reproduce_free ( void *data )
{
     reproduce_data * rd;

     rd = (reproduce_data *)data;

     FREE ( rd->sname );
     FREE ( rd );
}

/* operator_reproduce_start()
 *
 * gets the selection context for this phase. 
 */

void operator_reproduce_start ( population *oldpop, void *data )
{
     reproduce_data * rd;
     select_context_func_ptr select_con;

     rd = (reproduce_data *)data;
     
     select_con = get_select_context ( rd->sname );
     rd->sc = select_con ( SELECT_INIT, NULL, oldpop, rd->sname );
}

/* operator_reproduce_end()
 *
 * frees the selection context for this phase.
 */

void operator_reproduce_end ( void *data )
{
     reproduce_data * rd;

     rd = (reproduce_data *)data;
     rd->sc->context_method ( SELECT_CLEAN, rd->sc, NULL, NULL );
}


/* operator_reproduce()
 *
 * does the reproduction operation.
 */

void operator_reproduce ( population *oldpop, population *newpop,
                        void *data )
{
     int j;
     reproduce_data * rd;

     rd = (reproduce_data *)data;

     /* select an individual... */
     j = rd->sc->select_method ( rd->sc );

     /* ...and reproduce it into the new population. */
     duplicate_individual ( (newpop->ind)+newpop->next, (oldpop->ind)+j );
     newpop->ind[newpop->next].flags = FLAG_NONE;
     ++newpop->next;
}

