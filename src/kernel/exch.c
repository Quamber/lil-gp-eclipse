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
/* initialize_topology()
 *
 * reads the parameter database and builds the exchange table.
 */

void initialize_topology ( multipop *mpop )
{
     char pnamebuf[100], pnamebuf2[100];
     char *param, *param2, *param3;
     int i, j, k;
     int errors = 0;

     if ( mpop->size == 1 )
     {
	  /* singlepop problem -- no topology needed. */
          mpop->exch = NULL;
          mpop->exchanges = -1;
          return;
     }

     oprintf ( OUT_SYS, 30, "building subpopulation exchange topology:\n" );
     
     param = get_parameter ( "multiple.exchanges" );
     if ( param == NULL )
     {
	  /* multipop problem, but no exchanges specified. */
	  
          mpop->exch = NULL;
          mpop->exchanges = 0;
          return;
     }
     else
     {
          mpop->exchanges = atoi ( param );
          if ( mpop->exchanges < 0 )
               error ( E_FATAL_ERROR, "\"exchanges\" must be nonnegative." );
     }

     mpop->exch = (exchange *)MALLOC ( mpop->exchanges * sizeof ( exchange ) );
     
     for ( i = 0; i < mpop->exchanges; ++i )
     {
	  /** read the destination subpop. **/
	  
          sprintf ( pnamebuf, "exch[%d].to", i+1 );
          param = get_parameter ( pnamebuf );
          if ( param == NULL )
          {
               ++errors;
               error ( E_ERROR, "\"%s\" must be set.", pnamebuf );
          }
          else
          {
               mpop->exch[i].to = atoi ( param ) - 1;
               if ( mpop->exch[i].to < 0 || mpop->exch[i].to >= mpop->size )
               {
                    ++errors;
                    error ( E_ERROR, "\"%s\" is out of range.\n", pnamebuf );
               }
          }

	  /** read how the individuals to be replaced in the destination
	    subpop are selected. **/
	  
          sprintf ( pnamebuf, "exch[%d].toselect", i+1 );
          mpop->exch[i].tosc = get_parameter ( pnamebuf );
          if ( mpop->exch[i].tosc == NULL )
          {
               ++errors;
               error ( E_ERROR, "\"%s\" must be set.", pnamebuf );
          }
          else
          {
               if ( ! exists_select_method ( mpop->exch[i].tosc ) )
               {
                    ++errors;
                    error ( E_ERROR, "\"%s\": \"%s\" is not a selection method.",
                           pnamebuf, mpop->exch[i].tosc );
               }
          }

	  /** read how many individuals are to be exchanged in this
	    manner. **/
	  
          sprintf ( pnamebuf, "exch[%d].count", i+1 );
          param = get_parameter ( pnamebuf );
          if ( param == NULL )
          {
               ++errors;
               error ( E_ERROR, "\"%s\" must be set.", pnamebuf );
          }
          else
          {
               mpop->exch[i].count = atoi ( param );
               if ( mpop->exch[i].count < 0 )
               {
                    ++errors;
                    error ( E_ERROR, "\"%s\" must be nonnegative.", pnamebuf );
               }
          }

	  /** check to see if "from" is specified without a "tree[#]". **/
	  
          sprintf ( pnamebuf, "exch[%d].from", i+1 );
          param = get_parameter ( pnamebuf );
          if ( param )
          {
               /** if "from" is specified, then we're copying whole individuals
		 from one subpop to another. **/

	       /* these arrays are not needed. */
               mpop->exch[i].from = NULL;
               mpop->exch[i].as = NULL;
	       /* allocate an array of one string (to hold the selection
		  method). */
               mpop->exch[i].fromsc = (char **)MALLOC ( sizeof ( char * ) );

	       /* the subpop that individuals are taken from. */
               mpop->exch[i].copywhole = atoi ( param ) - 1;
               if ( mpop->exch[i].copywhole < 0 ||
                    mpop->exch[i].copywhole >= mpop->size )
               {
                    ++errors;
                    error ( E_ERROR, "\"%s\" is out of range.", pnamebuf );
               }

	       /* the selection method used to pick the individuals from the
		  source subpop. */
               sprintf ( pnamebuf, "exch[%d].fromselect", i+1 );
               mpop->exch[i].fromsc[0] = get_parameter ( pnamebuf );
               if ( mpop->exch[i].fromsc[0] == NULL )
               {
                    ++errors;
                    error ( E_ERROR, "\"%s\" must be set.", pnamebuf );
               }
               else
               {
                    if ( ! exists_select_method ( mpop->exch[i].fromsc[0] ) )
                    {
                         ++errors;
                         error ( E_ERROR, "\"%s\": \"%s\" is not a selection method.",
                                pnamebuf, mpop->exch[i].fromsc[0] );
                    }
               }
          }
          else
          {

               /** since "from" is not defined, we're taking trees from different
		 subpops and merging them to create a composite individual to place
		 in the destination subpop. **/
	       
               mpop->exch[i].copywhole = -1;
	       /* this array lists, for each tree, which subpop it comes from. */
               mpop->exch[i].from = (int *)MALLOC ( tree_count * sizeof ( int ) );
	       /* this array keeps track of when two trees are supposed to always
		  come from the same individual (not just the same subpop). */
               mpop->exch[i].as = (int *)MALLOC ( tree_count * sizeof ( int ) );
	       /* this array holds the selection method strings used for each
		  tree. */
               mpop->exch[i].fromsc = (char **)MALLOC ( tree_count * sizeof ( char * ) );

	       /* get the default selection method, if one is specified. */
               sprintf ( pnamebuf, "exch[%d].fromselect", i+1 );
               param3 = get_parameter ( pnamebuf );
               
               for ( j = 0; j < tree_count; ++j )
               {
		    /** for each tree, attempt to read the "from" and
		      "fromselect" parameters. **/
		    
                    sprintf ( pnamebuf, "exch[%d].from.tree[%d]", i+1, j );
                    param = get_parameter ( pnamebuf );
                    sprintf ( pnamebuf2, "exch[%d].fromselect.tree[%d]",
                             i+1, j );
                    param2 = get_parameter ( pnamebuf2 );

                    if ( param == NULL && param2 == NULL )
                    {
			 /* neither is set, we're supposed to leave this
			    tree untouched in the destination individual. */
			 
                         mpop->exch[i].from[j] = -1;
                         mpop->exch[i].as[j] = -1;
                         mpop->exch[i].fromsc[j] = NULL;
                    }
                    else if ( param2 == NULL )
                    {
                         /* only "from" is set, examine param3 for default
			    selection method. */

			 /* source subpop. */
                         mpop->exch[i].from[j] = atoi ( param ) - 1;
                         if ( mpop->exch[i].from[j] < 0 || mpop->exch[i].from[j] >= mpop->size )
                         {
                              ++errors;
                              error ( E_ERROR, "\"%s\" is out of range.", pnamebuf );
                         }

			 /* no default set, error. */
                         if ( param3 == NULL )
                         {
                              ++errors;
                              error ( E_ERROR, "\"%s\" must be set.", pnamebuf2 );
                         }
                         else
                         {
                              mpop->exch[i].as[j] = -1;
                              if ( ! exists_select_method ( param3 ) )
                              {
                                   ++errors;
                                   error ( E_ERROR, "\"%s\": \"%s\" is not a selection method.",
                                          pnamebuf, param3 );
                              }
                         }
                         mpop->exch[i].fromsc[j] = param3;
                    }
                    else if ( param == NULL )
                    {
                         /* only "fromselect" is set; it better be of the form
			    "as_#". */
                         
                         if ( strncmp ( param2, "as_", 3 ) == 0 )
                         {
                              mpop->exch[i].from[j] = -1;
                              mpop->exch[i].fromsc[j] = NULL;
			      /* "as" stores which tree this one comes from the
				 same subpop as. */
                              mpop->exch[i].as[j] = atoi ( param2 + 3 );
                              if ( mpop->exch[i].as[j] < 0 ||
                                  mpop->exch[i].as[j] >= tree_count )
                              {
                                   ++errors;
                                   error ( E_ERROR, "\"%s\" is out of range.", pnamebuf2 );
                              }
                         }
                         else
                         {
                              ++errors;
                              error ( E_ERROR, "\"%s\" must be \"as_#\".", pnamebuf2 );
                         }
                    }
                    else
                    {
                         /* they're both set. */

                         mpop->exch[i].as[j] = -1;
                         mpop->exch[i].from[j] = atoi ( param ) - 1;
                         if ( mpop->exch[i].from[j] < 0 || mpop->exch[i].from[j] >= mpop->size )
                         {
                              ++errors;
                              error ( E_ERROR, "\"%s\" is out of range.", pnamebuf );
                         }
                         mpop->exch[i].fromsc[j] = param2;
                         if ( ! exists_select_method ( param2 ) )
                         {
                              ++errors;
                              error ( E_ERROR, "\"%s\": \"%s\" is not a selection method.",
                                     pnamebuf2, param2 );
                         }
                    }
               }

	       /* now we need to resolve any chains of "as_" references: if
		  tree 2 comes from the same individual as tree 1, and tree 1
		  comes from the same individual as tree 0, we need to change that
		  to say that both 2 and 1 come from tree 0.

		  also detect circular references. */
	       
               for ( j = 0; j < tree_count; ++j )
               {
                    if ( mpop->exch[i].as[j] == -1 )
                         continue;
                    k = mpop->exch[i].as[j];
                    while ( k != -1 )
                    {
                         if ( k == j )
                         {
                              ++errors;
                              error ( E_ERROR, "Circular reference resolving \"exch[%d].fromselect.tree[%d]\".",
                                     i+1, j );
                              j = tree_count;
                              break;
                         }
                         mpop->exch[i].as[j] = k;
                         k = mpop->exch[i].as[k];
                    }
                    k = mpop->exch[i].as[j];
                    if ( mpop->exch[i].from[k] == -1 && mpop->exch[i].as[k] == -1 )
                         mpop->exch[i].as[j] = -1;
               }
          }

          
#ifdef DEBUG
	  /* print out information on this exchange. */
          printf ( "exchange %d:\n", i+1 );
          printf ( "to: %d; count: %d; select: %s\n", mpop->exch[i].to,
                  mpop->exch[i].count, mpop->exch[i].tosc );
          if ( mpop->exch[i].copywhole == -1 )
          {
               for ( j = 0; j < tree_count; ++j )
               {
                    param = mpop->exch[i].fromsc[j];
                    printf ( "   %3d:  from: %3d   as: %3d   select: %s\n",
                            j, mpop->exch[i].from[j], mpop->exch[i].as[j],
                            param==NULL?"NULL":param );
               }
          }
          else
          {
               param = mpop->exch[i].fromsc[0];
               printf ( "copywhole: %d   select: %s\n",
                       mpop->exch[i].copywhole, param==NULL?"NULL":param );
          }
#endif
     }

     /* if any errors occurred then stop now. */
     if ( errors )
          error ( E_FATAL_ERROR, "Errors occurred while building topology.  Aborting." );

     /* print out the summary of exchanges. */
     oprintf ( OUT_SYS, 30, "    %d exchange(s) total.\n", mpop->exchanges );
     for ( i = 0; i < mpop->exchanges; ++i )
     {
          oprintf ( OUT_SYS, 30, "    exchange %d:\n", i+1 );
          oprintf ( OUT_SYS, 30, "        replace %d individual(s) in subpop %d (selected by %s)\n",
                   mpop->exch[i].count, mpop->exch[i].to+1, mpop->exch[i].tosc );
          if ( mpop->exch[i].copywhole != -1 )
               oprintf ( OUT_SYS, 30, "        with individual(s) from subpop %d (selected by %s)\n",
                        mpop->exch[i].copywhole+1, mpop->exch[i].fromsc[0] );
          else
               for ( j = 0; j < tree_count; ++j )
               {
                    if ( mpop->exch[i].from[j] == -1 )
                    {
                         if ( mpop->exch[i].as[j] == -1 )
                              oprintf ( OUT_SYS, 30, "        tree %d: leaving original tree\n", j );
                         else
                              oprintf ( OUT_SYS, 30, "        tree %d: from same individual as tree %d\n", j, mpop->exch[i].as[j] );
                    }
                    else
                         oprintf ( OUT_SYS, 30, "        tree %d: from subpop %d (selected by %s)\n", j,
                                  mpop->exch[i].from[j]+1, mpop->exch[i].fromsc[j] );
               }
     }

}

/* free_topology()
 *
 * this frees the topology table.
 */

void free_topology ( multipop *mpop )
{
     int i;
     for ( i = 0; i < mpop->exchanges; ++i )
     {
          if ( mpop->exch[i].from )
               FREE ( mpop->exch[i].from );
          if ( mpop->exch[i].as )
               FREE ( mpop->exch[i].as );
          FREE ( mpop->exch[i].fromsc );
     }
     if ( mpop->exch )
          FREE ( mpop->exch );
}

/* exchange_subpopulations()
 *
 * this performs the actual exchanges, using the information stored
 * in the exchange table.
 */

void exchange_subpopulations ( multipop *mpop )
{
     int i, j, k;
     sel_context *tocon;
     sel_context **fromcon;
     select_context_func_ptr select_con;
     int tp, *fp;
     int ti, *fi;

     /** arrays used for composite individuals. **/

     /* fromcon[j] holds the selection context used to pick individual
	to take tree j from. */
     fromcon = (sel_context **)MALLOC ( tree_count * sizeof ( sel_context * ) );
     /* fp[j] holds the population from which to take tree j from. */
     fp = (int *)MALLOC ( tree_count * sizeof ( int ) );
     /* fi[j] holds the individual from which to take tree j from. */
     fi = (int *)MALLOC ( tree_count * sizeof ( int ) );

     for ( i = 0; i < mpop->exchanges; ++i )
     {
#ifdef DEBUG
          printf ( "working on exch[%d]\n", i+1 );
#endif

	  /* where individuals are going. */
          tp = mpop->exch[i].to;

	  /* set up selection method to pick individuals to be replaced. */
          select_con = get_select_context ( mpop->exch[i].tosc );
          tocon = select_con ( SELECT_INIT, NULL, mpop->pop[tp],
                              mpop->exch[i].tosc );

	  /* are we copying whole individuals or creating composites? */
          if ( mpop->exch[i].copywhole > -1 )
          {
	       /*** copying whole individuals. ***/

	       /* the source subpop. */
               fp[0] = mpop->exch[i].copywhole;

	       /* selection method for choosing individuals from source
		  subpop. */
               select_con = get_select_context ( mpop->exch[i].fromsc[0] );
               fromcon[0] = select_con ( SELECT_INIT, NULL, mpop->pop[fp[0]],
                                        mpop->exch[i].fromsc[0] );

               for ( k = 0; k < mpop->exch[i].count; ++k )
               {
		    /** pick an individual to be replaced that has not already
		      been replaced during this exchange cycle. **/
                    do
                    {
                         ti = tocon->select_method ( tocon );
                    }
                    while ( mpop->pop[tp]->ind[ti].flags & FLAG_NEWEXCH );

		    /* pick an individual from the source subpop. */
                    fi[0] = fromcon[0]->select_method ( fromcon[0] );
                    
#ifdef DEBUG
                    printf ( "COPYING WHOLE INDIVIDUAL: ind %d subpop %d --> ind %d subpop %d\n",
                            fi[0], fp[0], ti, tp );
#endif
                    
                    /** remove the old iondividual from the population. **/
                    for ( j = 0; j < tree_count; ++j )
                    {
			 /* always dereference ERCs when removing trees
			    from the population. */
                         reference_ephem_constants ( mpop->pop[tp]->ind[ti].tr[j].data, -1 );
                         free_tree ( mpop->pop[tp]->ind[ti].tr+j );
                    }

		    /* copy the individual. */
                    duplicate_individual ( mpop->pop[tp]->ind+ti,
                                           mpop->pop[fp[0]]->ind+fi[0] );

		    /* reference the ERCs in the new individual. */
                    for ( j = 0; j < tree_count; ++j )
                         reference_ephem_constants ( mpop->pop[tp]->ind[ti].tr[j].data, 1 );

		    /* mark the individual as just coming from an exchange. */
                    mpop->pop[tp]->ind[ti].flags = FLAG_NEWEXCH;
               }

	       /* all done with this exchange, delete the selection context. */
               fromcon[0]->context_method ( SELECT_CLEAN, fromcon[0],
                                           NULL, NULL );
          }
          else
          {
               /*** creating composite individuals. ***/

	       /** create selection contexts for each tree. **/
               for ( j = 0; j < tree_count; ++j )
               {
		    /* does this tree need a context? */
                    if ( mpop->exch[i].fromsc[j] )
                    {
#ifdef DEBUG
                         printf ( "getting selection context for tree %d (%s)\n",
                                 j, mpop->exch[i].fromsc[j] );
#endif
			 /* create it. */
                         select_con = get_select_context ( mpop->exch[i].fromsc[j] );
                         fromcon[j] = select_con ( SELECT_INIT, NULL,
                                                  mpop->pop[mpop->exch[i].from[j]],
                                                  mpop->exch[i].fromsc[j] );
                    }
                    else
			 /* don't need one. */
                         fromcon[j] = NULL;
               }

               for ( k = 0; k < mpop->exch[i].count; ++k )
               {
		    /** select an individual to be replaced that hasn't already
		      been during this exchange cycle. **/
                    do
                    {
                         ti = tocon->select_method ( tocon );
                    }
                    while ( mpop->pop[tp]->ind[ti].flags & FLAG_NEWEXCH );
                    
#ifdef DEBUG
                    printf ( "SELECTED TREE %d FOR REPLACEMENT.\n", ti );
                    print_individual ( mpop->pop[tp]->ind+ti, stdout );
#endif
		    /** now select the individuals that we will merge to
		      replace trees of the destination individual. **/
                    for ( j = 0; j < tree_count; ++j )
                    {
			 /* we don't need to do a selection for a particular
			    tree if (1) it uses the same individual as another
			    tree or (2) it doesn't get replaced in the destination
			    individual. */

                         fp[j] = mpop->exch[i].from[j];
                         if ( fp[j] != -1 )
                         {
                              fi[j] = fromcon[fp[j]]->select_method ( fromcon[fp[j]] );
#ifdef DEBUG
                              printf ( "selecting using (%s) from subpop %d (for tree %d): individual %d\n",
                                      mpop->exch[i].fromsc[j], fp[j], j, fi[j] );
                              print_individual ( mpop->pop[fp[j]]->ind+fi[j], stdout );
#endif
                         }
                    }

		    /** now resolve "as_" references in the fp and fi arrays. */
                    for ( j = 0; j < tree_count; ++j )
                         if ( fp[j] == -1 )
                         {
                              if ( mpop->exch[i].as[j] == -1 )
				   /* tree j doesn't get replaced, so set both
				      values to -1. */
                                   fp[j] = fi[j] = -1;
                              else
                              {
				   /* tree j comes from the same individual as
				      some other tree. */
                                   fp[j] = fp[mpop->exch[i].as[j]];
                                   fi[j] = fi[mpop->exch[i].as[j]];
                              }
                         }

#ifdef DEBUG
                    printf ( "the fp,fi arrays are:\n" );
                    for ( j = 0; j < tree_count; ++j )
                         printf ( "   %3d:  fp = %3d    fi = %4d\n", j, fp[j], fi[j] );
#endif

                    /** replace the appropriate parts of the old tree. **/
                    for ( j = 0; j < tree_count; ++j )
                    {
			 /* skip trees that don't get replaced. */
                         if ( fp[j] == -1 )
                              continue;

			 /* dereference ERCs in old tree. */
                         reference_ephem_constants ( mpop->pop[tp]->ind[ti].tr[j].data, -1 );
			 /* delete old tree. */
                         free_tree ( mpop->pop[tp]->ind[ti].tr+j );
			 /* copy new tree. */
                         copy_tree ( mpop->pop[tp]->ind[ti].tr+j, mpop->pop[fp[j]]->ind[fi[j]].tr+j );
			 /* reference ERCs in new tree. */
                         reference_ephem_constants ( mpop->pop[tp]->ind[ti].tr[j].data, 1 );
                    }
		    /* evaluate the fitness of the new composite individual. */
                    app_eval_fitness ( mpop->pop[tp]->ind+ti );
                    mpop->pop[tp]->ind[ti].flags = FLAG_NEWEXCH;

#ifdef DEBUG
                    printf ( "the new individual is:\n" );
                    print_individual ( mpop->pop[tp]->ind+ti, stdout );
#endif
               }

	       /* destroy source selection contexts. */
               for ( j = 0; j < tree_count; ++j )
                    if ( fromcon[j] )
                         fromcon[j]->context_method ( SELECT_CLEAN,
                                                     fromcon[j], NULL, NULL );
          }

	  /* destroy destination selection context. */
          tocon->context_method ( SELECT_CLEAN, tocon, NULL, NULL );
     }

     FREE ( fromcon );
     FREE ( fp );
     FREE ( fi );

     /* erase all the NEWEXCH flags. */
     for ( i = 0; i < mpop->size; ++i )
          for ( j = 0; j < mpop->pop[i]->size; ++j )
               mpop->pop[i]->ind[j].flags &= ~FLAG_NEWEXCH;
     
}
                   
/* rebuild_exchange_topology()
 *
 * rebuilds the exchange table.  called from user code after making changes
 * to the parameters governing exchanges.
 */

void rebuild_exchange_topology ( multipop *mpop )
{
     free_topology ( mpop );
     initialize_topology ( mpop );
}
