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
     int keep_trying;
     double internal;
     double external;
     double *tree;       /* probability that a given tree
			    will be selected for crossover. */
     double *treecumul;  /* running sum of "tree" field. */
     double treetotal;   /* total of all tree fields. */
     double *func;       /* probability that a given function
			    set will be selected for crossover. */
     char *sname;
     sel_context *sc;
     char *sname2;
     sel_context *sc2;
} crossover_data;

/* operator_crossver_init()
 *
 * called to parse crossover options and initialize one record
 * of a breedphase table appropriately.
 */

int operator_crossover_init ( char *options, breedphase *bp )
{
     int errors = 0;
     crossover_data *cd;
     int i, j, k, m;
     double r;
     char **argv, **targv;
     int internalset = 0, externalset = 0;
     char *cp;

     cd = (crossover_data *)MALLOC ( sizeof ( crossover_data ) );

     /* place values into the breedphase table record. */
     bp->operator = OPERATOR_CROSSOVER;
     bp->data = (void *)cd;
     bp->operator_free = operator_crossover_free;
     bp->operator_start = operator_crossover_start;
     bp->operator_end = operator_crossover_end;
     bp->operator_operate = operator_crossover;

     /* default values for all the crossover options. */
     cd->keep_trying = 0;
     cd->internal = 0.9;
     cd->external = 0.1;
     cd->tree = (double *)MALLOC ( tree_count * sizeof ( double ) );
     cd->treecumul = (double *)MALLOC ( tree_count * sizeof ( double ) );
     for ( j = 0; j < tree_count; ++j )
          cd->tree[j] = 0.0;
     cd->treetotal = 0.0;
     cd->func = (double *)MALLOC ( fset_count * sizeof ( double ) );
     cd->sname = NULL;
     cd->sname2 = NULL;

     /* break the options string into an argv-style array of strings. */
     j = parse_o_rama ( options, &argv );

     for ( i = 0; i < j; ++i )
     {
	  /* parse "keep_trying" option. */
          if ( strcmp ( "keep_trying", argv[i] ) == 0 )
          {
	       /* translate a string into a binary value.  returns -1 if
		  the string is not one of the valid strings meaning
		  yes or no. */
               cd->keep_trying = translate_binary ( argv[++i] );
               if ( cd->keep_trying == -1 )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"%s\" is not a valid setting for \"keep_trying\".",
                           argv[i] );
               }
          }
	  /* parse "internal" option. */
          else if ( strcmp ( "internal", argv[i] ) == 0 )
          {
               internalset = 1;
               cd->internal = strtod ( argv[++i], NULL );
               if ( cd->internal < 0.0 )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"internal\" must be nonnegative." );
               }
          }
	  /* parse "external" option. */
          else if ( strcmp ( "external", argv[i] ) == 0 )
          {
               externalset = 1;
               cd->external = strtod ( argv[++i], NULL );
               if ( cd->external < 0.0 )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"external\" must be nonnegative." );
               }
          }
	  /* parse "select" option. */
          else if ( strcmp ( "select", argv[i] ) == 0 )
          {
               if ( !exists_select_method ( argv[++i] ) )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"%s\" is not a known selection method.",
                           argv[i] );
               }
               FREE ( cd->sname );
               cd->sname = (char *)MALLOC ( (strlen(argv[i])+1) * sizeof ( char ) );
               strcpy ( cd->sname, argv[i] );
               if ( cd->sname2 == NULL )
                    cd->sname2 = cd->sname;
          }
	  /* parse "select2" option. */
          else if ( strcmp ( "select2", argv[i] ) == 0 )
          {
               if ( !exists_select_method ( argv[++i] ) )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"%s\" is not a known selection method.",
                           argv[i] );
               }
               if ( cd->sname2 && cd->sname != cd->sname2 )
                    FREE ( cd->sname2 );
               cd->sname2 = (char *)MALLOC ( (strlen(argv[i])+1) * sizeof ( char ) );
               strcpy ( cd->sname2, argv[i] );
          }
	  /* parse "tree" option. */
          else if ( strcmp ( "tree", argv[i] ) == 0 )
          {
               k = parse_o_rama ( argv[++i], &targv );
               if ( k != tree_count )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: wrong number of tree fields: \"%s\".",
                           argv[i] );
               }
               else
               {
                    for ( m = 0; m < k; ++m )
                    {
                         cd->tree[m] = strtod ( targv[m], &cp );
                         if ( *cp )
                         {
                              ++errors;
                              error ( E_ERROR, "crossover: \"%s\" is not a number.",
                                     targv[m] );
                         }
                    }
               }
               
               free_o_rama ( k, &targv );
          }
	  /* parse "tree#" option. */
          else if ( strncmp ( "tree", argv[i], 4 ) == 0 )
          {
               k = strtol ( argv[i]+4, &cp, 10 );
               if ( *cp )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: unknown option \"%s\".",
                           argv[i] );
               }
               if ( k < 0 || k >= tree_count )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"%s\" is out of range.",
                           argv[i] );
               }
               else
               {
                    cd->tree[k] = strtod ( argv[++i], &cp );
                    if ( *cp )
                    {
                         ++errors;
                         error ( E_ERROR, "crossover: \"%s\" is not a number.",
                                argv[i] );
                    }
               }
          }
          else
          {
               ++errors;
               error ( E_ERROR, "crossover: unknown option \"%s\".",
                      argv[i] );
          }
     }
     
     free_o_rama ( j, &argv );
     
     if ( internalset && !externalset )
          cd->external = 0.0;
     else if ( !internalset && externalset )
          cd->internal = 0.0;
     
     if ( cd->sname == NULL )
     {
          ++errors;
          error ( E_ERROR, "crossover: no selection method specified." );
     }

     /** compute "func" array from the "tree" array. **/
     
     for ( j = 0; j < tree_count; ++j )
          cd->treetotal += cd->tree[j];
     if ( cd->treetotal == 0.0 )
     {
          for ( j = 0; j < tree_count; ++j )
               cd->tree[j] = 1.0;
          cd->treetotal = tree_count;
     }
          
     for ( j = 0; j < fset_count; ++j )
          cd->func[j] = 0.0;
     for ( j = 0; j < tree_count; ++j )
          cd->func[tree_map[j].fset] += cd->tree[j];
     
     r = 0.0;
     for ( j = 0; j < fset_count; ++j )
          r = (cd->func[j] += r);

#ifdef DEBUG
     if ( !errors )
     {
          printf ( "crossover options:\n" );
          printf ( "   internal: %lf  external: %lf\n", cd->internal, cd->external );
          printf ( "   keep_trying: %d\n", cd->keep_trying );
          printf ( "   primary selection: %s\n", cd->sname==NULL?"NULL":cd->sname );
          printf ( "   second selection: %s\n", cd->sname2==NULL?"NULL":(cd->sname2==cd->sname?"same as primary":cd->sname2) );
          printf ( "   tree total: %lf\n", cd->treetotal );
          for ( j = 0; j < tree_count; ++j )
               printf ( "   tree %d: %lf\n", j, cd->tree[j] );
          for ( j = 0; j < fset_count; ++j )
               printf ( "   fset %d: %lf\n", j, cd->func[j] );
     }
#endif
     
     return errors;
}

/* operator_crossover_free()
 *
 * free the crossover-specific data structure.
 */

void operator_crossover_free ( void *data )
{
     crossover_data * cd;

     cd = (crossover_data *)data;

     FREE ( cd->sname );
     if ( cd->sname != cd->sname2 )
          FREE ( cd->sname2 );
     FREE ( cd->tree );
     FREE ( cd->treecumul );
     FREE ( cd->func );
     FREE ( cd );
}

/* operator_crossover_start()
 *
 * called at the start of the breeding process each generation.
 * initializes the selection contexts for this phase.
 */

void operator_crossover_start ( population *oldpop, void *data )
{
     crossover_data * cd;
     select_context_func_ptr select_con;

     cd = (crossover_data *)data;
     
     select_con = get_select_context ( cd->sname );
     cd->sc = select_con ( SELECT_INIT, NULL, oldpop, cd->sname );

     /* if there is a separate selection method specified for the
	second parent... */
     if ( cd->sname2 != cd->sname )
     {
	  /* ...then initialize it too. */
          select_con = get_select_context ( cd->sname2 );
          cd->sc2 = select_con ( SELECT_INIT, NULL, oldpop, cd->sname2 );
     }
     else
	  /* ...otherwise use the first context. */
          cd->sc2 = cd->sc;
}

/* operator_crossover_end()
 *
 * called when breeding is finished each generation.  frees up selection
 * contexts for this phase.
 */

void operator_crossover_end ( void *data )
{
     crossover_data * cd;

     cd = (crossover_data *)data;

     cd->sc->context_method ( SELECT_CLEAN, cd->sc, NULL, NULL );
     if ( cd->sname != cd->sname2 )
          cd->sc2->context_method ( SELECT_CLEAN, cd->sc2, NULL, NULL );
}

/* operator_crossover()
 *
 * performs the crossover, inserting one or both offspring into the
 * new population.
 */

void operator_crossover ( population *oldpop, population *newpop,
                         void *data )
{
     crossover_data * cd;
     int p1, p2;
     int ps1, ps2;
     int l1, l2;
     lnode *st[3];
     int sts1, sts2;
     int ns1, ns2;
     int badtree1, badtree2;
     double total;
     int forceany1, forceany2;
     int repcount;
     int f, t1, t2, j;
     double r, r2;
     int totalnodes1, totalnodes2;
     int i;

     /* get the crossover-specific data structure. */
     cd = (crossover_data *)data;
     total = cd->internal + cd->external;

     /* choose a function set. */
     r = random_double() * cd->treetotal;
     for ( f = 0; r >= cd->func[f]; ++f );

     /* fill in the "treecumul" array, zeroing all trees which
	don't use the selected function set. */
     r = 0.0;
     t1 = 0;
     for ( j = 0; j < tree_count; ++j )
     {
          if ( tree_map[j].fset == f )
               r = (cd->treecumul[j] = r + cd->tree[j]);
          else
               cd->treecumul[j] = r;
     }

     /* select the first and second trees. */
     r2 = random_double() * r;
     for ( t1 = 0; r2 >= cd->treecumul[t1]; ++t1 );
     r2 = random_double() * r;
     for ( t2 = 0; r2 >= cd->treecumul[t2]; ++t2 );

#ifdef DEBUG_CROSSOVER
     printf ( "selected function set %d --> t1: %d; t2: %d\n", f, t1, t2 );
#endif
     
     /* choose two parents */
     p1 = cd->sc->select_method ( cd->sc );
     ps1 = oldpop->ind[p1].tr[t1].nodes;
     /* if the tree only has one node, we obviously can't do
	fucntionpoint crossover.  use anypoint instead. */
     forceany1 = (ps1==1||total==0.0);
     
     p2 = cd->sc2->select_method ( cd->sc2 );
     ps2 = oldpop->ind[p2].tr[t2].nodes;
     forceany2 = (ps2==1||total==0.0);

#ifdef DEBUG_CROSSOVER
     fprintf ( stderr, "parent 1 is:\n" );
     print_individual ( oldpop->ind+p1, stderr );
     fprintf ( stderr, "parent 2 is:\n" );
     print_individual ( oldpop->ind+p2, stderr );
#endif

     while(1)
     {
          
          /* choose two crossover points */

          if ( forceany1 )
          {
	       /* choose any point. */
               l1 = random_int ( ps1 );
               st[1] = get_subtree ( oldpop->ind[p1].tr[t1].data, l1 );
          }
          else if ( total*random_double() < cd->internal )
          {
	       /* choose an internal point. */
               l1 = random_int ( tree_nodes_internal (oldpop->ind[p1].tr[t1].data) );
               st[1] = get_subtree_internal ( oldpop->ind[p1].tr[t1].data, l1 );
          }
          else
          {
	       /* choose an external point. */
               l1 = random_int ( tree_nodes_external (oldpop->ind[p1].tr[t1].data) );
               st[1] = get_subtree_external ( oldpop->ind[p1].tr[t1].data, l1 );
          }
                                
          if ( forceany2 )
          {
	       /* choose any point on second parent. */
               l2 = random_int ( ps2 );
               st[2] = get_subtree ( oldpop->ind[p2].tr[t2].data, l2 );
          }
          else if ( total*random_double() < cd->internal )
          {
	       /* choose internal point. */
               l2 = random_int ( tree_nodes_internal (oldpop->ind[p2].tr[t2].data) );
               st[2] = get_subtree_internal ( oldpop->ind[p2].tr[t2].data, l2 );
          }
          else
          {
	       /* choose external point. */
               l2 = random_int ( tree_nodes_external (oldpop->ind[p2].tr[t2].data) );
               st[2] = get_subtree_external ( oldpop->ind[p2].tr[t2].data, l2 );
          }

#ifdef DEBUG_CROSSOVER
          printf ( "subtree 1 is: " );
          print_tree ( st[1], stdout );
          printf ( "subtree 2 is: " );
          print_tree ( st[2], stdout );
#endif

	  /* count the nodes in the selected subtrees. */
          sts1 = tree_nodes ( st[1] );
          sts2 = tree_nodes ( st[2] );

	  /* calculate the sizes of the offspring. */
          ns1 = ps1 - sts1 + sts2;
          ns2 = ps2 - sts2 + sts1;

          totalnodes1 = ns1;
          totalnodes2 = ns2;

#ifdef DEBUG_CROSSOVER
          printf ( "newtree 1 has size %d; limit is %d\n",
                  ns1, tree_map[t1].nodelimit );
#endif

	  /** validate the first offspring against the tree node and depth
	    limits; set "badtree1" if any are violated. **/
	  
          badtree1 = 0;
          if ( tree_map[t1].nodelimit > -1 && ns1 > tree_map[t1].nodelimit )
               badtree1 = 1;
          else if ( tree_map[t1].depthlimit > -1 )
          {
               ns1 = tree_depth_to_subtree ( oldpop->ind[p1].tr[t1].data, st[1] ) +
                     tree_depth ( st[2] );
#ifdef DEBUG_CROSSOVER
               printf ( "newtree 1 has depth %d; limit is %d\n",
                       ns1, tree_map[t1].depthlimit );
#endif
               if ( ns1 > tree_map[t1].depthlimit )
                    badtree1 = 1;
          }

	  /* if we're supposed to keep trying, skip up and choose new crossover
	     points. */
          if ( cd->keep_trying && badtree1 )
               continue;

	  /** validate the second offspring against the tree node and depth
	    limits; set "badtree2" if any are violated. **/
	  
          badtree2 = 0;
          if ( tree_map[t2].nodelimit > -1 && ns2 > tree_map[t2].nodelimit )
               badtree2 = 1;
          else if ( tree_map[t2].depthlimit > -1 )
          {
               ns2 = tree_depth_to_subtree ( oldpop->ind[p2].tr[t2].data, st[2] ) +
                     tree_depth ( st[1] );
               if ( ns2 > tree_map[t2].depthlimit )
                    badtree2 = 1;
          }

	  /* if we're supposed to keep trying, skip up and choose new crossover
	     points. */
	  
          if ( cd->keep_trying && badtree2 )
               continue;

	  /* check both offspring against the individual node limits, set
	     badtree1 and/or badtree2 if either is too big. */
	  
          if ( ind_nodelimit > -1 )
          {
               for ( i = 0; i < tree_count; ++i )
               {
                    if ( i != t1 )
                         totalnodes1 += oldpop->ind[p1].tr[i].nodes;
                    if ( i != t2 )
                         totalnodes2 += oldpop->ind[p2].tr[i].nodes;
               }
               badtree1 |= (totalnodes1 > ind_nodelimit);
               badtree2 |= (totalnodes2 > ind_nodelimit);
#ifdef DEBUG_CROSSOVER
               printf ( "newind 1 has %d nodes; limit is %d\n",
                       totalnodes1, ind_nodelimit );
#endif
          }

	  /* choose new crossover points if either tree is too big. */
          if ( cd->keep_trying && (badtree1 || badtree2) )
               continue;
          
          /* copy the first parent to the first offspring position */
          duplicate_individual ( newpop->ind+newpop->next,
                                oldpop->ind+p1 );
          if ( !badtree1 )
          {
	       /* if the first offspring is allowable... */
	       
#ifdef DEBUG_CROSSOVER
               fprintf ( stderr, "offspring 1 is allowable.\n" );
#endif
	       
	       /* make a copy of the crossover tree, replacing the
		  selected subtree with the crossed-over subtree. */
               copy_tree_replace_many ( 0, oldpop->ind[p1].tr[t1].data,
                                       st+1, st+2, 1, &repcount );
               if ( repcount != 1 )
               {
		    /* this can't happen, but check anyway. */
                    error ( E_FATAL_ERROR,
                           "botched crossover:  this can't happen" );
               }

               /* free the appropriate tree of the new individual */
               free_tree ( newpop->ind[newpop->next].tr+t1 );
               /* copy the crossovered tree to the freed space */
               gensp_dup_tree ( 0, newpop->ind[newpop->next].tr+t1 );

	       /* the new individual's fitness fields are of course invalid. */
               newpop->ind[newpop->next].evald = EVAL_CACHE_INVALID;
               newpop->ind[newpop->next].flags = FLAG_NONE;
          }
          else
          {
	       /* offspring too big but keep_trying not set, just leave the copied
		  parent 1 in the offspring position. */
#ifdef DEBUG_CROSSOVER
               fprintf ( stderr, "offspring 1 is too big; copying parent 1.\n" );
#endif
          }

	  /* we've just filled in one member of the new population. */
          ++newpop->next;
          
#ifdef DEBUG_CROSSOVER
          fprintf ( stderr, "offspring 1:" );
          if ( newpop->ind[newpop->next-1].evald == EVAL_CACHE_VALID )
               fprintf ( stderr, "  (valid)\n" );
          else
               fprintf ( stderr, "  (invalid)\n" );
          print_individual ( newpop->ind+(newpop->next-1), stderr );
#endif          

	  /* if the new population needs another member (it's not full) */
          if ( newpop->next < newpop->size )
          {
	       /* copy the second parent to the second offspring position. */
               duplicate_individual ( newpop->ind+newpop->next,
                                     oldpop->ind+p2 );
               if ( !badtree2 )
               {
		    /* if the second offspring is allowable... */
#ifdef DEBUG_CROSSOVER
                    fprintf ( stderr, "offspring 2 is allowable.\n" );
#endif
		    /* then make a copy of the tree, replacing the crossover
		       subtree. */
                    copy_tree_replace_many ( 0, oldpop->ind[p2].tr[t2].data,
                                            st+2, st+1, 1, &repcount );
                    if ( repcount != 1 )
                    {
                         error ( E_FATAL_ERROR,
                                "bad crossover:  this can't happen" );
                    }

		    /* free the old tree in the new individual, and replace
		       it with the crossover tree. */
                    free_tree ( newpop->ind[newpop->next].tr+t2 );
                    gensp_dup_tree ( 0, newpop->ind[newpop->next].tr+t2 );
                    
                    newpop->ind[newpop->next].evald = EVAL_CACHE_INVALID;
                    newpop->ind[newpop->next].flags = FLAG_NONE;
               }
               else
               {
		    /* offspring too big but keep_trying not set; just leave
		       the copy of parent 2 where it is. */
#ifdef DEBUG_CROSSOVER
                    fprintf ( stderr, "offspring 2 is big; copying parent 2.\n" );
#endif
               }
               
               ++newpop->next;
               
#ifdef DEBUG_CROSSOVER
               fprintf ( stderr, "offspring 2:" );
               if ( newpop->ind[newpop->next-1].evald == EVAL_CACHE_VALID )
                    fprintf ( stderr, "  (valid)\n" );
               else
                    fprintf ( stderr, "  (invalid)\n" );
               print_individual ( newpop->ind+(newpop->next-1), stderr );
#endif          
          
          }
#ifdef DEBUG_CROSSOVER
          else
          {
	       /* the first offspring filled the population, discard the
		  second. */
               fprintf ( stderr, "offspring 2 not needed.\n\n" );
          }
#endif

          break;
     }

#ifdef DEBUG_CROSSOVER
     printf ( "CROSSOVER COMPLETE.\n\n\n" );
#endif
     
}

