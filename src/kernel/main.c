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

#ifdef MEMORY_LOG
FILE *mlog;
#endif

/* do we dup OUT_SYS to stdout? */
int quietmode = 0;

/* tree generation spaces. */
genspace gensp[GENSPACE_COUNT];

/* internal copy of function set(s). */
function_set *fset;
int fset_count;

/* information about each tree--which function set it uses,
   its name, size limits, etc. */
treeinfo *tree_map;
int tree_count;

/* maximum number of nodes per individual.  -1 if no limit is
   enforced. */
int ind_nodelimit;

int main ( int argc, char **argv )
{

     multipop *mpop;
     int startgen;
     char *param;
     event start, end, diff;
     event eval, breed;
     int startfromcheckpoint;

#ifdef MEMORY_LOG
     /* dump all memory allocations to a file. */
     mlog = fopen ( "memory.log", "w" );
#endif

     /* mark the start time, and zero the accumulators for evaluation
	and breeding time. */
     event_init();
     event_mark ( &start );
     event_zero ( &eval );
     event_zero ( &breed );

     if ( app_create_output_streams() )
          error ( E_FATAL_ERROR, "app_create_output_streams() failure." );
     initialize_output_streams();

     /* print copyright message and such. */
     initial_message();

     /* some initialization. */
     oprintf ( OUT_SYS, 30, "initialization:\n" );
     initialize_parameters();
     initialize_ephem_const();
     initialize_genspace();
     /* process the command line.  if starting from a checkpoint file, this
	function will load the population. */
     startfromcheckpoint = process_commandline ( argc, argv, &startgen, &mpop );

     /* open the files associated with each stream. */
     open_output_streams();
     
     /* make internal copies of function set(s), if it hasn't already been
	done. */
     if ( !startfromcheckpoint )
	  if ( app_build_function_sets() ) 
	       error ( E_FATAL_ERROR, "app_build_function_sets() failure." );

     /* read parameters limiting tree node count and/or depth. */
     read_tree_limits();

     /* if not starting from a checkpoint, seed the random number generator. */
     if ( !startfromcheckpoint )
	  initialize_random();
     
     if ( app_initialize ( startfromcheckpoint ) )
          error ( E_FATAL_ERROR, "app_initialize() failure." );

     /* if not starting from a checkpoint, create a random population. */
     if ( !startfromcheckpoint )
          mpop = initial_multi_population();

     /* build the breeding table and the subpop exchange table from
	the parameter database. */
     initialize_topology ( mpop );
     initialize_breeding ( mpop );

     /* do the GP. */
     run_gp ( mpop, startgen, &eval, &breed, startfromcheckpoint );

     /* free app stuff. */
     app_uninitialize();

     /* free lots of stuff. */
     free_breeding ( mpop );
     free_topology ( mpop );
     free_multi_population ( mpop );
     free_parameters();
     free_ephem_const();
     free_genspace();
     free_function_sets();

     /* mark the finish time. */
     event_mark ( &end );
     event_diff ( &diff, &start, &end );

     /* print memory/time statistics and close output files. */
     output_system_stats ( &diff, &eval, &breed );
     close_output_streams();

#ifdef MEMORY_LOG
     fclose ( mlog );
#endif

     /* all done. */
     return 0;
}

void localSort (void *base, unsigned long len, unsigned long width,
	  int (*compar) (const void *, const void *))
{
  int index1, index2;
  void *record1, *record2, *temp;

  temp = MALLOC (width);
  for (index1 = 0; index1 < len-1; index1++)
  {
    record1 = (void*) ((char*)base + width * index1);

    for (index2 = index1+1; index2 < len; index2++)
    {
      record2 = (void*) ((char*)base + width * index2);

      if (compar((const void*)record2, (const void*) record1) < 0)
      {
	/* swap two structures */

	memcpy (temp, record1, width);
	memcpy (record1, record2, width);
	memcpy (record2, temp, width);
      }
    }
  }

  FREE (temp);
}


/* function_sets_init()
 *
 * this function is called from user code and passed the function sets
 * and tree maps and names.  it makes internal copies and does some
 * validation.
 */



int function_sets_init ( function_set *user_fset, int user_fcount,
                        int *user_tree_map, char **user_tree_name,
                        int user_tcount )
{
     int i, j, k, m, n, p;
     int errors = 0;
     function *cur;

     /* allocate internal copies. */
     fset = (function_set *)MALLOC ( user_fcount * sizeof ( function_set ) );
     tree_map = (treeinfo *)MALLOC ( user_tcount * sizeof ( treeinfo ) );
     fset_count = user_fcount;
     tree_count = user_tcount;

     oprintf ( OUT_SYS, 30, "building function set(s):\n" );

     /* for each set of functions... */
     for ( i = 0; i < fset_count; ++i )
     {
          oprintf ( OUT_SYS, 30, "    set %d:", i );

          /* allocate memory for the set. */
          m = user_fset[i].size;
          fset[i].cset = (function *)MALLOC ( m * sizeof ( function ) );
          fset[i].num_args = 0;
          
          k = 0;
          for ( j = 0; j < m; ++j )
          {
               if ( user_fset[i].cset[j].type == FUNC_DATA ||
                    user_fset[i].cset[j].type == FUNC_EXPR ||
                    user_fset[i].cset[j].type == EVAL_DATA ||
                    user_fset[i].cset[j].type == EVAL_EXPR )
               {
		    /** functions and evaluation tokens **/
		    
                    cur = &(fset[i].cset[k]);
                    
		    /* copy some stuff over. */
                    cur->code = user_fset[i].cset[j].code;
                    cur->ephem_gen = user_fset[i].cset[j].ephem_gen;
                    cur->ephem_str = user_fset[i].cset[j].ephem_str;
                    cur->arity = user_fset[i].cset[j].arity;
                    cur->type = user_fset[i].cset[j].type;
                    cur->evaltree = user_fset[i].cset[j].evaltree;

		    /* copy the name string. */
                    n = strlen ( user_fset[i].cset[j].string );
                    cur->string = (char *)MALLOC ( n+1 );
		    for ( p = 0; p < n; ++p )
		    {
			 if ( isspace(user_fset[i].cset[j].string[p]) ||
			      user_fset[i].cset[j].string[p] == ':' ||
			      user_fset[i].cset[j].string[p] == ')' ||
			      user_fset[i].cset[j].string[p] == '(' ||
			      user_fset[i].cset[j].string[p] == '[' ||
			      user_fset[i].cset[j].string[p] == ']' )
			 {
			      error ( E_WARNING, "illegal character(s) in function name changed to '_'." );
			      cur->string[p] = '_';
			 }
			 else
			      cur->string[p] = user_fset[i].cset[j].string[p];
		    }
                    cur->string[n] = 0;

		    /* fill in the index field with this function's position in the
		       set. */
                    cur->index = k;

		    /* the ERC-related fields should be NULL. */
                    if ( cur->ephem_gen || cur->ephem_str )
                    {
                         ++errors;
                         error ( E_ERROR, "function has non-NULL ephem_gen and/or ephem_str field(s)." );
                    }

		    /* do some type-specific checking. */
                    switch ( cur->type )
                    {
                       case FUNC_DATA:
                       case FUNC_EXPR:
                         if ( cur->code == NULL )
                         {
                              ++errors;
                              error ( E_ERROR, "ordinary function has NULL code field." );
                         }
                         if ( cur->arity < 1 )
                         {
                              ++errors;
                              error ( E_ERROR, "ordinary function has arity of %d.",
                                     cur->arity );
                         }
                         if ( cur->evaltree != -1 )
                         {
                              error ( E_WARNING, "ordinary function has evaltree field of %d; this will be ignored.", cur->evaltree );
                         }
                         break;
                       case EVAL_DATA:
                       case EVAL_EXPR:
                         if ( cur->code != NULL )
                         {
                              ++errors;
                              error ( E_ERROR, "eval function function has non-NULL code field." );
                         }
                         if ( cur->arity != -1 )
                         {
                              error ( E_WARNING, "eval function has arity field of %d; this will be ignored.", cur->arity );
                         }
                         if ( cur->evaltree < 0 ||
                              cur->evaltree >= tree_count )
			      /* evaluation token refers to a tree that doesn't exist. */
                              error ( E_FATAL_ERROR, "eval function refers to nonexistent tree (%d).", cur->evaltree );
                         break;
                       default:
                         ++errors;
                         error ( E_ERROR, "unknown function type %d.", cur->type );
                    }
                    ++k;
               }
               else if ( user_fset[i].cset[j].type == TERM_NORM ||
                        user_fset[i].cset[j].type == TERM_ERC  ||
                        user_fset[i].cset[j].type == TERM_ARG )
               {
		    /** terminals (all kinds). **/

		    /* "cur" is so much easier to type.  :)  */
                    cur = &(fset[i].cset[k]);
                    
		    /* copy stuff. */
                    cur->code = user_fset[i].cset[j].code;
                    cur->ephem_gen = user_fset[i].cset[j].ephem_gen;
                    cur->ephem_str = user_fset[i].cset[j].ephem_str;
                    cur->arity = user_fset[i].cset[j].arity;
                    cur->type = user_fset[i].cset[j].type;
                    cur->evaltree = user_fset[i].cset[j].evaltree;

		    /* copy terminal name. */
                    n = strlen ( user_fset[i].cset[j].string );
                    cur->string = (char *)MALLOC ( n+1 );
                    strcpy ( cur->string, user_fset[i].cset[j].string );
                    cur->string[n] = 0;

		    /* fill in the index field. */
                    cur->index = k;

                    if ( cur->arity != 0 )
                    {
                         ++errors;
                         error ( E_ERROR, "terminal has nonzero arity." );
                    }

		    /* check for correctness of type-dependent fields. */
                    switch ( cur->type )
                    {
                       case TERM_NORM:
                         if ( cur->code == NULL )
                         {
                              ++errors;
                              error ( E_ERROR, "normal terminal has NULL code field." );
                         }
                         if ( cur->ephem_gen != NULL || cur->ephem_str != NULL )
                         {
                              ++errors;
                              error ( E_ERROR, "normal terminal has non-NULL ephem_gen and/or ephem_str field(s)." );
                         }
                         if ( cur->evaltree != -1 )
                         {
                              error ( E_WARNING, "normal terminal has evaltree field of %d; this will be ignored." );
                         }
                         break;
                       case TERM_ERC:
                         if ( cur->code != NULL )
                         {
                              ++errors;
                              error ( E_ERROR, "ERC terminal has non-NULL code field." );
                         }
                         if ( cur->ephem_gen == NULL || cur->ephem_str == NULL )
                         {
                              ++errors;
                              error ( E_ERROR, "ERC terminal has NULL ephem_hen and/or ephem_str field(s)." );
                         }
                         if ( cur->evaltree != -1 )
                         {
                              error ( E_WARNING, "ERC terminal has evaltree field of %d; this will be ignored." );
                         }
                         break;
                       case TERM_ARG:
                         ++fset[i].num_args;
                         if ( cur->code != NULL )
                         {
                              ++errors;
                              error ( E_ERROR, "argument terminal has non-NULL code field." );
                         }
                         if ( cur->ephem_gen != NULL || cur->ephem_str != NULL )
                         {
                              ++errors;
                              error ( E_ERROR, "argument terminal has non-NULL ephem_hen and/or ephem_str field(s)." );
                         }
                         if ( cur->evaltree < 0 )
                         {
                              ++errors;
                              error ( E_ERROR, "argument terminal should have nonnegative evaltree field." );
                         }
                         break;
                    }
                    ++k;
               }
	       oputs ( OUT_SYS, 30, " " );
	       oputs ( OUT_SYS, 30, fset[i].cset[k-1].string );
          }
          fset[i].size = k;
	  oputs ( OUT_SYS, 30, "\n" );
     }

     /* if there were any errors, stop now. */
     if ( errors )
     {
          error ( E_FATAL_ERROR, "error(s) occurred while processing function set(s)." );
     }

     /* build the internal tree map. */
     for ( i = 0; i < tree_count; ++i )
     {
	  /* the function set used for this tree. */
          tree_map[i].fset = user_tree_map[i];
          if ( tree_map[i].fset < 0 || tree_map[i].fset >= fset_count )
               error ( E_FATAL_ERROR, "tree %d uses a nonexistent function set.\n", i );

          oprintf ( OUT_SYS, 30, "    tree %d uses function set %d.\n", i, tree_map[i].fset );

	  /* these will be filled in by read_tree_limits(). */
          tree_map[i].nodelimit = -1;
          tree_map[i].depthlimit = -1;

	  /* copy the tree name. */
          j = strlen ( user_tree_name[i] );
          tree_map[i].name = (char *)MALLOC ( (j+1) * sizeof ( char ) );
          strcpy ( tree_map[i].name, user_tree_name[i] );
     }

     /* now some more processing on each function set. */
     for ( i = 0; i < fset_count; ++i )
     {
          fset[i].function_count = 0;
          fset[i].terminal_count = 0;
          for ( j = 0; j < fset[i].size; ++j )
          {
               if ( fset[i].cset[j].arity == -1 )
               {
		    /* change the arity of evaluation tokens from -1
		       to the number of argument tokens in the called tree. */
                    fset[i].cset[j].arity = fset[tree_map[fset[i].cset[j].evaltree].fset].num_args;
                    if ( fset[i].cset[j].arity == 0 )
			 /* if there are no argument tokens in the tree,
			    mark this as a terminal. */
                         fset[i].cset[j].type = EVAL_TERM;
               }

	       /* update count of functions and terminals. */
               if ( fset[i].cset[j].arity )
                    ++fset[i].function_count;
               else
                    ++fset[i].terminal_count;
          }

	  /* now sort the function set so that all the functions
	     come first. */
	  /*          qsort ( fset[i].cset, fset[i].size, sizeof ( function ),
                 function_compare ); */
	  localSort ( fset[i].cset, fset[i].size, sizeof ( function ),
                 function_compare );
     }
          

#ifdef DEBUG
     /* dump the function sets to stdout. */
     for ( i = 0; i < fset_count; ++i )
     {
          printf ( "FUNCTION SET %d\n", i );
          printf ( "   %d functions; %d terminals; %d arguments\n",
                  fset[i].function_count, fset[i].terminal_count,
                  fset[i].num_args );

          for ( j = 0; j < fset[i].size; ++j )
               printf ( "%10s %06x %06x %06x arity: %3d evaltree: %3d index: %3d type: %3d\n",
                       fset[i].cset[j].string,
                       fset[i].cset[j].code,
                       fset[i].cset[j].ephem_gen,
                       fset[i].cset[j].ephem_str,
                       fset[i].cset[j].arity,
                       fset[i].cset[j].evaltree,
                       fset[i].cset[j].index,
                       fset[i].cset[j].type );
     }
#endif
                       
     oprintf ( OUT_SYS, 30, "    function set complete.\n" );

     return 0;
}

/* function_compare()
 *
 * comparison function for qsort() that puts all functions ahead of
 * all terminals.
 */

int function_compare ( const void *a, const void *b )
{
     int aa, ba;
     aa = ((function *)a)->arity ? 1 : 0;
     ba = ((function *)b)->arity ? 1 : 0;
     return ba-aa;
}

/* free_function_sets()
 *
 * free up internal copies of function sets and tree maps.
 */

void free_function_sets ( void )
{
     int i, j;

     for ( i = 0; i < fset_count; ++i )
     {
          for ( j = 0; j < fset[i].function_count + fset[i].terminal_count; ++j )
               FREE ( fset[i].cset[j].string );
          FREE ( fset[i].cset );
     }
     FREE ( fset );

     for ( i = 0; i < tree_count; ++i )
          FREE ( tree_map[i].name );
     FREE ( tree_map );
     
     fset = NULL;
     tree_map = NULL;
}

/* read_tree_limits()
 *
 * read limits on tree node count and/or depth from the parameter
 * database and fill in the appropriate fields of the tree_map
 * array.
 */

void read_tree_limits ( void )
{
     int i, j;
     char pnamebuf[100];
     char *param;

     for ( i = 0; i < tree_count; ++i )
     {
	  /* read the node limit for this tree. */
          sprintf ( pnamebuf, "tree[%d].max_nodes", i );
          param = get_parameter ( pnamebuf );
          if ( param == NULL )
               tree_map[i].nodelimit = -1;
          else
               tree_map[i].nodelimit = atoi ( param );

	  /* read the depth limit for this tree. */
          sprintf ( pnamebuf, "tree[%d].max_depth", i );
          param = get_parameter ( pnamebuf );
          if ( param == NULL )
               tree_map[i].depthlimit = -1;
          else
               tree_map[i].depthlimit = atoi ( param );
     }

     /* read the node limit for the whole individual. */
     param = get_parameter ( "max_nodes" );
     if ( param == NULL )
          ind_nodelimit = -1;
     else
          ind_nodelimit = atoi ( param );

     /* read the depth limit for the whole individual.  note that
	this is implemented just as a cap on the maximum depth of
	any single tree in the individual. */
     param = get_parameter ( "max_depth" );
     if ( param )
     {
          j = atoi ( param );
          if ( j >= 0 )
               for ( i = 0; i < tree_count; ++i )
                    if ( tree_map[i].depthlimit < 0 ||
                        tree_map[i].depthlimit > j )
                         tree_map[i].depthlimit = j;
     }
}

/* initialize_random()
 *
 * initialize the random number generator.
 */
                    
void initialize_random ( void )
{
     char *param;
     long seed;

     /* look for a seed parameter. */
     param = get_parameter ( "random_seed" );
     if ( param == NULL )
     {
	  /* if it's not found... */
#ifdef RANDOMSEEDTIME
	  /* ...use the current time. */
	  seed = time(NULL);
#else
	  /* ...use 1. */
	  seed = 1;
#endif
	  /* print out what we're using. */
	  oprintf ( OUT_SYS, 20,
		   "    no random number seed specfied; using %d.\n",
		   seed );
     }
     else
     {
	  /* the parameter was found; use it. */
	  seed = atol ( param );
	  oprintf ( OUT_SYS, 20,
		   "    seeding random number generator with %d.\n",
		   seed );
     }
     random_seed ( seed );
}

/* pre_parameter_defaults()
 *
 * used to place values into the parameter database before any application
 * code is called or any command line options are processed.
 */

void pre_parameter_defaults ( void )
{
     add_parameter ( "output.basename",          "lilgp", PARAM_COPY_NONE );
     add_parameter ( "output.stt_interval",      "1", PARAM_COPY_NONE );
     add_parameter ( "output.detail",            "50", PARAM_COPY_NONE );
     add_parameter ( "output.bestn",             "1", PARAM_COPY_NONE );
     add_parameter ( "output.digits",            "4", PARAM_COPY_NONE );
     
     add_parameter ( "init.method",              "half_and_half",
                    PARAM_COPY_NONE );
     add_parameter ( "init.depth",               "2-6", PARAM_COPY_NONE );
     add_parameter ( "init.random_attempts",     "100", PARAM_COPY_NONE );
     
     add_parameter ( "checkpoint.filename",      "gp%06d.ckp",
                    PARAM_COPY_NONE );
     
     /* default problem uses a single population. */
     add_parameter ( "multiple.subpops", "1", PARAM_COPY_NONE );
}

/* post_parameter_defaults()
 *
 * add/change values in the parameter database after all command line options
 * are parsed.  can be used to set defaults based on values of other
 * parameters.
 */

void post_parameter_defaults ( void )
{
     binary_parameter ( "probabilistic_operators", 1 );
}

/* process_commandline()
 *
 * parses the command line.
 */

int process_commandline ( int argc, char **argv, int *gen,
                          multipop **mpop )
{
     int i;
     int errorflag = 0;
     int startfromcheckpoint = 0;

     *mpop = NULL;
     *gen = 0;

     /* if there are no arguments, print out a brief statement of usage
	and exit. */
     if ( argc < 2 )
     {
          fprintf ( stderr, "usage: %s options\nValid options are:\n", argv[0] );
	  fprintf ( stderr, "      [-f parameterfile]     read named parameter file\n" );
	  fprintf ( stderr, "      [-c checkpointfile]    restart from name checkpoint file\n" );
	  fprintf ( stderr, "      [-p name=value]        set parameter name to value\n" );
	  fprintf ( stderr, "      [-q]                   run in quiet mode\n" );
	  fprintf ( stderr, "      [-d symbol]            define symbol\n" );
	  fprintf ( stderr, "      [-u symbol]            undefine symbol\n" );
          exit(1);
     }
     
     /* load hardcoded defaults into database. */
     pre_parameter_defaults();

     for ( i = 1; i < argc; ++i )
     {
	  /* all options begin with '-' and have two characters,
	     except "-d" and "-u" which may have more. */
          if ( argv[i][0] != '-' || ( argv[i][1] != 'd' && argv[i][1] != 'u' && argv[i][2] != 0 ) )
          { 
               error ( E_ERROR, "unrecognized command line option: \"%s\".",
                      argv[i] );
               errorflag = 1;
               continue;
          }

          switch ( argv[i][1] )
          {
             case 'f':
	       /* load a parameter file, named in the next argument. */
               read_parameter_file ( argv[++i] );
               break;
             case 'p':
	       /* parse a single parameter, in the next argument. */
               if ( parse_one_parameter ( argv[++i] ) )
               {
                    errorflag = 1;
                    error ( E_ERROR, "malformed parameter: \"%s\".", argv[i] );
               }
               break;
             case 'c':
	       /* load a checkpoint file, named in the next argument. */
               if ( startfromcheckpoint )
               {
		    /* error if a checkpoint has already been loaded. */
                    error ( E_ERROR, "can't load multiple checkpoint files." );
                    errorflag = 1;
                    continue;
               }
               read_checkpoint ( argv[++i], gen, mpop );
               startfromcheckpoint = 1;
               break;
             case 'q':
	       /* turn on quiet mode (don't dup OUT_SYS to stdout). */
               quietmode = 1;
               break;
             case 'd':
	       /* define a symbol. */
               if ( argv[i][2] )
		    /* of the form "-dsymbol". */
                    define_directive ( argv[i]+2 );
               else
		    /* of the form "-d symbol". */
                    define_directive ( argv[++i] );
               break;
             case 'u':
	       /* undefine a symbol. */
               if ( argv[i][2] )
		    /* of the form "-usymbol". */
                    undefine_directive ( argv[i]+2 );
               else
		    /* of the form "-u symbol". */
                    undefine_directive ( argv[++i] );
               break;
             default:
               error ( E_ERROR, "unrecognized command line option: \"%s\".",
                      argv[i] );
               errorflag = 1;
               break;
          }
     }

     if ( errorflag )
          error ( E_FATAL_ERROR, "command line errors occurred.  dying." );

     if ( !startfromcheckpoint )
	  post_parameter_defaults();

     return startfromcheckpoint;
     
}

/* output_system_stats()
 *
 * print statistics about memory use, execution time, etc. to OUT_SYS
 * at conclusion of run.
 */

void output_system_stats ( event *t_total, event *t_eval, event *t_breed )
{
     int total, free, max, mallocc, reallocc, freec;
     int ercused, ercfree, ercblocks, ercalloc;
     int i;

     get_ephem_stats ( &ercused, &ercfree, &ercblocks, &ercalloc );
     
     oprintf ( OUT_SYS, 30, "\nSYSTEM STATISTICS\n" );

#ifdef TRACK_MEMORY
     /* if memory tracking available, then get and print the numbers. */
     get_memory_stats ( &total, &free, &max, &mallocc, &reallocc, &freec );
     oprintf ( OUT_SYS, 30, "\n------- memory -------\n" );
     oprintf ( OUT_SYS, 30, "           allocated:      %d\n", total );
     oprintf ( OUT_SYS, 30, "               freed:      %d\n", free );
     oprintf ( OUT_SYS, 30, "           not freed:      %d\n", total-free );
     oprintf ( OUT_SYS, 30, "       max allocated:      %d\n", max );
     oprintf ( OUT_SYS, 30, "    malloc'ed blocks:      %d\n", mallocc );
     oprintf ( OUT_SYS, 30, "   realloc'ed blocks:      %d\n", reallocc );  
     oprintf ( OUT_SYS, 30, "      free'ed blocks:      %d\n", freec );
#endif

#ifdef TIMING_AVAILABLE
     /* if timing is available, the get and print the numbers. */
     oprintf ( OUT_SYS, 30, "\n------- time -------\n" );
     oprintf ( OUT_SYS, 30, "             overall:      %s\n",
              event_string ( t_total ) );
     oprintf ( OUT_SYS, 30, "          evaluation:      %s\n",
              event_string ( t_eval ) );
     oprintf ( OUT_SYS, 30, "            breeding:      %s\n",
              event_string ( t_breed ) );
#endif     

     /* show how large the generation spaces grew. */
     oprintf ( OUT_SYS, 30, "\n------- generation spaces -------\n" );
     for ( i = 0; i < GENSPACE_COUNT; ++i )
          oprintf ( OUT_SYS, 30, "      space %3d size:      %d\n",
                   i, gensp[i].size );

     /* if any ERCs were used, then show these stats. */
     if ( ercused > 0 )
     {
          oprintf ( OUT_SYS, 30, "\n------- ephemeral random constants -------\n" );
          oprintf ( OUT_SYS, 30, "                used:      %d\n", ercused );
          oprintf ( OUT_SYS, 30, "               freed:      %d\n", ercfree );
          oprintf ( OUT_SYS, 30, "           allocated:      %d\n", ercalloc );
          oprintf ( OUT_SYS, 30, "              blocks:      %d\n", ercblocks );
     }
}   

/* initial_message()
 *
 * show startup and copyright messages.
 */

void initial_message ( void )
{
     oputs ( OUT_SYS, 0,
            "\n[ lil-gp Genetic Programming System.\n" );
     oputs ( OUT_SYS, 0,
            "[ Portions copyright (c) 1995 Michigan State University.  All rights reserved.\n" );
     oputs ( OUT_SYS, 0,
            "[ kernel version 1.0; 11 July 1995.\n\n\n" );
              
}
