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

#ifdef USEVFORK
extern char **environ;
#endif

/* read_checkpoint()
 *
 * reads a checkpoint file, placing the generation number in gen and filling
 * the population (and other structures) with information from the file
 */

void read_checkpoint ( char *filename, int *gen, multipop **mpop )
{
     FILE *f;
     char *buffer;
     ephem_const **eind;
     int random_state_bytes;
     int i, j;
     char *rand_state;

     /* miscellaneous buffer for reading. */
     buffer = (char *)MALLOC ( MAXCHECKLINELENGTH );

     /* open the file. */
     f = fopen ( filename, "rb" );
     if ( f == NULL )
     {
          error ( E_FATAL_ERROR, "couldn't read checkpoint \"%s\".",
                 filename );
     }

     oprintf ( OUT_SYS, 30, "reading from checkpoint \"%s\".\n",
              filename );
     
     /** confirm the magic word that starts every checkpoint file. **/
     fgets ( buffer, MAXCHECKLINELENGTH, f );
     if ( strcmp ( buffer, CK_MAGIC ) )
          error ( E_FATAL_ERROR,
                 "\"%s\" is not a lil-gp v1.0 checkpoint file.", filename );

     /* skip the human-readable id line. */
     fgets ( buffer, MAXCHECKLINELENGTH, f );
#ifdef DEBUG
     printf ( "id line: %s", buffer );
#endif

     /** read and print the timestamp. **/
     fscanf ( f, "%*s " );
     fgets ( buffer, MAXCHECKLINELENGTH, f );
     /* chop the newline. */
     buffer[strlen(buffer)-1] = 0;
     oprintf ( OUT_SYS, 30, "    checkpoint timestamp: [%s].\n", buffer );

     /* skip the "section: global" line. */
     fgets ( buffer, MAXCHECKLINELENGTH, f );
#ifdef DEBUG
     printf ( "should be global section: %s", buffer );
#endif

     /* read the generation number. */
     fscanf ( f, "%*s %d\n", gen );

     /** read the random number state encoded as a string of hex chars. **/

     /* first read the length. */
     fscanf ( f, "%*s %d ", &random_state_bytes );
#ifdef DEBUG
     fprintf ( stderr, "%d random state bytes.\n", random_state_bytes );
#endif
     /* allocate the buffer. */
     rand_state = (char *)MALLOC ( random_state_bytes+1 );
     /* read the hex data into the buffer. */
     read_hex_block ( rand_state, random_state_bytes, f );
     /* set the state. */
     random_set_state ( rand_state );
     /* free the buffer. */
     FREE ( rand_state );
     /* slurp the newline character following the hex data. */
     fgetc ( f );

     /** skip the "section: parameter" line. **/
     fgets ( buffer, MAXCHECKLINELENGTH, f );
#ifdef DEBUG
     printf ( "should be parameter section: %s", buffer );
#endif
     /* read the parameter database. */
     read_parameter_database ( f );

     /* make internal copies of function set(s). */
     if ( app_build_function_sets() ) 
          error ( E_FATAL_ERROR, "app_build_function_sets() failure." );

     /** skip the "section: erc" line. **/
     fgets ( buffer, MAXCHECKLINELENGTH, f );
#ifdef DEBUG
     printf ( "should be erc section: %s", buffer );
#endif
     
     /* read the list of ephemeral constants, and index them */
     eind = read_ephem_list ( f );

     /** skip the "section: erc" line. **/
     fgets ( buffer, MAXCHECKLINELENGTH, f );
#ifdef DEBUG
     printf ( "should be population section: %s", buffer );
#endif

     /** read the population **/

     /* allocate memory. */
     *mpop = (multipop *)MALLOC ( sizeof ( multipop ) );
     /* read number of subpops. */
     fscanf ( f, "%*s %d\n", &((**mpop).size) );
     /* allocate subpop list. */
     (**mpop).pop = (population **)MALLOC ( (**mpop).size *
                                           sizeof ( population * ) );
     for ( i = 0; i < (**mpop).size; ++i )
     {
	  /** skip each "subpop: #" line. **/
	  fgets ( buffer, MAXCHECKLINELENGTH, f );
#ifdef DEBUG
	  printf ( "should be subpop %d: %s", i, buffer );
#endif
	  /* read the population. */
          (**mpop).pop[i] = read_population ( eind, f );
     }

     /** skip the "section: application" line. **/
     fgets ( buffer, MAXCHECKLINELENGTH, f );
#ifdef DEBUG
     printf ( "should be application section: %s", buffer );
#endif
     /* read application-specific stuff. */
     app_read_checkpoint ( f );

     /** skip the "section: statistics" line. **/
     fgets ( buffer, MAXCHECKLINELENGTH, f );
#ifdef DEBUG
     printf ( "should be statistics section: %s", buffer );
#endif
     /* read the statistics. */
     read_stats_checkpoint ( *mpop, eind, f );

     /* close'n'free. */
     FREE ( eind );
     FREE ( buffer );
     fclose ( f );
     
     oprintf ( OUT_SYS, 30, "population read from checkpoint \"%s\".\n",
              filename );

}
     
/* write_checkpoint()
 *
 * checkpoints the population to the given file.
 */

void write_checkpoint ( int gen, multipop *mpop, char *filename )
{

     FILE *f;
     unsigned char *rand_state;
     ephem_index *eind;
     int i, j;
     int random_state_bytes;
     time_t now;
     char *param;
     char *compresscommand[4] = { NULL, NULL, NULL, NULL };

     /* open the file. */
     f = fopen ( filename, "w" );
     if ( f == NULL )
     {
          error ( E_ERROR, "couldn't write checkpoint \"%s\"; skipping.",
                 filename );
          return;
     }

     /* write magic number and id string. */
     fputs ( CK_MAGIC, f );
     fputs ( CK_IDSTRING, f );
     /* write timestamp. */
     time ( &now );
     fprintf ( f, "checkpoint-written: %s", ctime ( &now ) );

     /* global section. */
     fputs ( "section: global\n", f );
     fprintf ( f, "generation: %d\n", gen );
     
     /** write the state of the random number generator. **/
     rand_state = random_get_state ( &random_state_bytes );
     fprintf ( f, "random-state: %d ", random_state_bytes );
     /* store buffer as hex data. */
     write_hex_block ( rand_state, random_state_bytes, f );
     fputc ( '\n', f );
     FREE ( rand_state );

     /** write the parameter database. **/
     fprintf ( f, "section: parameter\n" );
     write_parameter_database ( f );

     /** write the list of ephemeral constants, and index them. **/
     fprintf ( f, "section: erc\n" );
     eind = write_ephem_list ( f );

     /** write the population. **/
     fprintf ( f, "section: population\n" );
     fprintf ( f, "subpop-count: %d\n", mpop->size );
     for ( i = 0; i < mpop->size; ++i )
     {
	  fprintf ( f, "subpop: %d\n", i );
          write_population ( mpop->pop[i], eind, f );
     }
     
     /** application-specific data. **/
     fprintf ( f, "section: application\n" );
     app_write_checkpoint ( f );

     /** statistics structures. **/
     fprintf ( f, "section: statistics\n" );
     write_stats_checkpoint ( mpop, eind, f );

     /** close'n'free. **/
     FREE ( eind );
     fclose ( f );

     oprintf ( OUT_SYS, 20, "    population checkpointed: \"%s\".\n",
              filename );

     /** do we compress the checkpoint file? **/
     param = get_parameter ( "checkpoint.compress" );
     if ( param )
     {
#if defined(USEVFORK) || defined(USESYSTEM)
	  /* allocate a string big enough to hold the command. */
	  compresscommand[2] = (char *)MALLOC ( 2*(strlen(param) +
						strlen(filename)) *
					    sizeof ( char ) );
	  /* create the command string. */
	  sprintf ( compresscommand[2], param, filename );
#ifdef DEBUG
	  oprintf ( OUT_SYS, 20, "    compression command is [%s]\n",
		   compresscommand[2] );
#endif
#ifdef USEVFORK
	  /* in unix (solaris at least), a system() call performs a fork(),
	   * then an exec().  the fork system call copies the entire address
	   * space of the parent process to the child process.  for large
	   * GP applications, this could be intolerably slow.
	   *
	   * vfork() does a fork without copying the address space.  it can
	   * be used when the child immediately exec()s following the
	   * vfork().
	   *
	   * we use exec() to do a "/bin/sh -c compresscommand" to parse
	   * and execute the compression command.
	   *
	   * we neither wait for the child to complete nor check the exit
	   * status to see if the compression was successful.
	   */

	  /** create the rest of the argv[] array to pass to the child. */
	  compresscommand[0] = "/bin/sh";
	  compresscommand[1] = "-c";
	  if ( !vfork() )
	  {
	       execve ( "/bin/sh", compresscommand, environ );
	       _exit(1);
	  }
#else
	  /* this is provided for non-unix systems which don't provide the
           * vfork() call but do have the system() call.
	   */

 	  system ( compresscommand[2] );
#endif
	  FREE ( compresscommand[2] );
	  oprintf ( OUT_SYS, 20, "    checkpoint compressed.\n" );
#else
	  /* neither vfork() nor system() is available;
	     can't do compression. */
	  oprintf ( OUT_SYS, 20, "    checkpoint compression unavailable.\n" );
#endif
     }
     
}

/* read_population()
 *
 * allocates a population structure, reads a population from a checkpoint
 * file into it, and returns it.  must be passed an index to look up ERCs
 * in.
 */

population *read_population ( ephem_const **eind, FILE *f )
{
     int i;
     lnode *l;
     char *buffer;
     population *pop;

     /* allocate. */
     pop = (population *)MALLOC ( sizeof ( population ) );
     /* read the "size" and "next" fields. */
     fscanf ( f, "%*s %d\n%*s %d\n", &(pop->size), &(pop->next) );
     /* allocate the individual array. */
     pop->ind = (individual *)MALLOC ( pop->size * sizeof ( individual ) );
     /* this buffer is used by read_individual for reading and parsing
	function names in trees.  we allocate it here, so that all calls
	to read_individual share the same buffer (we don't have to repeatedly
	allocate and free). */
     buffer = (char *)MALLOC ( MAXCHECKLINELENGTH );

     for ( i = 0; i < pop->size; ++i )
     {
	  read_individual ( pop->ind+i, eind, f, buffer );
     }

     FREE ( buffer );
     return pop;
}

/* read_individual()
 *
 * reads a single individual from a checkpoint file into the given individual
 * pointer.  it does NOT allocate the pointer.
 */

void read_individual ( individual *ind, ephem_const **eind, FILE *f,
		      char *buffer )
{
     int j, k[3];

     /* read the evald and flags fields. */
     fscanf ( f, "%d %d ", &(ind->evald), &(ind->flags) );
     if ( ind->evald == EVAL_CACHE_VALID )
     {
	  /** if the individual has valid fitness values saved in the
	    file, read them. **/

	  /* skip over the human-readable fitness values and read the
	     hits count. */
	  fscanf ( f, "%*f %*f %*f %d ", &(ind->hits) );
	  /** the fitness values, which are double precision, are dumped out
	    in hex so that no significant digits are lost. **/
	  read_hex_block ( &(ind->r_fitness), sizeof(double), f );
	  fgetc ( f );
	  read_hex_block ( &(ind->s_fitness), sizeof(double), f );
	  fgetc ( f );
	  read_hex_block ( &(ind->a_fitness), sizeof(double), f );
	  fgetc ( f );
     }
     
#ifdef DEBUG
     fprintf ( stderr, "%d: %lf %lf %lf %d %d %d\n", i,
	      ind->r_fitness,
	      ind->s_fitness,
	      ind->a_fitness,
	      ind->hits,
	      ind->evald,
	      ind->flags );
#endif

     /* allocate the array of trees. */
     ind->tr = (tree *)MALLOC ( tree_count * sizeof ( tree ) );
     for ( j = 0; j < tree_count; ++j )
     {
	  /* read the tree number, tree size, and tree node count. */
	  fscanf ( f, "%d %d %d ", k+0, k+1, k+2 );

	  /** read the tree into a generation space, then copy it to
	    it's final location. **/
	  gensp_reset ( 0 );
	  read_tree_recurse ( 0, eind, f, j, buffer );
	  gensp_dup_tree ( 0, ind->tr+j );
	  
#ifdef DEBUG_READTREE
	  fprintf ( stderr, "file: %d %d %d    here: %d %d %d\n",
		   k[0], k[1], k[2], j, ind->tr[j].size,
		   ind->tr[j].nodes );
	  print_tree ( ind->tr[j].data, stderr );
#endif
	  if ( k[0] != j ||
	       k[1] != ind->tr[j].size ||
	       k[2] != ind->tr[j].nodes )
	  {
	       /** if the values in the checkpoint file don't match the
		 values of the tree we read, this is a problem.  this, of
		 course should never happen. **/
	       error ( E_FATAL_ERROR, "checkpoint file corrupted in population section." );
	  }
     }
}

/* read_tree_recurse()
 *
 * function to recursively read a tree from a checkpoint file.
 */

void read_tree_recurse ( int space, ephem_const **eind, FILE *fil, int tree,
			char *string )
{
     function *f;
     int i, j;
     ephem_const *ep;

     /* read up until a nonwhitespace character in file.   the nonwhitespace
      character is saved in string[0]. */
     while ( isspace(string[0]=fgetc(fil)) );
     /* get the next character. */
     i = fgetc ( fil );
     if ( isspace(i) )
	  /* if the next character is whitespace, then string[0] is a
	     one-character function name.  null-terminate the string. */
	  string[1] = 0;
     else
     {
	  /** if the next character is not whitespace, then string[0]
	    is either an open parenthesis or the first character of a
	    multi-character function name. **/
	  
	  /* push the next character back. */
	  ungetc ( i, fil );
	  /* read the function name.  skip over an open parenthesis, if there
	     is one. */
	  fscanf ( fil, "%s ", string+(string[0]!='(') );
     }
#ifdef DEBUG_READTREE
     fprintf ( stderr, "function name is [%s]\n", string );
#endif

     /* look up the function name in this tree's function set.  if the
	function is an ERC terminal (the name is of the form "name:ERCindex"),
	then place the ERC address in ep. */
     f = get_function_by_name ( tree, string, &ep, eind );
     /* add an lnode to the tree. */
     gensp_next(space)->f = f;
     
     switch ( f->type )
     {
	case TERM_NORM:
	case TERM_ARG:
	case EVAL_TERM:
	  break;
	case TERM_ERC:
	  /* record the ERC address as the next lnode in the array. */
	  gensp_next(space)->d = ep;
	  break;
	case FUNC_DATA:
	case EVAL_DATA:
	  /** recursively read child functions, no skip nodes needed. **/
	  for ( i = 0; i < f->arity; ++i )
	       read_tree_recurse ( space, eind, fil, tree, string );
	  break;
	case FUNC_EXPR:
	case EVAL_EXPR:
	  /** recursively read child functions, recording skip values. **/
	  for ( i = 0; i < f->arity; ++i )
	  {
	       /* save an lnode for the skip value. */
	       j = gensp_next_int ( space );
	       /* read the child tree. */
	       read_tree_recurse ( space, eind, fil, tree, string );
	       /* figure out how big the child tree was, and save that
		  number in the skip node. */
	       gensp[space].data[j].s = gensp[space].used-j-1;
	  }
	  break;
     }
}

/* get_function_by_name()
 *
 * looks up a function name in the function set for the given tree.  if
 * the function is an ERC, looks up the index (encoded in the name)
 * and stores the ERC address in ep.
 */

function * get_function_by_name ( int tree, char *string, ephem_const **ep,
				 ephem_const **eind )
{
     int i, j, k;
     function_set *fs = fset+tree_map[tree].fset;

     k = strlen ( string );
     for ( i = 0; i < k; ++i )
     {
	  if ( string[i] == ':' )
	  {
	       /* names of the form "name:index" are chopped at the colon,
		  and the value of the index saved. */
	       string[i] = 0;
	       j = atoi ( string+i+1 );
	       break;
	  }
	  else if ( string[i] == ')' )
	  {
	       /* chop the name at the first closing parenthesis, since we
		  could be passed a string like "function))))" */
	       string[i] = 0;
	       break;
	  }
     }

     /* find the string in the function set. */
     for ( i = 0; i < fs->size; ++i )
	  if ( strcmp ( string, fs->cset[i].string ) == 0 )
	  {
	       if ( fs->cset[i].type == TERM_ERC )
	       {
		    /* if this is an ERC, lookup the saved index in the
		       eind table, and store the looked-up address in ep. */
		    *ep = eind[j];
		    (*ep)->f = fs->cset+i;
	       }
	       /* return a pointer to the function. */
	       return fs->cset+i;
	  }

     /* this, of course, should never happen. */
     return NULL;
}
     
/* write_population()
 *
 * writes a population to a checkpoint file.
 */

void write_population ( population *pop, ephem_index *eind, FILE *f )
{
     int i, j;

     /* write size and next fields. */
     fprintf ( f, "size: %d\nnext: %d\n", pop->size, pop->next );

     /* write each individual. */
     for ( i = 0; i < pop->size; ++i )
     {
	  write_individual ( pop->ind+i, eind, f );
     }
}

/* write_individual()
 *
 * writes an individual to a checkpoint file.  uses eind to change ERC
 * addresses to integer indices.
 */

void write_individual ( individual *ind, ephem_index *eind, FILE *f )
{
     int j;
     lnode *l;

     /* write evald and flags fields. */
     fprintf ( f, "%d %d ", ind->evald, ind->flags );
     if ( ind->evald == EVAL_CACHE_VALID )
     {
	  /** if the fitness values are valid... **/

	  /* ...write them in human-readable form. */
	  fprintf ( f, "%lf %lf %lf %d ",
		   ind->r_fitness, ind->s_fitness,
		   ind->a_fitness, ind->hits );
	  /** then write the double-precision values as hex blocks, so
	     as not to lose significant digits. **/
	  write_hex_block ( &(ind->r_fitness), sizeof(double), f );
	  fputc ( ' ', f );
	  write_hex_block ( &(ind->s_fitness), sizeof(double), f );
	  fputc ( ' ', f );
	  write_hex_block ( &(ind->a_fitness), sizeof(double), f );
     }
     fputc ( '\n', f );

     /** now write the trees of the individual. **/
     for ( j = 0; j < tree_count; ++j )
     {
	  /* write tree number, size, nodes. */
	  fprintf ( f, "%d %d %d ", j, ind->tr[j].size,
		   ind->tr[j].nodes );

	  /** write tree data. **/
	  l = ind->tr[j].data;
	  write_tree_recurse ( &l, eind, f );
	  fputc ( '\n', f );
     }
}

/* write_tree_recurse()
 *
 * function to recursively write trees to a checkpoint file.  the same
 * as print_tree_recurse(), except that ERC nodes are written as
 * "name:index" rather than the value, using eind to translate addresses
 * to indices.
 */

void write_tree_recurse ( lnode **l, ephem_index *eind, FILE *fil )
{
     function *f;
     int i;

     /* remember which function we are. */
     f = (**l).f;

     /* a space, then an open-paren if this function is not a terminal. */
     fputc ( ' ', fil );
     if ( f->arity != 0 )
          fprintf ( fil, "(" );
     
     ++*l;
     if ( f->type == TERM_ERC )
     {
	  /* ERCs printed as "name:index". */
	  fprintf ( fil, "%s:%d", f->string,
		   lookup_ephem ( eind, (**l).d ) );
          ++*l;
     }
     else
	  /* everything else printed normally. */
          fprintf ( fil, "%s", f->string );
     
     switch ( f->type )
     {
        case FUNC_DATA:
        case EVAL_DATA:
	  /** recursively print children. **/
          for ( i = 0; i < f->arity; ++i )
               write_tree_recurse ( l, eind, fil );
          break;
        case FUNC_EXPR:
        case EVAL_EXPR:
	  /** recursive print children, ignoring the skip nodes. **/
          for ( i = 0; i < f->arity; ++i )
          {
               ++*l;
               write_tree_recurse ( l, eind, fil );
          }
          break;
     }

     if ( f->arity != 0 )
          fprintf ( fil, ")" );
}

/* write_hex_block()
 *
 * writes a block of memory to a file, as a string of hex characters.
 */

void write_hex_block ( void *buf, int n, FILE *f )
{
     int i;
     unsigned char *b = (unsigned char *)buf;
     
     for ( i = 0; i < n; ++i )
	  fprintf ( f, "%02x", b[i] );
}

/* read_hex_block()
 *
 * reads hex characters into a block of memory, the inverse of
 * write_hex_block().
 */

void read_hex_block ( void *buf, int n, FILE *f )
{
     int i;
     unsigned char *b = (unsigned char *)buf;
     int c[2] = { 0, 0 };
     
     for ( i = 0; i < n; ++i )
     {
	  c[0] = fgetc ( f );
	  c[1] = fgetc ( f );
	  
	  /* convert hex chars to base 10. */
	  c[0] = c[0]>'9' ? c[0]-'a'+10 : c[0]-'0';
	  c[1] = c[1]>'9' ? c[1]-'a'+10 : c[1]-'0';
	  
	  b[i] = c[0] * 16 + c[1];
     }
}

	  
	  
