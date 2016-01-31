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

#ifndef _PROTOS_H
#define _PROTOS_H

/*** exch.c ***/

void exchange_subpopulations ( multipop *mpop );
void initialize_topology ( multipop *mpop );
void free_topology ( multipop *mpop );
void rebuild_exchange_topology ( multipop *mpop );

/*** change.c ***/

population *change_population ( population *pop, breedphase * );
void show_population ( population *p );
breedphase * initialize_one_breeding ( char *prefix );
void initialize_breeding ( multipop * );
void free_one_breeding ( breedphase * );
void free_breeding ( multipop * );
void rebuild_breeding ( multipop * );
char *get_breed_parameter ( char *prefix, char *format, ... );


/*** ckpoint.c ***/

void read_checkpoint ( char *filename, int *gen, multipop **mpop );
void write_checkpoint ( int gen, multipop *mpop, char *filename );

population *read_population ( ephem_const **eind, FILE *f );
void read_individual ( individual *ind, ephem_const **eind, FILE *f,
		      char *buffer );
void write_individual ( individual *ind, ephem_index *eind, FILE *f );

void read_tree_recurse ( int space, ephem_const **eind, FILE *fil, int tree,
			char *string );
function * get_function_by_name ( int tree, char *string, ephem_const **ep,
				 ephem_const **eind );
void write_population ( population *pop, ephem_index *eind, FILE *f );
void write_tree_recurse ( lnode **l, ephem_index *eind, FILE *fil );

void write_hex_block ( void *, int, FILE * );
void read_hex_block ( void *, int, FILE * );


/*** ephem.c ***/

void initialize_ephem_const ( void );
void free_ephem_const ( void );
void enlarge_ephem_space ( void );
void ephem_const_gc ( void );
ephem_const *new_ephemeral_const ( function *f );
int ephem_index_comp ( const void *a, const void *b );
ephem_index *write_ephem_list ( FILE *f );
int lookup_ephem ( ephem_index *ind, ephem_const *e );
ephem_const **read_ephem_list ( FILE *f );
void get_ephem_stats ( int *used, int *free, int *blocks, int *alloc );



/*** eval.c ***/

void set_current_individual ( individual * );
DATATYPE evaluate_tree ( lnode *, int );
DATATYPE evaluate_tree_recurse ( lnode **, int );


/*** gp.c ***/

void run_gp ( multipop *mpop, int startgen,
             event *t_eval, event *t_breed, int startfromcheckpoint );
int generation_information ( int gen, multipop *mpop, int stt_interval,
                            int bestn );
void evaluate_pop ( population *pop );
int accumulate_pop_stats ( popstats *total, popstats *n );
void calculate_pop_stats ( popstats *s, population *pop, int gen, int subpop );
void saved_individual_gc ( void );
saved_ind ** write_saved_individuals ( ephem_index *eind, FILE *f );
void write_stats_checkpoint ( multipop *mpop, ephem_index *eind, FILE *f );
saved_ind ** read_saved_individuals ( ephem_const **eind, FILE *f );
void read_stats_checkpoint ( multipop *mpop, ephem_const **eind, FILE *f );


/*** main.c ***/

int function_sets_init ( function_set *, int, int *, char **, int );
int function_compare ( const void *a, const void *b );
void free_function_sets ( void );
void read_tree_limits ( void );
void initialize_random ( void );
void pre_parameter_defaults ( void );
void post_parameter_defaults ( void );
int process_commandline ( int argc, char **argv, int *gen, multipop ** );
void output_system_stats ( event *t_total, event *t_eval, event *t_breed );
void initial_message ( void );



/*** memory.c ***/

#ifdef TRACK_MEMORY
#define MALLOC track_malloc
#define FREE track_free
#define REALLOC track_realloc
#else
#define MALLOC malloc
#define FREE free
#define REALLOC realloc
#endif

void *track_malloc ( int );
void track_free ( void * );
void *track_realloc ( void *, int );
void get_memory_stats ( int *total, int *free, int *max,
                       int *mallocc, int *reallocc, int *freec );





/*** output.c ***/

int create_output_stream ( int id, char *ext, int reset, char *mode,
                          int autoflush );
void initialize_output_streams ( void );
void open_output_streams ( void );
void oputs ( int streamid, int detail, char *string );
void oprintf ( int streamid, int detail, char *format, ... );
FILE *output_filehandle ( int streamid );
void output_stream_close ( int streamid );
void output_stream_open ( int streamid );
void output_stream_flush ( int streamid );
void close_output_streams ( void );
void error ( int severity, char *format, ... );
void set_detail_level ( int );
int test_detail_level ( int );
void flush_output_streams ( void );



/*** params.c ***/

void initialize_parameters ( void );
void free_parameters ( void );
void add_parameter ( char *name, char *value, int copyflags );
int delete_parameter ( char *name );
char *get_parameter ( char *name );
void print_parameters ( void );
void write_parameter_database ( FILE *f );
void read_parameter_database ( FILE *f );
void read_parameter_file ( char * );
int delete_comment ( char * );
int check_continuation ( char * );
int parse_one_parameter ( char *buffer );
int trim_string ( char *string );
int translate_binary ( char *string );
void binary_parameter ( char *name, int value );
void define_directive ( char *string );
void undefine_directive ( char *string );
int test_directive ( char *string );


/*** populate.c ***/

void generate_random_population ( population *p, int *mindepth,
                                 int *maxdepth, int *method );
population *allocate_population ( int size );
void free_population ( population *p );
void free_multi_population ( multipop *mp );
population *initial_population ( int *, int *, int * );
multipop *initial_multi_population ( void );




/*** postscript.c ***/

int postscript_recurse ( lnode **, FILE *, int, int, int );
void make_postscript_tree ( lnode *, char *, int );



/*** random.c ***/

void random_seed ( long );
int random_int ( int );
double random_double ( void );
void *random_get_state ( int * );
void random_set_state ( void * );



/*** select.c ***/

int exists_select_method ( char *string );
select_context_func_ptr get_select_context ( char *string );
void free_o_rama ( int, char *** );
int parse_o_rama ( char *string, char *** argv );
int rev_ind_compare ( const void *a, const void *b );
int select_interval ( sel_context *sc );

/*** fitness.c ***/

sel_context *select_afit_context ( int op, sel_context *sc,
                                  population *p, char *string );
int select_afit ( sel_context *sc );
sel_context *select_inverse_afit_context ( int op, sel_context *sc,
                                          population *p, char *string );
int select_inverse_afit ( sel_context *sc );
sel_context *select_afit_overselect_context ( int op, sel_context *sc,
                                             population *p, char *string );
int select_afit_overselect ( sel_context *sc );

/*** tournament.c ***/

sel_context *select_tournament_context ( int op, sel_context *sc,
                                        population *p, char *string );
int select_tournament ( sel_context *sc );

/*** bestworst.c ***/

int select_bestworst ( sel_context *sc );
sel_context *select_best_context ( int op, sel_context *sc,
                                  population *p, char *string );
int select_best_compare ( const void *a, const void *b );
sel_context *select_worst_context ( int op, sel_context *sc,
                                   population *p, char *string );
int select_worst_compare ( const void *a, const void *b );
sel_context *select_random_context ( int op, sel_context *sc,
                                    population *p, char *string );
int select_random ( sel_context *sc );

/*** tree.c ***/

int tree_nodes ( lnode *tree );
int tree_nodes_recurse ( lnode ** );

int tree_nodes_internal ( lnode * );
int tree_nodes_internal_recurse ( lnode ** );

int tree_nodes_external ( lnode * );
int tree_nodes_external_recurse ( lnode ** );

void generate_random_full_tree ( int space, int depth, function_set * );
void generate_random_grow_tree ( int space, int depth, function_set * );

int tree_depth ( lnode * );
int tree_depth_recurse ( lnode ** );

int tree_depth_to_subtree ( lnode *, lnode * );
int tree_depth_to_subtree_recurse ( lnode **, lnode *, int );

void print_tree ( lnode *, FILE * );
void print_tree_recurse ( lnode **, FILE * );
void print_tree_array ( lnode * );
void print_tree_array_recurse ( lnode **, int * );

lnode *get_subtree ( lnode *, int );
lnode *get_subtree_recurse ( lnode **, int * );

lnode *get_subtree_internal ( lnode *, int );
lnode *get_subtree_internal_recurse ( lnode **, int * );

lnode *get_subtree_external ( lnode *, int );
lnode *get_subtree_external_recurse ( lnode **, int * );

void copy_tree ( tree *to, tree *from );
void free_tree ( tree * );

int tree_size ( lnode * );
int tree_size_recurse ( lnode ** );

void copy_tree_replace_many ( int space, lnode *parent, lnode **replace,
                            lnode **with, int count, int *repcount );
void copy_tree_replace_many_recurse ( int space, lnode **lp, lnode **lr,
                                    lnode **lw, int count, int *repcount );
void skip_over_subtree ( lnode ** );

void reference_ephem_constants ( lnode *, int );
void reference_ephem_constants_recurse ( lnode **, int );


/*** pretty.c ***/

void gen_indents ( lnode **l, int **is, int start, int sameline );
void pretty_print_tree_recurse ( lnode **l, int **is, FILE *fil );
void pretty_print_tree ( lnode *data, FILE *fil );

/*** genspace.c ***/

void initialize_genspace ( void );
void free_genspace ( void );
lnode * gensp_next ( int space );
int gensp_next_int ( int space );
void gensp_dup_tree ( int space, tree *t );
void gensp_reset ( int space );

/*** individ.c ***/

void print_individual ( individual *ind, FILE *f );
void pretty_print_individual ( individual *ind, FILE *f );
int individual_size ( individual *ind );
int individual_depth ( individual *ind );
void duplicate_individual ( individual *to, individual *from );

/*** crossover.c ***/

int operator_crossover_init ( char *options, breedphase *bp );
void operator_crossover_free ( void * );
void operator_crossover_start ( population *oldpop, void *data );
void operator_crossover_end ( void *data );
void operator_crossover ( population *oldpop, population *newpop, void *data );

/*** reproduce.c ***/

int operator_reproduce_init ( char *options, breedphase *bp );
void operator_reproduce_free ( void * );
void operator_reproduce_start ( population *oldpop, void *data );
void operator_reproduce_end ( void *data );
void operator_reproduce ( population *oldpop, population *newpop, void *data );

/*** mutate.c ***/

int operator_mutate_init ( char *options, breedphase *bp );
void operator_mutate_free ( void * );
void operator_mutate_start ( population *oldpop, void *data );
void operator_mutate_end ( void *data );
void operator_mutate ( population *oldpop, population *newpop, void *data );


extern genspace gensp[GENSPACE_COUNT];
extern function_set *fset;
extern int fset_count;
extern treeinfo *tree_map;
extern int tree_count;
extern int ind_nodelimit;

#endif
