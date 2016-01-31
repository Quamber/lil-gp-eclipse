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

#include <stdio.h>
#include <math.h>

#include "kernel/lilgp.h"

globaldata g;

static int fitness_cases = -1;
static double *app_fitness_cases[2];
static double value_cutoff;

int app_build_function_sets ( void )
{
     function_set fset;
     int tree_map;
     char *tree_name;
     function sets[10] = 
     { { f_multiply,   NULL,      NULL,         2, "*",     FUNC_DATA, -1, 0 },
       { f_protdivide, NULL,      NULL,         2, "/",     FUNC_DATA, -1, 0 },
       { f_add,        NULL,      NULL,         2, "+",     FUNC_DATA, -1, 0 },
       { f_subtract,   NULL,      NULL,         2, "-",     FUNC_DATA, -1, 0 },
       { f_sin,        NULL,      NULL,         1, "sin",   FUNC_DATA, -1, 0 },
       { f_cos,        NULL,      NULL,         1, "cos",   FUNC_DATA, -1, 0 },
       { f_exp,        NULL,      NULL,         1, "exp",   FUNC_DATA, -1, 0 },
       { f_rlog,       NULL,      NULL,         1, "rlog",  FUNC_DATA, -1, 0 },
       { f_indepvar,   NULL,      NULL,         0, "X",     TERM_NORM, -1, 0 },
       { NULL,         f_erc_gen, f_erc_print,  0, "R",     TERM_ERC,  -1, 0 } };
     
     binary_parameter ( "app.use_ercs", 1 );
     if ( atoi ( get_parameter ( "app.use_ercs" ) ) )
          fset.size = 10;
     else
          fset.size = 9;
     fset.cset = sets;

     tree_map = 0;
     tree_name = "TREE";

     return function_sets_init ( &fset, 1, &tree_map, &tree_name, 1 );
}

void app_eval_fitness ( individual *ind )
{

     int i;
     double v, dv;
     double disp;

     set_current_individual ( ind );
     
     ind->r_fitness = 0.0;
     ind->hits = 0;
     
     for ( i = 0; i < fitness_cases; ++i )
     {
          g.x = app_fitness_cases[0][i];
          v = evaluate_tree ( ind->tr[0].data, 0 );
          dv = app_fitness_cases[1][i];
          disp = fabs ( dv-v );
          
          if ( disp < value_cutoff )
          {
               ind->r_fitness += disp;
               if ( disp <= 0.01 )
                    ++ind->hits;
          }
          else
          {
               ind->r_fitness += value_cutoff;
          }
     }
     ind->s_fitness = ind->r_fitness;
     ind->a_fitness = 1/(1+ind->s_fitness);
     ind->evald = EVAL_CACHE_VALID;
}

int app_end_of_evaluation ( int gen, multipop *mpop, int newbest,
                           popstats *gen_stats, popstats *run_stats )
{
     int i;
     double v;
                    
     if ( newbest )
     {
          output_stream_open ( OUT_USER );
          
          for ( i = -100; i <= 100; ++i )
          {
               g.x = (double)i*.01;
               v = evaluate_tree ( run_stats[0].best[0]->ind->tr[0].data, 0 );
               oprintf ( OUT_USER, 50, "%lf %lf\n", g.x, v );
          }
          
          output_stream_close ( OUT_USER );
          
          if ( run_stats[0].best[0]->ind->hits == fitness_cases )
               return 1;
     }
     
     return 0;
}

void app_end_of_breeding ( int gen, multipop *mpop )
{
     return;
}

int app_create_output_streams ( void )
{
     if ( create_output_stream ( OUT_USER, ".fn", 1, "w", 0 ) !=
         OUTPUT_OK )
          return 1;

     return 0;
}

int app_initialize ( int startfromcheckpoint )
{
     int i;
     double x, y;
     char *param;

     if ( !startfromcheckpoint )
     {
	  oprintf ( OUT_PRG, 50, "not starting from checkpoint file.\n" );
	  
	  param = get_parameter ( "app.fitness_cases" );
	  if ( param == NULL )
	       fitness_cases = 200;
	  else
	  {
	       fitness_cases = atoi ( param );
	       if ( fitness_cases < 0 )
		    error ( E_FATAL_ERROR, "invalid value for \"app.fitness_cases\"." );
	  }
	  
	  app_fitness_cases[0] = (double *)MALLOC ( fitness_cases *
						   sizeof ( double ) );
	  app_fitness_cases[1] = (double *)MALLOC ( fitness_cases *
						   sizeof ( double ) );
	  
	  oprintf ( OUT_PRG, 50, "%d fitness cases:\n", fitness_cases );
	  for ( i = 0; i < fitness_cases; ++i )
	  {
	       x = (random_double()*2.0)-1.0;
	       
	       /* change this line to modify the goal function. */
	       y = x*x*x*x + x*x*x + x*x + x;
	       
	       app_fitness_cases[0][i] = x;
	       app_fitness_cases[1][i] = y;
     
	       oprintf ( OUT_PRG, 50,
			"    x = %12.5lf, y = %12.5lf\n", x, y );
	  }
     }
     else
     {
	  oprintf ( OUT_PRG, 50, "started from checkpoint file.\n" );
     }

     param = get_parameter ( "app.value_cutoff" );
     if ( param == NULL )
          value_cutoff = 1.e15;
     else
          value_cutoff = strtod ( param, NULL );
     
     return 0;
}

void app_uninitialize ( void )
{
     FREE ( app_fitness_cases[0] );
     FREE ( app_fitness_cases[1] );
}

void app_write_checkpoint ( FILE *f )
{
     int i;
     fprintf ( f, "fitness-cases: %d\n", fitness_cases );
     for ( i = 0; i < fitness_cases; ++i )
     {
	  write_hex_block ( app_fitness_cases[0]+i, sizeof(double), f );
	  fputc ( ' ', f );
	  write_hex_block ( app_fitness_cases[1]+i, sizeof(double), f );
	  fprintf ( f, " %.5lf %.5lf\n",
		   app_fitness_cases[0][i], app_fitness_cases[1][i] );
     }
}

void app_read_checkpoint ( FILE *f )
{
     int i;
     
     fscanf ( f, "%*s %d\n", &fitness_cases );
     
     app_fitness_cases[0] = (double *)MALLOC ( fitness_cases *
                                              sizeof ( double ) );
     app_fitness_cases[1] = (double *)MALLOC ( fitness_cases *
                                              sizeof ( double ) );

     for ( i = 0; i < fitness_cases; ++i )
     {
	  read_hex_block ( app_fitness_cases[0]+i, sizeof ( double ), f );
	  fgetc ( f );
	  read_hex_block ( app_fitness_cases[1]+i, sizeof ( double ), f );
	  fscanf ( f, " %*f %*f\n" );
	  fprintf ( stderr, "%.5lf %.5lf\n", app_fitness_cases[0][i],
		   app_fitness_cases[1][i] );
     }
}
