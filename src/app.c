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
#include <stdbool.h>

#include "kernel/lilgp.h"

globaldata g;

 int fitness_cases = -1;
 double *app_fitness_cases[3];
static int *app_fitness_importance;
static double value_cutoff;
multipop *mpop;
int startgen;
event start, end, diff;
event eval, breed;
int startfromcheckpoint;
int population_No = 0;
int generation_No = 0;
int populationSIZE = 1;
float current_top = 1000;
int generationSIZE = 1;
//int minimum=10000;
float *optimal_in_generation;
int *optimal_index_in_generation;
float **error_array;
int same_optimal_count = 1;
int current_max_importance = 1;
int best_starting = -1;
int best_ending = 1;
//static int fittestTreeNo=0;
//static int change_counter=0;

void init() {
	populationSIZE = atoi(get_parameter("pop_size"));
	generationSIZE = atoi(get_parameter("max_generations"));
	best_starting = atoi(get_parameter("fn_start"));
	best_ending = atoi(get_parameter("fn_end"));
	print_parameters();
	printf("Best %d : %d", best_starting, best_ending);

	if (best_starting < best_ending) {
		best_starting = -1;
		best_ending = 1;
	}
	printf("%d : %d", populationSIZE, generationSIZE);
	error_array = malloc(generationSIZE * sizeof(float *));
	optimal_in_generation = malloc(generationSIZE * sizeof(float));
	optimal_index_in_generation = malloc(generationSIZE * sizeof(int));
	for (int i = 0; i < generationSIZE; i++) {
		optimal_in_generation[i] = -1;
	}
	//memset(optimal_in_generation,-1,generationSIZE*sizeof(float));
	optimal_in_generation[0] = 1000;
	//
	int i;
	for (i = 0; i < generationSIZE; i++) {
		error_array[i] = malloc(populationSIZE * sizeof(float));
		memset(error_array[i], 0.0f, populationSIZE * sizeof(float));

	}
}
int app_build_function_sets(void) {
	function_set fset;
	int tree_map;
	char *tree_name;
	function sets[10] = { { f_multiply, NULL, NULL, 2, "*", FUNC_DATA, -1, 0 },
			{ f_protdivide, NULL, NULL, 2, "/", FUNC_DATA, -1, 0 }, { f_add,
			NULL, NULL, 2, "+", FUNC_DATA, -1, 0 }, { f_subtract, NULL,
			NULL, 2, "-", FUNC_DATA, -1, 0 }, { f_sin, NULL, NULL, 1, "sin",
			FUNC_DATA, -1, 0 }, { f_cos, NULL, NULL, 1, "cos",
			FUNC_DATA, -1, 0 }, { f_exp, NULL, NULL, 1, "exp",
			FUNC_DATA, -1, 0 }, { f_rlog, NULL, NULL, 1, "rlog",
			FUNC_DATA, -1, 0 }, { f_indepvar, NULL, NULL, 0, "X",
			TERM_NORM, -1, 0 }, { NULL, f_erc_gen, f_erc_print, 0, "R",
			TERM_ERC, -1, 0 } };

	binary_parameter("app.use_ercs", 1);
	if (atoi(get_parameter("app.use_ercs")))
		fset.size = 10;
	else
		fset.size = 9;
	fset.cset = sets;

	tree_map = 0;
	tree_name = "TREE";

	return function_sets_init(&fset, 1, &tree_map, &tree_name, 1);
}

void app_eval_fitness(individual *ind) {

	int i;
	double v, dv;
	double disp;
	float error = 0.0f;
	set_current_individual(ind);
	ind->r_fitness = 0.0;
	ind->hits = 0;

	for (i = 0; i < fitness_cases; ++i) {
		//	if (app_fitness_importance[i] <= current_max_importance&&app_fitness_importance[i] !=0) {
		g.x = app_fitness_cases[0][i];
		v = evaluate_tree(ind->tr[0].data, 0);
		dv = app_fitness_cases[1][i];
		disp = fabs(dv - v);
		error += disp;
		if (disp < value_cutoff) {
			ind->r_fitness += disp;
			if (disp <= 0.01)
				++ind->hits;
		} else {
			ind->r_fitness += value_cutoff;
		}
		//}
	}
	error = error / fitness_cases;
	//error = error/
	//  error_array[(generation_No*50)+population] = error;
	error_array[generation_No][population_No] = error;

	if (optimal_in_generation[generation_No] > error) {
		optimal_in_generation[generation_No] = error;
		optimal_index_in_generation[generation_No] = population_No;
	}
	ind->s_fitness = ind->r_fitness;
	ind->a_fitness = 1 / (1 + ind->s_fitness);
	ind->evald = EVAL_CACHE_VALID;

}

int app_end_of_evaluation(int gen, multipop *mpop, int newbest,
		popstats *gen_stats, popstats *run_stats) {
	int i;
	double v;

	if (newbest) {
		output_stream_open( OUT_USER);

		for (i = (best_starting * 100); i <= (100 * best_ending); ++i) {
			g.x = (double) i * .01;
			v = evaluate_tree(run_stats[0].best[0]->ind->tr[0].data, 0);
			oprintf( OUT_USER, 50, "%lf %lf\n", g.x, v);
		}

		output_stream_close( OUT_USER);

		if (run_stats[0].best[0]->ind->hits == fitness_cases)
			return 1;
	}

	return 0;
}

void app_end_of_breeding(int gen, multipop *mpop) {
	return;
}

int app_create_output_streams(void) {
	if (create_output_stream( OUT_USER, ".fn", 1, "w", 0) !=
	OUTPUT_OK)
		return 1;

	return 0;
}
bool in_between(double start, double end, double value) {
	if (value > start && value < end) {
		return true;
	}
	return false;
}
int checkImportance(double x) {

	if (in_between(-3, -2.75, x) || in_between(-2.59, -2.2, x)
			|| in_between(-1.39, -0.4, x) || in_between(1.8, 3.0, x)
			|| in_between(3.81, 4.0, x)) {
		return 1;
	} else if (in_between(-2.74, -2.6, x) || in_between(-2.1, -2, x)
			|| in_between(-0.09, 1.7, x) || in_between(3.1, 3.8, x)) {
		return 2;
	} else if (in_between(-1.9, -1.4, x) || in_between(-0.3, -0.1, x)) {
		return 3;
	}
	return 0;
}

int app_initialize(int startfromcheckpoint) {
	int i;
	double x, y;
	char *param;
	init();
	if (!startfromcheckpoint) {
		oprintf( OUT_PRG, 50, "not starting from checkpoint file.\n");

		param = get_parameter("app.fitness_cases");
		if (param == NULL)
			fitness_cases = 200;
		else {
			fitness_cases = atoi(param);
			if (fitness_cases < 0)
				error( E_FATAL_ERROR,
						"invalid value for \"app.fitness_cases\".");
		}
		FILE *in_file = fopen("500_XSquare.csv", "r");
		fscanf(in_file, "%d", &fitness_cases);
		app_fitness_cases[0] = (double *) MALLOC(
				fitness_cases * sizeof(double));
		app_fitness_cases[1] = (double *) MALLOC(
				fitness_cases * sizeof(double));
		app_fitness_importance = (int *) MALLOC(fitness_cases * sizeof(int));
		//Asim Code
		float x, y;
		for (i = 0; i < fitness_cases; ++i) {
			fscanf(in_file, "%f", &x);
			fscanf(in_file, "%f", &y);
			app_fitness_cases[0][i] = x;
			app_fitness_cases[1][i] = y;
			//app_fitness_importance[i] = checkImportance(x);
		}
		/*oprintf( OUT_PRG, 50, "%d fitness cases:\n", fitness_cases);
		 for (i = 0; i < fitness_cases; ++i) {
		 x = (random_double() * 2.0) - 1.0;

		 // change this line to modify the goal function.
		 y = x * x * x * x + x * x * x + x * x + x;

		 app_fitness_cases[0][i] = x;
		 app_fitness_cases[1][i] = y;

		 // oprintf( OUT_PRG, 50, "    x = %12.5lf, y = %12.5lf\n", x, y);
		 }*/
	} else {
		oprintf( OUT_PRG, 50, "started from checkpoint file.\n");
	}

	param = get_parameter("app.value_cutoff");
	if (param == NULL)
		value_cutoff = 1.e15;
	else
		value_cutoff = strtod(param, NULL);

	return 0;
}

void app_uninitialize(void) {

	free(optimal_index_in_generation);
	free(optimal_in_generation);
	free(app_fitness_importance);
	free(app_fitness_cases);
	//free(app_fitness_cases[1]);
	//free(app_fitness_cases[2]);
	//int i = 0;
	//for (; i < generationSIZE; i++) {
	free(error_array);
	//}

}

void app_write_checkpoint(FILE *f) {
	int i;
	fprintf(f, "fitness-cases: %d\n", fitness_cases);
	for (i = 0; i < fitness_cases; ++i) {
		write_hex_block(app_fitness_cases[0] + i, sizeof(double), f);
		fputc(' ', f);
		write_hex_block(app_fitness_cases[1] + i, sizeof(double), f);
		fprintf(f, " %.5lf %.5lf\n", app_fitness_cases[0][i],
				app_fitness_cases[1][i]);
	}
}

void app_read_checkpoint(FILE *f) {
	int i;

	fscanf(f, "%*s %d\n", &fitness_cases);

	app_fitness_cases[0] = (double *) MALLOC(fitness_cases * sizeof(double));
	app_fitness_cases[1] = (double *) MALLOC(fitness_cases * sizeof(double));

	for (i = 0; i < fitness_cases; ++i) {
		read_hex_block(app_fitness_cases[0] + i, sizeof(double), f);
		fgetc(f);
		read_hex_block(app_fitness_cases[1] + i, sizeof(double), f);
		fscanf(f, " %*f %*f\n");
		fprintf( stderr, "%.5lf %.5lf\n", app_fitness_cases[0][i],
				app_fitness_cases[1][i]);
	}
}
