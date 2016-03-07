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

#ifndef _DEFINES_H
#define _DEFINES_H

/* remove this #define to turn off memory tracking. */
#define TRACK_MEMORY

/* define this to write a file "memory.log" with a record of all MALLOC()s,
   REALLOC()s, and FREE()s.  useful for debugging and finding memory leaks. */
/*#define MEMORY_LOG*/

/* when this symbol is defined, the current time will be used instead of
   1 as the default random seed. */
/*#define RANDOMSEEDTIME*/

#define EXTRAMEM              8
#define EPHEM_METABLOCKSIZE   10
#define EPHEM_STARTSIZE       1000
#define EPHEM_GROWSIZE        500

#define MAXPARAMLINELENGTH    255
#define MAXCHECKLINELENGTH    255

#define FUNC_DATA  1
#define FUNC_FUNC  1
#define FUNC_EXPR  2
#define FUNC_MACRO 2
#define TERM_NORM  3
#define TERM_ERC   4
#define TERM_ARG   5
#define EVAL_DATA  6
#define EVAL_FUNC  6
#define EVAL_EXPR  7
#define EVAL_MACRO 7
#define EVAL_TERM  8

#define EVAL_CACHE_INVALID   1
#define EVAL_CACHE_VALID     0

#define SELECT_INIT    1
#define SELECT_CLEAN   3

#define GENERATE_FULL           1
#define GENERATE_GROW           2
#define GENERATE_HALF_AND_HALF  3

#define MAXDETAILLEVEL   100
#define DEFDETAILLEVEL   50
#define MINDETAILLEVEL   1

#define E_WARNING        0
#define E_ERROR          1
#define E_FATAL_ERROR    2

#define OUTPUT_OK        0
#define OUTPUT_DUP_ID    1
#define OUTPUT_DUP_EXT   2
#define OUTPUT_TOOMANY   3
#define OUTPUT_TOOLATE   4

#define OUT_SYS    0
#define OUT_GEN    1
#define OUT_PRG    2
#define OUT_STT    3
#define OUT_BST    4
#define OUT_HIS    5
#define OUT_USER   6
#define OUT_ERROR  7

#define PARAM_COPY_NONE   0
#define PARAM_COPY_NAME   1
#define PARAM_COPY_VALUE  2

#define MAXMESSAGELENGTH 4096
#define MAXOUTPUTSTREAMS 25
#define SYSOUTPUTSTREAMS 6

#define PARAMETER_MINSIZE       31
#define PARAMETER_CHUNKSIZE     16

#define OPERATOR_CROSSOVER      1
#define OPERATOR_REPRODUCE      2
#define OPERATOR_MUTATE         3

#define FLAG_NONE               0
#define FLAG_NEWEXCH            1

#define GENSPACE_COUNT          2

#define GENSPACE_START          100
#define GENSPACE_GROW           100

#define CK_MAGIC                "lilgp1.0\n"
#define CK_IDSTRING             "id: lilgp v1.0 checkpoint file\n"

#endif
