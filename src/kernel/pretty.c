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
/* pretty_print_tree()
 *
 * print a tree as an S-expression, indenting functions to show
 * structure.
 */

void pretty_print_tree ( lnode *data, FILE *fil )
{
     lnode *l = data;
     int *indent;
     int *istep;

     /* indent is a table, parallel to the tree's array, of indentation
	values. */
     indent = (int *)MALLOC ( (tree_nodes(data)+1) * sizeof(int) );
     istep = indent;
     /* fill the indentation table. */
     gen_indents ( &l, &istep, 0, 1 );
     
     istep = indent;
     l = data;

     /* now print the tree using the indentation table. */
     pretty_print_tree_recurse ( &l, &istep, fil );
     
     fprintf ( fil, "\n" );
     
     FREE ( indent );
     
}

/* pretty_print_tree_recurse()
 *
 * recursive tree printer, using a table of indentations.
 */

void pretty_print_tree_recurse ( lnode **l, int **is, FILE *fil )
{
     int i;
     function *f;

     f = (**l).f;

     /** a positive indentation value means move to the next line and
       print that many spaces before the function name. */
     if ( **is >= 0 )
     {
          fprintf ( fil, "\n" );
          for ( i = 0; i < **is; ++i )
               fprintf ( fil, " " );
     }
     ++*l;
     ++*is;

     /** for terminals, don't print parens. **/
     if ( f->arity == 0 )
     {
          if ( f->ephem_gen )
          {
	       /* show value of ERCs. */
               fprintf ( fil, " %s", (f->ephem_str)((**l).d->d) );
               ++*l;
          }
          else
	       /* show name of other terminals. */
               fprintf ( fil, " %s", f->string );
          
          return;
     }

     /* print function name with a parenthesis. */
     fprintf ( fil, " (%s", f->string );
     
     switch ( f->type )
     {
	case FUNC_DATA:
	case EVAL_DATA:
	  /* recursively print children. */
	  for ( i = 0; i < f->arity; ++i )
	       pretty_print_tree_recurse ( l, is, fil );
	  break;
	case FUNC_EXPR:
	case EVAL_EXPR:
	  /* recursively print children, ignoring skip nodes. */
	  for ( i = 0; i < f->arity; ++i )
	  {
	       ++*l;
	       pretty_print_tree_recurse ( l, is, fil );
	  }
	  break;
     }

     /* print matching parenthesis. */
     fprintf ( fil, ")" );
}

/* gen_indents()
 *
 * generates a table of indentations -- positive numbers indicate skipping
 * to next line and indenting that much, negative indicate continuing current
 * line.  the result should look like:
 *
 * (function terminal
 *           (function terminal
 *                     terminal)
 *           terminal)
 */

void gen_indents ( lnode **l, int **is, int start, int sameline )
{
     function *f = (**l).f;
     int i;

     /** sameline is true for the first child of a function.   first
       children and terminals always go on the same line as their
       parent. **/
     if ( sameline || f->arity == 0 )
          **is = -start;
     else
          **is = start;
     
     ++*is;
     ++*l;
     
     if ( f->arity == 0 )
     {
          if ( f->ephem_gen )
	       /* skip the value of an ERC. */
               ++*l;
          return;
     }

     /* move forward the length of the function name plus a space
	plus a '('. */
     start += strlen ( f->string ) + 2;
     
     switch ( f->type )
     {
	case FUNC_DATA:
	case EVAL_DATA:
	  /* generate children's indents. */
	  for ( i = 0; i < f->arity; ++i )
	       gen_indents ( l, is, start, i==0 );
	  break;
	case FUNC_EXPR:
	case EVAL_EXPR:
	  /* generate children's indents, ignoring skip nodes. */
	  for ( i = 0; i < f->arity; ++i )
	  {
	       ++*l;
	       gen_indents ( l, is, start, i==0 );
	  }
	  break;
     }
}

