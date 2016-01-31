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
/* set_current_individual()
 *
 * the current_individual variable is used so that evaluation tokens know
 * where to find their target trees.  set_current_individual() is used to
 * set this value from the application code.
 */

static individual * current_individual;

void set_current_individual ( individual *ind )
{
     current_individual = ind;
}

/* evaluate_tree()
 *
 * this is the wrapper which sets up a traversal pointer for doing the
 * evaluation.  the whichtree argument is needed so that ARG terminals
 * know where to find their values.
 */

DATATYPE evaluate_tree ( lnode *tree, int whichtree )
{
     lnode *l = tree;
#ifdef DEBUG_EVAL
     printf ( "call to evaluate_tree in context %d\n", whichtree );
#endif
     return evaluate_tree_recurse ( &l, whichtree );
}

/* evaluate_tree_recurse()
 *
 * the recursive part of the tree evaluator.  returns the value of the
 * evaluated tree.
 */

DATATYPE evaluate_tree_recurse ( lnode **l, int whichtree )
{

     farg arg[MAXARGS];
     int i;
     function *f = (**l).f;
     treeinfo savearg;

     /* step the traversal pointer forward, now that we've saved which
	function we're at. */
     ++*l;
     
     switch ( f->type )
     {
        case TERM_NORM:
	  /* normal terminal:  just call the user code. */
          return (f->code)(whichtree, NULL);
          break;
        case TERM_ERC:
	  /* ERC terminal:  traversal pointer points to ERC structure.
	     pull the value out, and step the pointer forward. */
          return (*((*l)++)).d->d;
          break;
        case FUNC_DATA:
	  /* function (DATA type):  recursively evaluate each subtree,
	     saving the values.  pass these values to the user code. */
          for ( i = 0; i < f->arity; ++i )
               arg[i].d = evaluate_tree_recurse ( l, whichtree );
          return (f->code)(whichtree, arg);
          break;
        case EVAL_DATA:
	  /* evaluation token (DATA type):  first evaluate each child,
	     saving the returned values.  these values will be returned
	     by the appropriate ARG tokens in the called tree. */
          for ( i = 0; i < f->arity; ++i )
               arg[i].d = evaluate_tree_recurse ( l, whichtree );

	  /* the arguments are stored using three fields in the global
	     tree_map.  we save whatever was in these three fields in
	     local variables. */
          savearg.arguments = tree_map[f->evaltree].arguments;
          savearg.argtype = tree_map[f->evaltree].argtype;
          savearg.evaluatedfrom = tree_map[f->evaltree].evaluatedfrom;
	       
	  /** now we store the new values in the global structure. **/
	  /* first, the argument list. */
          tree_map[f->evaltree].arguments = arg;
	  /* next, the type of this eval token. */
          tree_map[f->evaltree].argtype = EVAL_DATA;
	  /* now the tree number which we are currently evaluating.  this
	     is necessary so that nested ADFs evaluate correctly -- we must
	     remember where we are so that ARG tokens know which tree's
	     arguments to look at. */
          tree_map[f->evaltree].evaluatedfrom = whichtree;
	  
	  /* finally call evaluate_tree to evaluate the target tree. */
          arg->d = evaluate_tree ( current_individual->tr[f->evaltree].data,
                                  f->evaltree );
	  
	  /* restore the old values in the global structure. */
          tree_map[f->evaltree].arguments = savearg.arguments;
          tree_map[f->evaltree].argtype = savearg.argtype;
          tree_map[f->evaltree].evaluatedfrom = savearg.evaluatedfrom;
	  
	  /* return the final value. */
          return arg->d;
          break;
        case FUNC_EXPR:
	  /* function (EXPR type):  save the address of each child tree,
	     using the skip nodes to quickly move the traversal pointer
	     past the child. */
          for ( i = 0; i < f->arity; ++i )
          {
               arg[i].t = (*l+1);
               *l += (**l).s;
               ++*l;
          }
	  /* now pass this array of saved trees to the user code. */
          return (f->code)(whichtree, arg);
          break;
        case EVAL_EXPR:
	  /* evaluation token (EXPR type):  works just like the DATA
	     type, only the saved arguments are tree address instead
	     of values. */
          for ( i = 0; i < f->arity; ++i )
          {
               arg[i].t = (*l+1);
               *l += (**l).s;
               ++*l;
          }
          savearg.arguments = tree_map[f->evaltree].arguments;
          savearg.argtype = tree_map[f->evaltree].argtype;
          savearg.evaluatedfrom = tree_map[f->evaltree].evaluatedfrom;
          tree_map[f->evaltree].arguments = arg;
          tree_map[f->evaltree].argtype = EVAL_EXPR;
          tree_map[f->evaltree].evaluatedfrom = whichtree;
          arg->d = evaluate_tree ( current_individual->tr[f->evaltree].data,
                             f->evaltree );
          tree_map[f->evaltree].arguments = savearg.arguments;
          tree_map[f->evaltree].argtype = savearg.argtype;
          tree_map[f->evaltree].evaluatedfrom = savearg.evaluatedfrom;
          return arg->d;
          break;
        case EVAL_TERM:
	  /* evaluation token (TERM type):  works just like the DATA
	     type, only we pass NULL as an argument list. */
          savearg.arguments = tree_map[f->evaltree].arguments;
          savearg.argtype = tree_map[f->evaltree].argtype;
          savearg.evaluatedfrom = tree_map[f->evaltree].evaluatedfrom;
          tree_map[f->evaltree].arguments = NULL;
          tree_map[f->evaltree].argtype = EVAL_TERM;
          tree_map[f->evaltree].evaluatedfrom = whichtree;
          arg->d = evaluate_tree ( current_individual->tr[f->evaltree].data,
                             f->evaltree );
          tree_map[f->evaltree].arguments = savearg.arguments;
          tree_map[f->evaltree].argtype = savearg.argtype;
          tree_map[f->evaltree].evaluatedfrom = savearg.evaluatedfrom;
          return arg->d;
          break;
        case TERM_ARG:
	  /* an ARG terminal. */
          if ( tree_map[whichtree].argtype == EVAL_DATA )
	       /* if the EVAL token calling this tree is of type DATA, then
		  just pull the value out of the argument list and return it. */
               return tree_map[whichtree].arguments[f->evaltree].d;
          else
	       /* if the EVAL token calling this tree is of type EXPR, then
		  evaluate the tree pointer in the argument list and return
		  the value. */
               return evaluate_tree ( tree_map[whichtree].arguments[f->evaltree].t, tree_map[whichtree].evaluatedfrom );
          break;
     }

}



