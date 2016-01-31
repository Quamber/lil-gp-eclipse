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
/*
 * note:  that almost all these functions that recurse through the tree
 * come in two parts:  the wrapper and the recursive part.  the "standard"
 * method of traversing a tree like this is to set up a pointer to the
 * first node in the tree, then subsequent recursive calls step the
 * pointer forward, traversing the tree in preorder.  the wrapper function
 * sets up this pointer, then passes the address of that pointer [that's
 * an (lnode **)] to the recursive part.
 *
 * the first function in this file [tree_size()] is commented in detail,
 * others (which are just variations on this theme) are rather glossed over.
 *
 * note:  when I refer to a tree's "size" I almost always mean "how many
 * lnodes the tree takes to store", NOT how many nodes are in the tree.
 * watch out for this.  (the major exception is, of course, "tree_size()"
 * which returns the node count.)
 */

/*
 * tree_nodes:  return the number of nodes in the tree.
 */

int tree_nodes ( lnode *tree )
{
     lnode *l = tree;
     return tree_nodes_recurse ( &l );
}

int tree_nodes_recurse ( lnode **l )
{
     /*
      *  *l always points at a function node here; save the function.
      */
     function *f = (**l).f;
     int i, j = 1;

     /* step the pointer over the function. */
     ++*l;

     /* if the function is a terminal, then the recursion bottoms out. */
     if ( f->arity == 0 )
     {
          /* skip the pointer over the ERC value if this node has one. */
          if ( f->ephem_gen )
               ++*l;

          /* return this subtree size. */
          return 1;
     }
     else
     {
          /* function is not a terminal, so add up its subtrees and return
           * the total.
           */

          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
                    j += tree_nodes_recurse ( l );
               break;
               
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    /* skip the pointer over the skipsize node. */
                    ++*l;
                    /* add this subtree's size to the total. */
                    j += tree_nodes_recurse ( l );
               }
               break;
          }
     }
     
     return j;
}

/*
 * tree_nodes_internal:  return number of internal (function) nodes in the
 *     tree.
 */

int tree_nodes_internal ( lnode *data )
{
     lnode *l = data;
     return tree_nodes_internal_recurse ( &l );
}

int tree_nodes_internal_recurse ( lnode **l )
{
     function *f = (**l).f;
     int i, j = 1;

     ++*l;
     
     if ( f->arity == 0 )
     {
          if ( f->ephem_gen )
               ++*l;
          return 0;
     }
     else
     {
          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
                    j += tree_nodes_internal_recurse ( l );
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    ++*l;
                    j += tree_nodes_internal_recurse ( l );
               }
               break;
          }
     }

     return j;
}

/*
 * tree_nodes_external:  return number of external (terminal) nodes in the
 *     tree.
 */

int tree_nodes_external ( lnode *data )
{
     lnode *l = data;
     return tree_nodes_external_recurse ( &l );
}

int tree_nodes_external_recurse ( lnode **l )
{
     function *f = (**l).f;
     int i, j = 0;

     ++*l;
     
     if ( f->arity == 0 )
     {
          if ( f->ephem_gen )
               ++*l;
          return 1;
     }
     else
     {
          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
                    j += tree_nodes_external_recurse ( l ) ;
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    ++*l;
                    j += tree_nodes_external_recurse ( l ) ;
               }
               break;
          }
     }
     return j;
}

/*
 * print_tree_array:  print the tree, as an array, on stdout.  primarily
 *     for debugging.
 */

void print_tree_array ( lnode *data )
{
     lnode *l = data;
     int i = 0;
     print_tree_array_recurse ( &l, &i );
}

void print_tree_array_recurse ( lnode **l, int *index )
{
     function *f = (**l).f;
     int i;

     printf ( "%3d: function: \"%s\"\n", *index, f->string );
     
     ++*l;
     ++*index;

     if ( f->arity == 0 )
     {
          if ( f->ephem_gen )
          {
               printf ( "%3d:    value: %s\n", *index, (f->ephem_str)((**l).d->d) );
               ++*l;
               ++*index;
          }
     }
     else
     {
          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
               {
                    print_tree_array_recurse ( l, index );
               }
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    printf ( "%3d:    {skip %d}\n", *index, (**l).s );
                    ++*l;
                    ++*index;
                    print_tree_array_recurse ( l, index );
               }
               break;
          }
     }
}

/*
 * print_tree:  print the tree, as a LISP S-expression, to the given FILE *.
 */

void print_tree ( lnode *data, FILE *fil )
{
     lnode *c = data;
     print_tree_recurse ( &c, fil );
     fprintf ( fil, "\n" );
}

void print_tree_recurse ( lnode **l, FILE *fil )
{
     function *f;
     int i;

     f = (**l).f;

     fprintf ( fil, " " );
     if ( f->arity != 0 )
          fprintf ( fil, "(" );
     
     ++*l;
     if ( f->ephem_gen )
     {
          fprintf ( fil, "%s", (f->ephem_str)((**l).d->d) );
#ifdef DEBUG
          fprintf ( fil, " <%d>", (**l).d->refcount );
#endif
          ++*l;
     }
     else
          fprintf ( fil, "%s", f->string );
     
     switch ( f->type )
     {
        case FUNC_DATA:
        case EVAL_DATA:
          for ( i = 0; i < f->arity; ++i )
          {
               print_tree_recurse ( l, fil );
          }
          break;
        case FUNC_EXPR:
        case EVAL_EXPR:
          for ( i = 0; i < f->arity; ++i )
          {
#ifdef DEBUG               
               fprintf ( fil, " {skip %d}", (**l).s );
#endif               
               ++*l;
               print_tree_recurse ( l, fil );
          }
          break;
     }

     if ( f->arity != 0 )
          fprintf ( fil, ")" );
}

/*
 * generate_random_full_tree:  generates, in the given generation space,
 *     a random full tree of
 *     the specified depth. 
 */

void generate_random_full_tree ( int space, int depth, function_set *fset )
{
     int i, j;
     int save;

     if ( depth == 0 )
     {
          /* we've reached depth 0, so choose a terminal node for this
           * subtree.
           */
          
          i = random_int(fset->terminal_count)+(fset->function_count);
          gensp_next(space)->f = (fset->cset)+i;

          /* if this terminal is an ERC, then generate one and store it. */
          if ( (fset->cset)[i].ephem_gen )
               gensp_next(space)->d = new_ephemeral_const ( (fset->cset)+i );

          return;
     }

     /* we're not a depth 0, so generate a non-terminal node. */
     i = random_int(fset->function_count);
     gensp_next(space)->f = (fset->cset)+i;

     /* now generate the node's children. */
          switch ( (fset->cset)[i].type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( j = 0; j < (fset->cset)[i].arity; ++j )
               {
                    /* generate a full tree of depth [depth-1]. */
                    generate_random_full_tree ( space, depth-1, fset );
               }
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( j = 0; j < (fset->cset)[i].arity; ++j )
               {
                    /* yes, so save where the skip value needs to be placed
                     * since we don't know what it is yet.
                     */
                    save = gensp_next_int ( space );
                    
                    /* generate a full tree of depth [depth-1]. */
                    generate_random_full_tree ( space, depth-1, fset );
                    
                    /* save the size of the subtree in the skip node, if we need to. */
                    gensp[space].data[save].s = gensp[space].used-save-1;
                    
               }
               break;
          }
     
     return;
     
}

/*
 * generate_random_grow_tree:  grow a random tree, of (maximum) depth [depth].
 *     works just like genereate_random_full_tree, except that a node with
 *     depth > 0 can be a terminal.
 */
       
void generate_random_grow_tree ( int space, int depth, function_set *fset )
{
     int i, j;
     int save;

     if ( depth == 0 )
     {
          i = random_int(fset->terminal_count)+(fset->function_count);
          gensp_next(space)->f = (fset->cset)+i;

          if ( (fset->cset)[i].ephem_gen )
               gensp_next(space)->d = new_ephemeral_const ( (fset->cset)+i );
          
          return;
     }

     /* note that here we can generate a terminal node. */
     i = random_int(fset->function_count+fset->terminal_count);
     gensp_next(space)->f = (fset->cset)+i;
     
     if ( (fset->cset)[i].arity )
     {
          switch ( (fset->cset)[i].type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( j = 0; j < (fset->cset)[i].arity; ++j )
                    generate_random_grow_tree ( space, depth-1, fset );
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( j = 0; j < (fset->cset)[i].arity; ++j )
               {
                    save = gensp_next_int(space);
                    generate_random_grow_tree ( space, depth-1, fset );
                    gensp[space].data[save].s = gensp[space].used-save-1;
               }
               break;
          }
     }
     else
     {
          if ( (fset->cset)[i].ephem_gen )
               gensp_next(space)->d = new_ephemeral_const ( (fset->cset)+i );
     }
     return;
}

/*
 * tree_depth:  return the depth of the tree.  a tree with one node has
 *     depth 0.
 */

int tree_depth ( lnode *data )
{
     lnode *l = data;
     return tree_depth_recurse ( &l );
}

int tree_depth_recurse ( lnode **l )
{

     function *f = (**l).f;
     int i, j, k = 0;

     ++*l;
     if ( f->arity == 0 )
     {
          /* a terminal node; advance the pointer and return 0. */
          
          if ( f->ephem_gen )
               ++*l;
          return 0;
     }

     /* a nonterminal; find the deepest child and return its depth plus one. */
     
          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
               {
                    j = tree_depth_recurse ( l );
                    if ( j > k )
                         k = j;
               }
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    ++*l;
                    j = tree_depth_recurse ( l );
                    if ( j > k )
                         k = j;
               }
               break;
          }
     
     return k+1;
}

int tree_depth_to_subtree ( lnode *data, lnode *sub )
{
     lnode *l = data;
     return tree_depth_to_subtree_recurse ( &l, sub, 0 );
}

int tree_depth_to_subtree_recurse ( lnode **l, lnode *sub, int depth )
{
     function *f = (**l).f;
     int i, j;

     if ( *l == sub )
          return depth;

     ++*l;
     if ( f->arity == 0 )
     {
          if ( f->ephem_gen )
               ++*l;
          return -1;
     }

          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
               {
                    j = tree_depth_to_subtree_recurse ( l, sub, depth+1 );
                    if ( j != -1 )
                         return j;
               }
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    ++*l;
                    j = tree_depth_to_subtree_recurse ( l, sub, depth+1 );
                    if ( j != -1 )
                         return j;
               }
               break;
          }

     return -1;
}
     

/*
 * get_subtree:  returns the start'th subtree of the tree, where start
 *     ranges from 0..(nodecount-1).  the zeroth subtree is the tree itself.
 *     the subtrees are returned in preorder.
 */

lnode *get_subtree ( lnode *data, int start )
{
     lnode *l = data;
     /* initialize a counter with start.  the counter is passed by address,
      * and so is shared across all recursive calls.  it is decremented once
      * for every node in the tree -- the current node when it reaches zero
      * is returned.
      */
     int c = start;
     
     return get_subtree_recurse ( &l, &c );
}

lnode *get_subtree_recurse ( lnode **l, int *c )
{
     function *f = (**l).f;
     lnode *r;
     int i;

     /* if the counter is zero, return the current subtree. */
     if ( *c == 0 )
          return *l;

     ++*l;
     --*c;

     if ( f->arity == 0 )
     {
          /* skip over the terminal nodes. */
          if ( f->ephem_gen )
               ++*l;
     }
     else
     {
          /* recurse into this node's children.  if one of them returns
           * non-NULL, it means that the subtree has been found and the
           * return value is immediately propagated up the call chain.
           * if all of them return NULL, then the desired subtree is
           * not in this subtree and we return NULL. */
          
          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
               {
                    r = get_subtree_recurse ( l, c );
                    if ( r != NULL )
                         return r;
               }
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    ++*l;
                    r = get_subtree_recurse ( l, c );
                    if ( r != NULL )
                         return r;
               }
               break;
          }
     }

     return NULL;
}

/*
 * get_subtree_internal:  just like get_subtree, but only selects nonterminal
 *     points.  start should range from 0..(internalnodecount-1).
 */

lnode *get_subtree_internal ( lnode *data, int start )
{
     lnode *l = data;
     int c = start;
     return get_subtree_internal_recurse ( &l, &c );
}

lnode *get_subtree_internal_recurse ( lnode **l, int *c )
{
     function *f = (**l).f;
     int i;
     lnode *r;

     /* if the current subtree is a terminal node, skip it immediately. */
     if ( f->arity == 0 )
     {
          if ( f->ephem_gen )
               ++*l;
          ++*l;
     }
     else
     {
          if ( *c == 0 )
               return *l;

          ++*l;
          --*c;

          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
               {
                    r = get_subtree_internal_recurse ( l, c );
                    if ( r != NULL )
                         return r;
               }
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    ++*l;
                    r = get_subtree_internal_recurse ( l, c );
                    if ( r != NULL )
                         return r;
               }
               break;
          }
     }

     return NULL;
}

/*
 * get_subtree_external:  just like get_subtree, but only selects terminal
 *     nodes (that is, 1-node subtrees).  start ranges from
 *     0..(externalnodecount-1).
 */

lnode *get_subtree_external ( lnode *data, int start )
{
     lnode *l = data;
     int c = start;
     return get_subtree_external_recurse ( &l, &c );
}

lnode *get_subtree_external_recurse ( lnode **l, int *c )
{
     function *f = (**l).f;
     int i;
     lnode *r;

     if ( f->arity == 0 )
     {
          if ( *c == 0 )
               return *l;

          ++*l;
          --*c;

          if ( f->ephem_gen )
               ++*l;
     }
     else
     {
          ++*l;

          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
               {
                    r = get_subtree_external_recurse ( l, c );
                    if ( r != NULL )
                         return r;
               }
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    ++*l;
                    r = get_subtree_external_recurse ( l, c );
                    if ( r != NULL )
                         return r;
               }
               break;
          }
     }

     return NULL;
}

/*
 * copy_tree:  allocates space for and makes a copy of a tree.
 */

void copy_tree ( tree *to, tree *from )
{
     to->data = (lnode *)MALLOC ( from->size * sizeof ( lnode ) );
     to->size = from->size;
     to->nodes = from->nodes;
     memcpy ( to->data, from->data, from->size * sizeof ( lnode ) );
}

/*
 * free_tree:  frees the memory allocated by a tree, and resets variables.
 */

void free_tree ( tree *t )
{
     FREE ( t->data );
     t->data = NULL;
     t->size = -1;
     t->nodes = -1;
}

/*
 * tree_size:  returns the number of lnodes used to store the tree.
 *     this is equal to
 *        (node count)+(ERC count)+(conditionally evaluated subtree count).
 */

int tree_size ( lnode *data )
{
     lnode *l = data;
     return tree_size_recurse ( &l );
}

int tree_size_recurse ( lnode **l )
{
     function *f = (**l).f;
     int i, j = 1;

     ++*l;

     if ( f->arity == 0 )
     {
          if ( f->ephem_gen )
          {
               ++*l;
               ++j;
          }
     }
     else
     {
          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
                    j += tree_size_recurse ( l );
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    ++*l;
                    ++j;
                    j += tree_size_recurse ( l );
               }
               break;
          }
     }

     return j;

}

/*
 * copy_tree_replace_many:  copies a tree, replacing some of its subtrees
 *     with other subtrees.  arguments:
 *
 *         space      - the generation space to put the new tree in
 *         parent     - the tree to copy from
 *         replace    - a list of subtrees to replace
 *         with       - a list of subtrees to replace them with
 *         count      - the size of the replace and width arrays
 *         repcount   - returns the number of subtrees replaced
 *
 *     this function starts to recursively copy the "parent" tree to "dest".
 *     when the recursive copy hits any subtree in the "replace" array,
 *     the corresponding subtree in the "with" array is copied in its
 *     place.  repcount returns the number of subtrees that were hit and
 *     replaced.  this can be less than count, as some subtrees may never
 *     be found (if, for instance, one of the subtrees in the "replace"
 *     array is the subtree of another).
 */

void copy_tree_replace_many ( int space, lnode *parent, lnode **replace,
                            lnode **with, int count, int *repcount )
{
     lnode *lp = parent;

     gensp_reset ( space );
     *repcount = 0;
     copy_tree_replace_many_recurse ( space, &lp, replace,
                                     with, count, repcount );
}

void copy_tree_replace_many_recurse ( int space, lnode **lp, lnode **lr,
                                    lnode **lw, int count, int *repcount )
{
     function *f = (**lp).f;
     int i;
     int save;
     lnode *new;

     /* only do the comparison if the lr is non-NULL. */
     if ( lr )
     {
          /* check the current subtree against everything in the lr array. */
          for ( i = 0; i < count; ++i )
               if ( *lp == lr[i] )
               {
                    /* we have a match! */

                    /* increment the replacement count. */
                    ++*repcount;
                    
                    /* copy the new tree into the destination.  note that the
                     * lr and lw arguments are passed as NULL to prevent
                     * further replacement within the replacement tree.
                     */
                    new = lw[i];
                    copy_tree_replace_many_recurse ( space, &new, NULL,
                                                    NULL, count, repcount );
                    
                    /* now skip the lp pointer over the lr subtree. */
                    skip_over_subtree ( lp );

                    return;
               }
     }

     /* copy the node from the parent to the destination. */

     gensp_next(space)->f = (**lp).f;
     ++*lp;

     if ( f->arity == 0 )
     {
          if ( f->ephem_gen )
          {
               /* copy the ERC pointer. */

               gensp_next(space)->d = (**lp).d;
               ++*lp;
          }
     }
     else
     {
          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
               {
                    copy_tree_replace_many_recurse ( space, lp, lr, lw,
                                                    count, repcount );
               }
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    /* note that we just can't copy the skip value, since
                     * replacement within the subtree may change it.  we
                     * have to save where it is needed and fill it in 
                     * afterwards, just like generate_random_full_tree(). */
                    
                    save = gensp_next_int ( space );
                    ++*lp;
                    copy_tree_replace_many_recurse ( space, lp, lr, lw,
                                                    count, repcount );
                    gensp[space].data[save].s = gensp[space].used-save-1;
               }
               break;
          }
     }
     return;
}

/*
 * skip_over_subtree:  takes a traversal pointer and skips it over the
 *     subtree it points to.
 */

void skip_over_subtree ( lnode **l )
{
     function *f = (**l).f;
     int i;

     ++*l;
     if ( f->arity == 0 )
     {
          if ( f->ephem_gen )
               ++*l;
     }
     else
     {
          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
                    skip_over_subtree ( l );
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    ++*l;
                    skip_over_subtree ( l );
               }
               break;
          }
     }
     return;
}

/*
 * reference_ephem_constants:  traverses a tree, adding "count" to the
 *     refcount of each ERC that is used within the tree.
 */

void reference_ephem_constants ( lnode *data, int count )
{
     lnode *l = data;
     reference_ephem_constants_recurse ( &l, count );
}

void reference_ephem_constants_recurse ( lnode **l, int count )
{
     function *f = (**l).f;
     int i;

     ++*l;

     if ( f->arity == 0 )
     {
          if ( f->ephem_gen )
          {
               if ( (**l).d )
               {
                    (**l).d->refcount += count;
                    ++*l;
               }
               else
                    error ( E_FATAL_ERROR, "aarg: this can't happen." );
          }
     }
     else
     {
          switch ( f->type )
          {
             case FUNC_DATA:
             case EVAL_DATA:
               for ( i = 0; i < f->arity; ++i )
                    reference_ephem_constants_recurse ( l, count );
               break;
             case FUNC_EXPR:
             case EVAL_EXPR:
               for ( i = 0; i < f->arity; ++i )
               {
                    ++*l;
                    reference_ephem_constants_recurse ( l, count );
               }
               break;
          }
     }

}


