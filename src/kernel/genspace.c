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
/* initialize_genspace()
 *
 * allocates each genspace with GENSPACE_START lnodes.
 */

void initialize_genspace ( void )
{
     int i;

     oputs ( OUT_SYS, 30, "    generation spaces.\n" );
     
     for ( i = 0; i < GENSPACE_COUNT; ++i )
     {
          gensp[i].size = GENSPACE_START;
          gensp[i].data = (lnode *)MALLOC ( gensp[i].size * sizeof ( lnode ) );
          memset ( gensp[i].data, 0, gensp[i].size * sizeof ( lnode ) );
          gensp[i].used = 0;
#ifdef DEBUG
          printf ( "genspace %d initialized with %d nodes.\n",
                  i, gensp[i].size );
#endif
     }
}

/* free_genspace()
 *
 * frees all the genspaces.
 */

void free_genspace ( void )
{
     int i;
     for ( i = 0; i < GENSPACE_COUNT; ++i )
     {
          FREE ( gensp[i].data );
          gensp[i].data = NULL;
     }
}

/* gensp_next()
 *
 * returns the address of the next free lnode in the given generation
 * space.  enlarges the generation space by GENSPACE_GROW lnodes if
 * there is no free space.
 */

lnode * gensp_next ( int space )
{
     while ( gensp[space].used >= gensp[space].size )
     {
          int oldsize = gensp[space].size;
          gensp[space].size += GENSPACE_GROW;
          gensp[space].data = (lnode *)REALLOC ( gensp[space].data,
                                                gensp[space].size *
                                                sizeof ( lnode ) );
          memset ( gensp[space].data+oldsize, 0,
                   (gensp[space].size-oldsize) * sizeof ( lnode ) );
#ifdef DEBUG
          printf ( "next: genspace %d grown to %d nodes.\n",
                  space, gensp[space].size );
#endif
     }

     return gensp[space].data+(gensp[space].used++);
}

/* gensp_next_int()
 *
 * like gensp_next(), but returns the position (not the address) of the
 * next free lnode.
 */

int gensp_next_int ( int space )
{
     while ( gensp[space].used >= gensp[space].size )
     {
          int oldsize = gensp[space].size;
          gensp[space].size += GENSPACE_GROW;
          gensp[space].data = (lnode *)REALLOC ( gensp[space].data,
                                                gensp[space].size *
                                                sizeof ( lnode ) );
          memset ( gensp[space].data+oldsize, 0,
                   (gensp[space].size-oldsize) * sizeof ( lnode ) );
#ifdef DEBUG
          printf ( "next_int: genspace %d grown to %d nodes.\n",
                  space, gensp[space].size );
#endif
     }

     return gensp[space].used++;
}

/* gensp_dup_tree()
 *
 * copies a completed tree out of a generation space into the tree
 * pointer passed.
 */

void gensp_dup_tree ( int space, tree *t )
{
     t->size = gensp[space].used;
     t->nodes = tree_nodes ( gensp[space].data );
     t->data = (lnode *)MALLOC ( t->size * sizeof ( lnode ) );
     memcpy ( t->data, gensp[space].data, t->size * sizeof ( lnode ) );
}

/* gensp_reset()
 *
 * marks a genspace as being empty.
 */

void gensp_reset ( int space )
{
     gensp[space].used = 0;
     memset ( gensp[space].data, 0, gensp[space].size * sizeof ( lnode ) );
}
