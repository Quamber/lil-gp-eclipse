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
/* total counts of ERCs used and freed */
int ercused = 0;
int ercfree = 0;
int ercalloc = 0;

/* the active list. */
ephem_const *active_head;
int active_count;

/* the free list. */
ephem_const *free_head;
int free_count;

/* pointers for the allocated blocks for ERCs. */
ephem_const **block_list;
int block_list_size;
int block_count;

/* initialize_ephem_const()
 *
 * allocate and set up the first block of ERCs, the free list,
 * etc.
 */

void initialize_ephem_const ( void )
{
     int i;
     int size;

     oputs ( OUT_SYS, 30, "    ephemeral random constants.\n" );

     active_head = (ephem_const *)MALLOC ( sizeof ( ephem_const ) );
     active_head->refcount = 1;
     active_head->next = NULL;
     
     free_head = (ephem_const *)MALLOC ( sizeof ( ephem_const ) );
     free_head->refcount = 1;
     free_head->next = NULL;
     
     /** how many block >>pointers<< to allocate (not blocks). **/
     block_list_size = EPHEM_METABLOCKSIZE;
     block_list = (ephem_const **)MALLOC ( block_list_size *
					  sizeof ( ephem_const * ) );

     /* allocate the first block. */
     size = EPHEM_STARTSIZE;
     block_count = 1;
     block_list[0] = (ephem_const *)MALLOC ( size *
                                                  sizeof ( ephem_const ) );
     free_count = ercalloc = size;

     /* chain all the ERC records in the block together. */
     for ( i = 0; i < size-1; ++i )
          block_list[0][i].next = block_list[0]+i+1;
     block_list[0][size-1].next = NULL;

     /* add the chain to the free list. */
     free_head->next = block_list[0];
}

/* free_ephem_const()
 *
 * free all the memory allocated to hold ERCs.
 */

void free_ephem_const ( void )
{
     int i;

     ephem_const_gc();
     
     for ( i = 0; i < block_count; ++i )
          FREE ( block_list[i] );
     FREE ( block_list );

     FREE ( active_head );
     FREE ( free_head );
}

/* enlarge_ephem_space()
 *
 * allocate a new block of ERCs, and add all the records in it
 * to the free list.
 */

void enlarge_ephem_space ( void )
{
     int i;
     ephem_const *p;
     int size = EPHEM_GROWSIZE;

     /** if we've allocated too many blocks, we need to lengthen the
       list that holds block pointers. **/
     if ( block_count == block_list_size )
     {
          block_list_size += EPHEM_METABLOCKSIZE;
          block_list = (ephem_const **)REALLOC ( block_list,
						block_list_size *
						sizeof ( ephem_const *));
     }

     /* allocate the new block. */
     block_list[block_count] =
          (ephem_const *)MALLOC ( size * sizeof ( ephem_const ) );
     free_count += size;
     ercalloc += size;

     /* chain together all the records in it. */
     for ( i = 0; i < size-1; ++i )
          block_list[block_count][i].next =
               block_list[block_count]+i+1;
     
     block_list[block_count][size-1].next = free_head->next;
     free_head->next = block_list[block_count];
     
     ++block_count;
}
     
/* new_ephemeral_const()
 *
 * create a new ERC, corresponding to the given function.
 */

ephem_const *new_ephemeral_const ( function *f )
{
     ephem_const *p;

     /* make sure we have enough space. */
     while ( free_count <= 0 )
          enlarge_ephem_space();

     /* take the next record off the free list. */
     p = free_head->next;
     free_head->next = free_head->next->next;
     --free_count;

     /* call user code to generate the constant, placing
	the value in the new record. */
     f->ephem_gen ( &(p->d) );
     p->f = f;

     /* no references yet. */
     p->refcount = 0;

     /* add this record to the linked list of active ERCs. */
     p->next = active_head->next;
     active_head->next = p;
     ++active_count;

     ++ercused;
     
     return p;
}
     
/* ephem_const_gc()
 *
 * traverse the linked list of active ERCs, removing those with
 * a reference count of 0.
 */

void ephem_const_gc ( void )
{
     ephem_const *p = active_head->next;
     ephem_const *m = active_head;

     while ( p != NULL )
     {
          if ( p->refcount == 0 )
          {
	       /* patch the linked list to skip this record. */
               m->next = p->next;

	       /* add this record to the free list. */
               p->next = free_head->next;
	       free_head->next = p;
               ++free_count;
               ++ercfree;
	       --active_count;
	       
               p = m->next;
          }
          else
          {
	       /* move past this record. */
               m = p;
               p = p->next;
          }
     }

}
               
/* read_ephem_list()
 *
 * read list of ERCs from a checkpoint file. 
 */

ephem_const **read_ephem_list ( FILE *f )
{
     ephem_const **ind;
     ephem_const *p;
     int count;
     int i, j;
     char *buffer;

     /* read the count. */
     fscanf ( f, "%*s %d\n", &count );

     /** if no ERCs, do nothing and return a NULL pointer for the index. **/
     if ( count == 0 )
	  return NULL;

     /* somewhere to place (and ignore) the human-readable forms of
	the ERCs after the hex blocks. */
     buffer = (char *)MALLOC ( MAXCHECKLINELENGTH );

     /* allocate the index translating integers --> addresses. */
     ind = (ephem_const **)MALLOC ( count * sizeof ( ephem_const * ) );
     
     /** allocate a new block to hold all the ERCs we read. **/

     /* we shouldn't EVER need to lengthen the block list while reading
	a checkpoint, but check anyway... */
     if ( block_count == block_list_size )
     {
          block_list_size += EPHEM_METABLOCKSIZE;
          block_list = (ephem_const **)REALLOC ( block_list,
						block_list_size *
						sizeof ( ephem_const *));
     }

     /* allocate the new block. */
     block_list[block_count] =
	  (ephem_const *)MALLOC ( count * sizeof ( ephem_const ) );
     ercalloc += count;
     
     /* read the checkpointed ERCs into the new block. */
     for ( i = 0; i < count; ++i )
     {
	  p = block_list[block_count]+i;
	  fscanf ( f, "%d %d ", &j, &(p->refcount) );
	  ind[j] = p;
	  read_hex_block ( &(p->d), sizeof ( DATATYPE ), f );
	  fgets ( buffer, MAXCHECKLINELENGTH, f );

	  /* chain together all the records. */
	  block_list[block_count][i].next = block_list[block_count]+i+1;
     }
     /* add the chained block to the active list. */
     block_list[block_count][count-1].next = active_head->next;
     active_head->next = block_list[block_count];
     active_count += count;
     ercused += count;
     
     ++block_count;

     FREE ( buffer );
     
     return ind;
     
}
     
/* write_ephem_list()
 *
 * write the active list of ERCs to a checkpoint file.  returns an
 * index listing for translating an ERC address to a unique integer.
 * (we can't store the address directly in the checkpoint, since
 * the ERCs won't be loaded in the same spot in memory on restart.)
 */

ephem_index *write_ephem_list ( FILE *f )
{
     ephem_index *ind;
     ephem_const *p = active_head->next;
     int j;

     ind = (ephem_index *)MALLOC ( active_count * sizeof ( ephem_index ) );
     fprintf ( f, "erc-count: %d\n", active_count );
     
     j = 0;
     while ( p )
     {
	  /* store the index entry. */
          ind[j].e = p;
          ind[j].i = j;

	  /* write the reference count and the value. */
	  fprintf ( f, "%d %d ", j, p->refcount );
	  write_hex_block ( &(p->d), sizeof(DATATYPE), f );
	  fprintf ( f, " %s %s\n", p->f->string, p->f->ephem_str ( p->d ) );

	  /* move down the list. */
          p = p->next;
          ++j;
     }

     /* sort the index by address, so we can use binary searching
	on it. */
     qsort ( ind, active_count, sizeof(ephem_index), ephem_index_comp );

     return ind;
}

/* lookup_ephem()
 *
 * look up an ERC (by address) in an index returned by write_ephem_list()
 * and return its integer index.
 */

int lookup_ephem ( ephem_index *ind, ephem_const *e )
{
     int low = 0;
     int high = (ercused-ercfree);
     int mid;

     while ( low < high-1 )
     {
          mid = (low+high)/2;
          if ( e >= ind[mid].e )
               low = mid;
          else
               high = mid;
     }

     return ind[low].i;
}
     
/* ephem_index_comp()
 *
 * comparison function for using qsort() to order the ERC index by
 * address.
 */

int ephem_index_comp ( const void *a, const void *b )
{
     if ( ((ephem_index *)a)->e < ((ephem_index *)b)->e )
          return -1;
     else
          return 1;
}

/* get_ephem_stats()
 *
 * return ERC statistics.
 */

void get_ephem_stats ( int *used, int *free, int *blocks, int *alloc )
{
     *used = ercused;
     *free = ercfree;
     *blocks = block_count;
     *alloc = ercalloc;
}
