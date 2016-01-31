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

#include <lilgp.h>

#ifdef MEMORY_LOG
extern FILE *mlog;
#endif

static int totalalloc = 0;
static int maxalloc   = 0;
static int freealloc  = 0;
static int curalloc   = 0;
static int malloccalls = 0;
static int freecalls = 0;
static int realloccalls = 0;

/* get_memory_stats()
 *
 * returns the memory statistics stored in static global
 * variables.
 */

void get_memory_stats ( int *total, int *free, int *max,
                       int *mallocc, int *reallocc, int *freec )
{
     *total = totalalloc;
     *free = freealloc;
     *max = maxalloc;
     *mallocc = malloccalls;
     *reallocc = realloccalls;
     *freec = freecalls;
}

/* track_malloc()
 *
 * like malloc(), but tracks memory block size.
 */

void *track_malloc ( int size )
{
     unsigned char *p;

     if ( size == 0 )
          return NULL;
     
     p = (unsigned char *)malloc ( size+EXTRAMEM );
     if ( p == NULL )
          return NULL;
     ++malloccalls;
     totalalloc += size;
     curalloc += size;
     if ( curalloc > maxalloc )
          maxalloc = curalloc;
     *(int *)p = size;
#ifdef MEMORY_LOG
     fprintf ( mlog, "MALLOC %d %08x\n", size, (void *)(p+EXTRAMEM) );
     fflush ( mlog );
#endif
     return (void *)(p+EXTRAMEM);
}

/* track_free()
 *
 * like free(), but tracks memory block size.  use with pointers
 * returned by track_malloc().
 */

void track_free ( void *p )
{
     int size;

     if ( p == NULL )
          return;

#ifdef MEMORY_LOG
     fprintf ( mlog, "FREE %08x\n", p );
     fflush ( mlog );
#endif
     
     size = *(int *)((unsigned char *)p-EXTRAMEM);

     ++freecalls;
     curalloc -= size;
     freealloc += size;

     free ( (unsigned char *)p-EXTRAMEM );
}
     
/* *track_realloc()
 *
 * like realloc(), but tracks memory block size.  use with pointers
 * returned by track_malloc().
 */

void *track_realloc ( void *p, int newsize )
{
     int size, change;

     if ( p == NULL )
          return MALLOC ( newsize );
     
     size = *(int *)((unsigned char *)p-EXTRAMEM);

     change = newsize-size;
     curalloc += change;
     if ( curalloc > maxalloc )
          maxalloc = curalloc;
     if ( change > 0 )
          totalalloc += change;
     else
          freealloc -= change;
     ++realloccalls;

#ifdef MEMORY_LOG
     fprintf ( mlog, "REALLOC %08x", p );
#endif
     p = realloc ( (unsigned char *)p-EXTRAMEM, newsize+EXTRAMEM );
     *(int *)p = newsize;
#ifdef MEMORY_LOG
     fprintf ( mlog, " %d %08x\n", newsize, (void *)((unsigned char *)p+EXTRAMEM) );
     fflush ( mlog );
#endif
     return (void *)((unsigned char *)p+EXTRAMEM);
}


