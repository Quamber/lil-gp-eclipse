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

#include <sys/times.h>
#include <limits.h>

#include "event.h"

static int tickspersec;
static char string[200];

void event_init ( void )
{
     tickspersec = CLK_TCK;
}
     
void event_mark ( event *e )
{
     struct tms b;

     e->wall = times(&b);
     e->user = b.tms_utime;
     e->sys = b.tms_stime;
}

char *event_string ( event *z )
{
     sprintf ( string, "%.2lfs wall, %.2lfs user, %.2lfs sys",
              (double)(z->wall)/tickspersec,
              (double)(z->user)/tickspersec,
              (double)(z->sys)/tickspersec );
     return string;
}

void event_zero ( event *z )
{
     z->wall = z->user = z->sys = 0;
}

void event_diff ( event *z, event *one, event *two )
{
     z->wall = two->wall - one->wall;
     z->user = two->user - one->user;
     z->sys = two->sys - one->sys;
}

void event_accum ( event *z, event *one )
{
     z->wall += one->wall;
     z->user += one->user;
     z->sys += one->sys;
}

void event_accumdiff ( event *z, event *one, event *two )
{
     z->wall += two->wall - one->wall;
     z->user += two->user - one->user;
     z->sys += two->sys - one->sys;
}
     
     
