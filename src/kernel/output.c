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
char buffer[MAXMESSAGELENGTH];

typedef struct
{
     int id;
     char *ext;
     int reset;
     char *mode;
     int autoflush;
     
     FILE *f;
     int valid;

     char *buffer;
} output;

output streams[MAXOUTPUTSTREAMS] =
     { { OUT_SYS, ".sys", 0, "w", 1, NULL, 0 },
       { OUT_GEN, ".gen", 0, "w", 0, NULL, 0 },
       { OUT_PRG, ".prg", 0, "w", 0, NULL, 0 },
       { OUT_STT, ".stt", 0, "w", 0, NULL, 0 },
       { OUT_BST, ".bst", 1, "w", 0, NULL, 0 },
       { OUT_HIS, ".his", 0, "w", 0, NULL, 0 } };
     
int output_stream_count = SYSOUTPUTSTREAMS;
int toolate = 0;
char *global_basename = NULL;
int detail_level = DEFDETAILLEVEL;

char error_type[3][20] = { "WARNING: ", "ERROR: ", "FATAL ERROR: " };

extern int quietmode;

/* create_output_stream()
 *
 * make a new entry in the outputstream table.
 */

int create_output_stream ( int id, char *ext, int reset, char *mode,
                          int autoflush )
{
     int i;

     /* can't create a new stream if they've already been opened. */
     if ( toolate )
          return OUTPUT_TOOLATE;

     /* is the stream table full? */
     if ( output_stream_count >= MAXOUTPUTSTREAMS )
          return OUTPUT_TOOMANY;

     /* is the id or extension the same as an existing stream? */
     for ( i = 0; i < output_stream_count; ++i )
     {
          if ( id == streams[i].id )
               return OUTPUT_DUP_ID;
          if ( strcmp ( ext, streams[i].ext ) == 0 )
               return OUTPUT_DUP_EXT;
     }

     /** save stuff in the table. **/
     
     streams[output_stream_count].id = id;
     streams[output_stream_count].ext = (char *)malloc ( strlen(ext)+1 );
     strcpy ( streams[output_stream_count].ext, ext );
     streams[output_stream_count].reset = reset;
     streams[output_stream_count].mode = (char *)malloc ( strlen(mode)+1 );
     strcpy ( streams[output_stream_count].mode, mode );
     streams[output_stream_count].autoflush = autoflush;
     
     streams[output_stream_count].f = NULL;
     streams[output_stream_count].valid = 0;
     ++output_stream_count;

     return OUTPUT_OK;
     
}

/* initialize_output_streams()
 *
 * initialize the output streams.  until they are opened, anything printed
 * to one (via oputs() or oprintf()) is saved in memory.
 */

void initialize_output_streams ( void )
{
     int i;
     
     toolate = 1;
     
     for ( i = 0; i < output_stream_count; ++i )
	  streams[i].buffer = NULL;

}

/* open_output_streams()
 *
 * open files associated with each output stream.  dump anything buffered
 * in memory to the file.
 */

void open_output_streams ( void )
{
     int i;
     char *fn;
     char *basename;

     basename = get_parameter ( "output.basename" );
     
     fn = (char *)malloc ( strlen(basename)+50 );

     for ( i = 0; i < output_stream_count; ++i )
     {
          strcpy ( fn, basename );
          strcat ( fn, streams[i].ext );
          streams[i].f = fopen ( fn, streams[i].mode );
          if ( streams[i].f == NULL )
          {
	       error ( E_ERROR, "can't open output file \"%s\".", fn );
          }
          else
          {
               streams[i].valid = 1;
	       if ( streams[i].buffer )
		    fputs ( streams[i].buffer, streams[i].f );
	       FREE ( streams[i].buffer );
	       streams[i].buffer = NULL;
          }
     }

     global_basename = (char *)malloc ( strlen(basename)+1 );
     strcpy ( global_basename, basename );

     set_detail_level ( atoi ( get_parameter ( "output.detail" ) ) );
}

/* oputs()
 *
 * prints a string to an output stream.
 */

void oputs ( int streamid, int detail, char *string )
{
     int i, j;

     if ( streamid == OUT_SYS && !quietmode )
     {
          printf ( string );
          fflush ( stdout );
     }

     if ( detail_level < detail )
          return;
     
     for ( i = 0; i < output_stream_count; ++i )
     {
          if ( streamid == streams[i].id )
          {
               if ( streams[i].valid )
               {
                    fputs ( string, streams[i].f );
                    if ( streams[i].autoflush )
                         fflush ( streams[i].f );
                    break;
               }
	       else
	       {
		    if ( streams[i].buffer )
		    {
			 j = strlen ( streams[i].buffer ) +
			      strlen ( string ) + 1;
			 streams[i].buffer = (char *)REALLOC ( streams[i].buffer, j );
			 strcat ( streams[i].buffer, string );
		    }
		    else
		    {
			 streams[i].buffer = (char *)MALLOC ( strlen ( string ) + 1 );
			 strcpy ( streams[i].buffer, string );
		    }
	       }
          }
     }
}

/* oprintf()
 *
 * prints a formatted string to an output stream.
 */

void oprintf ( int streamid, int detail, char *format, ... )
{
     va_list ap;

     if ( detail_level < detail )
          return;
     
     va_start ( ap, format );
     vsprintf ( buffer, format, ap );
     va_end ( ap );

     oputs ( streamid, detail, buffer );
}

/* output_filehandle()
 *
 * returns the filehandle associated with a given stream.
 */

FILE *output_filehandle ( int streamid )
{
     int i;
     for ( i = 0; i < output_stream_count; ++i )
          if ( streamid == streams[i].id )
               if ( streams[i].valid )
                    return ( streams[i].f );
     return NULL;
}

/* output_stream_close()
 *
 * closes an output stream, if it has been marked as resettable.
 */

void output_stream_close ( int streamid )
{
     int i;
     for ( i = 0; i < output_stream_count; ++i )
          if ( streamid == streams[i].id )
               if ( streams[i].reset )
                    if ( streams[i].valid )
                    {
                         fclose ( streams[i].f );
                         streams[i].valid = 0;
                         break;
                    }
}

/* output_stream_open()
 *
 * reopens a resettable output stream.
 */

void output_stream_open ( int streamid )
{
     char *fn;
     int i;
     
     for ( i = 0; i < output_stream_count; ++i )
          if ( streamid == streams[i].id )
               if ( streams[i].reset )
                    if ( !streams[i].valid )
                    {
                         fn = (char *)malloc ( strlen(global_basename)+50 );
                         
                         strcpy ( fn, global_basename );
                         strcat ( fn, streams[i].ext );
                         streams[i].f = fopen ( fn, streams[i].mode );
                         if ( streams[i].f == NULL )
                         {
                              /* an error. */
                         }
                         else
                              streams[i].valid = 1;
                    }
}

/* output_stream_flush()
 *
 * flushes all the resettable output streams.
 */

void output_stream_flush ( int streamid )
{
     int i;
     
     for ( i = 0; i < output_stream_count; ++i )
          if ( streamid == streams[i].id )
               if ( streams[i].reset )
                    if ( streams[i].valid )
                         fflush ( streams[i].f );
}

/* close_output_streams()
 *
 * closes all the output streams.
 */

void close_output_streams ( void )
{
     int i;

     for ( i = 0; i < output_stream_count; ++i )
     {
          if ( streams[i].valid )
          {
               fclose ( streams[i].f );
               streams[i].valid = 0;
          }
     }

     for ( i = SYSOUTPUTSTREAMS; i < MAXOUTPUTSTREAMS; ++i )
     {
          free ( streams[i].ext );
          free ( streams[i].mode );
     }
     
     free ( global_basename );
}

/* flush_output_streams()
 *
 * flushes all the output streams.
 */

void flush_output_streams ( void )
{
     int i;
     for ( i = 0; i < output_stream_count; ++i )
          if ( streams[i].valid )
               fflush ( streams[i].f );
}

/* error()
 *
 * prints an error message to the system output stream, with a given
 * level of severity.  If the error is fatal, then the error is also
 * printed to stderr and the program exits with status 1.
 */

void error ( int severity, char *format, ... )
{
     va_list ap;

     va_start ( ap, format );
     vsprintf ( buffer, format, ap );
     va_end ( ap );
     strcat ( buffer, "\n" );

     oputs ( OUT_SYS, 0, error_type[severity] );
     oputs ( OUT_SYS, 0, buffer );

     if ( severity == E_FATAL_ERROR )
     {
          fprintf ( stderr, "exiting due to fatal error.\n" );
          exit(1);
     }

}

/* set_detail_level()
 *
 * sets the detail level to a given value.
 */

void set_detail_level ( int newlevel )
{
     if ( newlevel >= MINDETAILLEVEL &&
         newlevel <= MAXDETAILLEVEL )
          detail_level = newlevel;
}

/* test_detail_level()
 *
 * returns 1 if the argument is less than or equal to the current
 * detail level.
 */

int test_detail_level ( int detail )
{
     if ( detail_level < detail )
          return 0;
     return 1;
}
