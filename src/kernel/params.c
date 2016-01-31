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
/** array of parameters, along with the size of the array and the number
  of parameters used. **/
parameter *param;
int param_size = 0, param_alloc = 0;

/* read_parameter_file()
 *
 * reads and parses a parameter file.
 */

void read_parameter_file ( char *filename )
{
     FILE *f;
     char buffer[MAXPARAMLINELENGTH+1];
     char *buf2;
     int buf2size;
     int line = 0;
     int errors = 0;
     int including = 1;
     int flag;
     int buflen, buf2len;
     int cont;

     /** open it. **/
     f = fopen ( filename, "r" );
     if ( f == NULL )
          error ( E_FATAL_ERROR, "can't open parameter file \"%s\".",
                 filename );

     /** lines are read into buffer.  they are appended to buf2 until a
       line which is not continued is hit, then buf2 is parsed and the
       parameter added. **/

     /* initial allocation of buf2. */
     buf2 = (char *)MALLOC ( MAXPARAMLINELENGTH * sizeof ( char ) );
     buf2size = MAXPARAMLINELENGTH;
     buf2[0] = 0;
     buf2len = 0;
     cont = 0;

     while ( fgets ( buffer, MAXPARAMLINELENGTH, f ) != NULL )
     {
          ++line;

          /* remove comments, skip line if it's blank after comment
	     removal. */
          if ( delete_comment ( buffer ) )
               continue;

	  /** if this is a %ifdef, %ifndef, or %endif directive, then
	    set the including variable accordingly. **/
          if ( buffer[0] == '%' )
          {
               flag = 0;
               if ( strncmp ( buffer+1, "ifdef", 5 ) == 0 )
                    including = test_directive ( buffer+6 );
               else if ( strncmp ( buffer+1, "ifndef", 6 ) == 0 )
                    including = !test_directive ( buffer+7 );
               else if ( strncmp ( buffer+1, "endif", 5 ) == 0 )
                    including = 1;
               else
                    flag = 1;

               if ( !flag )
                    continue;
          }

	  /* skip this line unless including is set. */
          if ( !including )
               continue;

	  /** if this is a %define or %undefine directive, then process
	    it. **/
          if ( buffer[0] == '%' )
          {
               if ( strncmp ( buffer+1, "define", 6 ) == 0 )
                    define_directive ( buffer+7 );
               else if ( strncmp ( buffer+1, "undefine", 8 ) == 0 )
                    undefine_directive ( buffer+9 );
               else
                    error ( E_ERROR, "%s line %d: unknown directive.",
                           filename, line );
               continue;
          }

	  /* is this line continued?  (is the last nonwhitespace character 
	     a backslash?) */
          cont = check_continuation ( buffer );

	  /** make buf2 bigger if necessary. **/
          buflen = strlen ( buffer );
          while ( buf2len + buflen > buf2size + 10 )
          {
               buf2size += MAXPARAMLINELENGTH;
               buf2 = (char *)REALLOC ( buf2, buf2size * sizeof ( char ) );
          }

	  /** tack the current line onto buf2. **/
          strcat ( buf2, buffer );
          buf2len += buflen;

	  /** if this line is not continued further, then parse buf2. and
	   add the parameter. **/
          if ( !cont )
          {
               if ( parse_one_parameter ( buf2 ) )
               {
                    errors = 1;
                    error ( E_ERROR, "%s line %d: syntax error.",
                           filename, line );
               }

	       /** reset buf2 as empty. **/
               buf2len = 0;
               buf2[0] = 0;
          }
     }

     /* close file. */
     fclose ( f );

     /** if the last line was continued, and nothing followed it, that's
       an error. **/
     if ( buf2len != 0 )
     {
          errors = 1;
          error ( E_ERROR, "unexpected EOF." );
     }

     /** if any errors occurred during parsing, stop now. **/
     if ( errors )
          error ( E_FATAL_ERROR,
                 "some syntax errors occurred parsing \"%s\".", filename );

     FREE ( buf2 );

}

/* check_continuation()
 *
 * returns true/false depending on whether the given string is continued
 * (the last nonwhitespace character is a backslash).  If so, the backslash
 * and everything after it is chopped.
 */

int check_continuation ( char *buffer )
{
     int i, l;
     l = strlen ( buffer );
     for ( i = l-1; i >= 0; --i )
          if ( !isspace(buffer[i]) )
          {
               if ( buffer[i] == '\\' )
               {
                    buffer[i] = '\n';
                    buffer[i+1] = 0;
                    return 1;
               }
               else
                    return 0;
          }
     
     /* a blank line was passed. */
     return 0;
}

/* delete_comment()
 *
 * this searches the string for a '#' or ';' and chops everything
 * following, if one is found.  returns 1 if the resulting line is blank
 * (all whitespace), 0 otherwise.
 */

int delete_comment ( char *buffer )
{
     int i, l;

     l = strlen ( buffer );
     /* zero-length lines are considered blank. */
     if ( l == 0 )
          return 1;
     /* if the last character is a newline, chop it. */
     if ( buffer[--l] == '\n' )
          buffer[l] = 0;
     for ( i = 0; i < l; ++i )
          if ( buffer[i] == '#' || buffer[i] == ';' )
          {
	       /* chop the line at a '#' or ';'. */
               buffer[i] = 0;
               break;
          }
     /* look for a nonwhitespace character. */
     l = strlen ( buffer );
     for ( i = 0; i < l; ++i )
          if ( !isspace(buffer[i]) )
	       /* found one, return 0. */
               return 0;

     /* blank line, return 1. */
     return 1;
}

/* translate_binary()
 *
 * this translates all of the valid strings representing binary values
 * to the corresponding integer.
 */

int translate_binary ( char *string )
{
     if ( strcmp ( string, "true" ) == 0 ||
         strcmp ( string, "t" ) == 0 ||
         strcmp ( string, "on" ) == 0 ||
         strcmp ( string, "yes" ) == 0 ||
         strcmp ( string, "y" ) == 0 ||
         strcmp ( string, "1" ) == 0 )
          return 1;
     else if ( strcmp ( string, "false" ) == 0 ||
              strcmp ( string, "f" ) == 0 ||
              strcmp ( string, "off" ) == 0 ||
              strcmp ( string, "no" ) == 0 ||
              strcmp ( string, "n" ) == 0 ||
              strcmp ( string, "0" ) == 0 )
          return 0;
     else
          return -1;
}

/* read_parameter_database()
 *
 * this reads parameters from a checkpoint file.
 */

void read_parameter_database ( FILE *f )
{
     int i, j, k;
     int count;
     char *name, *value;
     char *buf1, *buf2;
     int buf2len, buf2alloc;

     /* how many parameters are we supposed to find? */
     fscanf ( f, "%*s %d\n", &count );
     if ( fgetc ( f ) != '#' )
	  error ( E_FATAL_ERROR, "error in parameter section of checkpoint file." );

     buf1 = (char *)MALLOC ( MAXPARAMLINELENGTH );
     buf2 = (char *)MALLOC ( MAXPARAMLINELENGTH );
     buf2alloc = MAXPARAMLINELENGTH;
     buf2[0] = 0;
     buf2len = 0;
     
     for ( i = 0; i < count; )
     {
	  /* get a line in buf1. */
	  fgets ( buf1, MAXPARAMLINELENGTH, f );

	  /** lengthen buf2 if necessary. */
	  while ( buf2len + strlen ( buf1 ) >= buf2alloc )
	  {
	       buf2 = (char *)REALLOC ( buf2, buf2alloc + MAXPARAMLINELENGTH );
	       buf2alloc += MAXPARAMLINELENGTH;
	  }
	  /** tack line onto buf2. **/
	  strcat ( buf2, buf1 );
	  buf2len += strlen ( buf1 );

	  /* get the first character of the next line. */
	  k = fgetc ( f );
	  if ( k != '+' )
	  {
	       /** not a '+', so the line is not continued. **/

	       /* chop the final newline. */
	       buf2[buf2len-1] = 0;
	       for ( j = 0; j < buf2len; ++j )
		    /* look for a " = " substring, and break it into name/value
		       there. */
		    if ( buf2[j] == ' ' && buf2[j+1] == '=' &&
			buf2[j+2] == ' ' )
		    {
			 /** add the parameter. **/
			 buf2[j] = 0;
			 add_parameter ( buf2, buf2+j+3,
					PARAM_COPY_NAME|PARAM_COPY_VALUE );
#ifdef DEBUG
			 fprintf ( stderr, "name = [%s]\nvalue = [%s]\n",
				  buf2, buf2+j+3 );
#endif
			 break; /* bug fix, rsteele@ist.flinders.edu.au */
		    }

	       /** reset buf2. **/
	       buf2[0] = 0;
	       buf2len = 0;

	       /* count of how many we've found. */
	       ++i;
	  }
     }
     /* put the extra character we read back. */
     ungetc ( k, f );

     FREE ( buf1 );
     FREE ( buf2 );
}

/* write_parameter_database()
 *
 * this writes all the parameters to a checkpoint file, as "name = value\n".
 * since parameters can have embedded newlines, we begin each line of the
 * file with a "#" to indicate the start of a new name/value pair or a "+"
 * to indicate a continuation of the previous line.
 */

void write_parameter_database ( FILE *f )
{
     int i, j;

     /* write the total count. */
     fprintf ( f, "parameter-count: %d\n", param_size );
     for ( i = 0; i < param_size; ++i )
     {
	  /* start the pair with a '#'. */
	  fputc ( '#', f );
	  /** write the name, adding '+' after newlines. */
	  for ( j = 0; j < strlen ( param[i].n ); ++j )
	  {
	       fputc ( param[i].n[j], f );
	       if ( param[i].n[j] == '\n' )
		    fputc ( '+', f );
	  }
	  /* write " = ". */
	  fputs ( " = ", f );
	  /** write the value, adding '+' after newlines. */
	  for ( j = 0; j < strlen ( param[i].v ); ++j )
	  {
	       fputc ( param[i].v[j], f );
	       if ( param[i].v[j] == '\n' )
		    fputc ( '+', f );
	  }
	  /* end the pair. */
	  fputc ( '\n', f );
     }
}

/* initialize_parameters()
 *
 * initializes the parameter database. */

void initialize_parameters ( void )
{
     oputs ( OUT_SYS, 30, "    parameter database.\n" );
     
     param = (parameter *)MALLOC ( PARAMETER_MINSIZE * sizeof ( parameter ) );
     param_alloc = PARAMETER_MINSIZE;
     param_size = 0;
}

/* free_parameters()
 *
 * frees all the parameters.
 */

void free_parameters ( void )
{
     int i;

     for ( i = 0; i < param_size; ++i )
     {
	  /* if add_parameter made a copy of the name, then free it. */
          if ( param[i].copyflags & PARAM_COPY_NAME )
               FREE ( param[i].n );
	  /* if add_parameter make a copy of the value, then free it. */
          if ( param[i].copyflags & PARAM_COPY_VALUE )
               FREE ( param[i].v );
     }
     
     FREE ( param );
     param = NULL;
     param_alloc = 0;
     param_size = 0;
}

/* add_parameter()
 *
 * adds the given name/value pair to the database.  the flags indicate
 * which if any of the strings need to be copied.
 */

void add_parameter ( char *name, char *value, int copyflags )
{

     /* erase any existing parameter of the same name. */
     delete_parameter ( name );

     /** if the database is full, make it bigger. **/
     while ( param_alloc < param_size+1 )
     {
          param_alloc += PARAMETER_CHUNKSIZE;
          param = (parameter *)REALLOC ( param,
                                       param_alloc * sizeof ( parameter ) );
     }

     /** add the name. **/
     if ( copyflags & PARAM_COPY_NAME )
     {
	  /* make a copy of the string if requested. */
          param[param_size].n = (char *)MALLOC ( strlen(name)+1 );
          strcpy ( param[param_size].n, name );
     }
     else
	  /* just store the pointer passed to us. */
          param[param_size].n = name;

     /** add the value. **/
     if ( copyflags & PARAM_COPY_VALUE )
     {
	  /* make a copy of the string if requested. */
          param[param_size].v = (char *)MALLOC ( strlen(value)+1 );
          strcpy ( param[param_size].v, value );
     }
     else
	  /* just store the pointer passed to us. */
          param[param_size].v = value;

     /* record whether our values are copies or not. */
     param[param_size].copyflags = copyflags;
     
     ++param_size;

}

/* delete_parameter()
 *
 * deletes a parameter from the database.
 */

int delete_parameter ( char *name )
{
     int i;

     for ( i = 0; i < param_size; ++i )
          if ( strcmp ( name, param[i].n ) == 0 )
          {
	       /** free any copies make by add_parameter. **/
               if ( param[i].copyflags & PARAM_COPY_NAME )
                    FREE ( param[i].n );
               if ( param[i].copyflags & PARAM_COPY_VALUE )
                    FREE ( param[i].v );

	       /** move the last value in the database to the position
		 of the deleted one. **/
               if ( param_size-1 != i )
               {
                    param[i].n = param[param_size-1].n;
                    param[i].v = param[param_size-1].v;
                    param[i].copyflags = param[param_size-1].copyflags;
               }
               --param_size;
               return 1;
          }

     return 0;
}

/* get_parameter()
 *
 * looks up a parameter in the database.
 */

char *get_parameter ( char *name )
{
     int i;

     for ( i = 0; i < param_size; ++i )
          if ( strcmp ( name, param[i].n ) == 0 )
               return param[i].v;
     return NULL;
}

/* print_parameters()
 *
 * dumps parameter database to stdout.
 */

void print_parameters ( void )
{
     int i;

     for ( i = 0; i < param_size; ++i )
          printf ( "name: \"%s\"  value: \"%s\"  copy: %d\n",
                  param[i].n, param[i].v, param[i].copyflags );
     
}

/* parse_one_parameter()
 *
 * breaks a string at the first equals sign into name and value parts.
 * removes leading and trailing whitespace from both parts. inserts the
 * resulting pair into the parameter database.
 */

int parse_one_parameter ( char *buffer )
{
     char name[MAXPARAMLINELENGTH+1];
     char data[MAXPARAMLINELENGTH+1];
     int i, j, k, l;
     int n, d;

     k = -1;
     j = 0;
     l = strlen ( buffer );
     /** scan for a equals sign. **/ 
     for ( i = 0; i < l; ++i )
     {
	  /* j records whether or not we have found a nonwhitespace
	    character. */
          j += (buffer[i] != ' ' && buffer[i] != '\t' && buffer[i] != '\n');
          if ( buffer[i] == '=' )
          {
               k = i;
	       /* copy the name part. */
               strncpy ( name, buffer, k );
               name[k] = 0;
	       /* copy the value part. */
               strcpy ( data, buffer+k+1 );
               break;
          }
     }

     /* if we found no '=', return an error unless the line was
	completely blank. */
     if ( k == -1 )
          return !!j;

     /* trim leading and trailing whitespace. */
     n = trim_string ( name );
     d = trim_string ( data );

     /** if either section is blank, return an error, otherwise add
       the pair as a parameter. **/
     if ( n == 0 || d == 0 )
          return 1;
     else
          add_parameter ( name, data, PARAM_COPY_NAME|PARAM_COPY_VALUE );

     return 0;
}

/* trim_string()
 *
 * trims leading and trailing whitespace from a string, overwriting the
 * argument with the result.  returns number of characters in result.
 */

int trim_string ( char *string )
{
     int i, j, l;

     j = -1;
     l = strlen ( string );
     for ( i = 0; i < l; ++i )
     {
          if ( j == -1 )
          {
               if ( string[i] != ' ' && string[i] != '\t' &&
                    string[i] != '\n' )
               {
                    j = i;
                    --i;
               }
          }
          else
               string[i-j] = string[i];
     }
     if ( j == -1 )
     {
          string[0] = 0;
          return 0;
     }
     
     string[i-j] = 0;
     l = i-j;
     j = -1;
     for ( i = 0; i < l; ++i )
     {
          if ( string[i] != ' ' && string[i] != '\t' && string[i] != '\n' )
               j = i;
     }
     string[j+1] = 0;

     return j+1;
}

/* define_directive()
 *
 * defines a directive "SYMBOL", which is just a parameter called
 * "__define:SYMBOL".  trims leading and trailing whitespace from SYMBOL.
 */

void define_directive ( char *string )
{
     char *buffer;
     int i;

     for ( i = 0; i < strlen(string) && isspace(string[i]); ++i );
     
     buffer = (char *)MALLOC ( (20 + strlen(string)) * sizeof ( char ) );
     strcpy ( buffer, "__define:" );
     strcat ( buffer, string+i );

     for ( i = strlen(buffer)-1; i >= 0 && isspace(buffer[i]); --i )
          buffer[i] = 0;
     
     add_parameter ( buffer, "1", PARAM_COPY_NAME );
     FREE ( buffer );
}

/* undefine_directive()
 *
 * undefines a directive "SYMBOL".
 */

void undefine_directive ( char *string )
{
     char *buffer;
     int i;

     for ( i = 0; i < strlen(string) && isspace(string[i]); ++i );

     buffer = (char *)MALLOC ( (20 + strlen(string)) * sizeof ( char ) );
     strcpy ( buffer, "__define:" );
     strcat ( buffer, string+i );

     for ( i = strlen(buffer)-1; i >= 0 && isspace(buffer[i]); --i )
          buffer[i] = 0;
     
     delete_parameter ( buffer );
     FREE ( buffer );
}

/* test_directive()
 *
 * returns 1 iff a given directive is defined.
 */

int test_directive ( char *string )
{
     char *buffer;
     int ret;
     int i;

     for ( i = 0; i < strlen(string) && isspace(string[i]); ++i );

     buffer = (char *)MALLOC ( (20 + strlen(string)) * sizeof ( char ) );
     strcpy ( buffer, "__define:" );
     strcat ( buffer, string+i );

     for ( i = strlen(buffer)-1; i >= 0 && isspace(buffer[i]); --i )
          buffer[i] = 0;
     
     if ( get_parameter ( buffer )  )
          ret = 1;
     else
          ret = 0;

     FREE ( buffer );
     return ret;
}

/* binary_parameter()
 *
 * checks for the existence of a parameter.  if it exists, then it is changed
 * to the string "0" or "1" using lilgp's list of strings representing
 * binary values.  if the value is not on the list, or the parameter is
 * not found, then the parameter is set according to the value argument (it
 * acts as a default).
 */

void binary_parameter ( char *name, int value )
{
     char *param = get_parameter ( name );
     char string[2];
     char *i, *is;
     int v;

     if ( param != NULL )
     {
	  /* copy the value and lowercase it. */
          v = strlen ( param );
          i = (char *)MALLOC ( (v+1)*sizeof ( char ) );
          strcpy ( i, param );
          for ( is = i; *is; ++is )
               *is = tolower(*is);

	  /* translate to a binary integer. */
          v = translate_binary ( i );
          
          if ( v == -1 )
          {
	       /* translation failed, use the value argument. */
               error ( E_ERROR,
                      "\"%s\" is not a legal value for \"%s\"; assuming default.",
                      i, name );
               v = value;
          }

          FREE ( i );
               
     }
     else
	  /* parameter not found, use the value argument. */
          v = value;

     /** print the value to a string and put it in the parameter database. */
     sprintf ( string, "%d", !!v );
     add_parameter ( name, string, PARAM_COPY_VALUE|PARAM_COPY_NAME );

}


               

