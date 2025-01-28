/*
Program:   fastaselecth.c
Version:   1.0.11
Date:      20-MAY-2019
Author:    David Mathog, Biology Division, Caltech
email:     mathog@caltech.edu
Copyright: 2019 David Mathog and California Institute of Technology (Caltech)

Description:

  Similar to FASTARANGE except it reads a list of fasta entry
  names from an input file and then in a single pass through the
  fasta file emits those entries which match in the order
  specified in the list file.  The name is the
  part following > and up to the first space, tab or |.

  In order to conserve memory the program emits matching
  strings whenever the next one in order is available.
  In the worst possible case though it would have to
  hold all strings in memory until the very last one
  was found, and that could cause an out of memory condition.

  See also fastaselecti which performs a similar function but does
  it by matching entry number 1->N in the file.


Changes:

  1.0.12 20-MAY-2019
         Added -cod, duplicates not fatal in select list, still fatal in
         fasta file.
  1.0.11 07-JAN-2019
         Added check for missing group with -frag[ac].
         Improved explanation of -frag[ac] -in -h.
  1.0.9  06-NOV-2018
         Added -reject mode.  
  1.0.8  15-OCT-2018
         Added -frag mode.  
  1.0.7  10-OCT-2018
         -ht only applies to -sel input, added -hi to apply to -input headers.
  1.0.6  24-SEP-2018
         Replaced final search loop with a call to a binary search.
         Previous method (Nkeys x Mseqs) worked well when one of the numbers was small
         even if the other was very large.  This works well even when both are large. Using
         a hash instead might be slightly faster, but this is fast enough.
  1.0.5  25-APR-2018
         Fixed -dl/-ht confusion.
  1.0.4  18-JAN-2017
         Added command line switches.  Old interface is deprecated.
  1.0.3  02-NOV-2010, increased buffer size to 10M, hopefully keep ahead of exploding
         NCBI header sizes.
  1.0.2  13-JUL-2009, added an optional 2nd argument which contains the list of characters which
         terminate a header.  "| :", \t,\n and \0 always terminate, whether they are added to list or not.
  1.0.1  15-MAY-2007, modified string comparison slightly so that an entry like "foo" in list
         will match >foo followed by space, tab, null (EOL) or \1 (control A).
  1.0.0  19-MAY-2003, David Mathog, initial release.
  
License terms:
    You may run this program on any platform. You may
    redistribute the source code of this program subject to
    the condition that you do not first modify it in any way.
    You may  distribute binary versions of this program so long
    as they were compiled from unmodified source code.  There
    is no charge for using this software.  You may not charge
    others for the use of this software.
License:  GNU General Public License 2
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

Miscellaneous:
    This should be portable.  Compile like:
    
    gcc -O3 -Wall -std=c99 -pedantic -o fastaselecth fastaselecth.c
    

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <ctype.h>

/* definitions and enums */
#define EXVERSTRING "1.0.11  20-MAY-2019"
#define COPYSTRING  "2019 David Mathog and California Institute of Technology"
#define BUGSTRING   "mathog@caltech.edu"
#define LICSTRING   "GNU General Public License 2"

#define MYMAXSTRING 10000000
#define DEFENTRIES     32000

#define FRAG_NONE   0
#define FRAG_NEW    1
#define FRAG_APPEND 2

/*function prototypes */
int  bin_search(char *find, char **list, int size );
int  convert_escape(char *string);
void emit_help(void);
void emit_hhead(void);
int  get_entries(char *bigstring, char ***header_name_list, char ***group_name_list);
void insane(char *string);
int  lcl_strcasecmp(const char *s1, const char *s2);
char *lcl_strdup(const char *string);
void remove_dups(char **header_name_list, char **group_name_list, int *emit_order, int *entrynum);
void setirangenumeric(int *val,int *numarg, int lower, int upper, int argc,char **argv,char * label);
void sort_entries(char **header_name_list, char **group_name_list, int *order, int entrynum);
void process_command_line_args(int argc,char **argv);

/* global variables */
char *gbl_hs;
char *gbl_hi;
char *gbl_in;
char *gbl_sel;
char *gbl_out;
int   gbl_frag;
int   gbl_com;
int   gbl_cod;
int   gbl_wl;
int   gbl_reject;


/* functions */

void insane(char *string){
 (void) fprintf(stderr,"%s\n",string);
 exit(EXIT_FAILURE);
}

void setirangenumeric(int *val,int *numarg, int lower, int upper, int argc,char **argv,char * label){
  (*numarg)++;
  if( ( *numarg >= argc ) || (argv[*numarg] == NULL)){
    (void) fprintf( stderr, "%s: missing argument\n",label);
    exit(EXIT_FAILURE);
  }
  if(sscanf(argv[*numarg],"%d",val) != 1){
    (void) fprintf(stderr,"Bad integer argument/parameter [%s %s] \n",label,argv[*numarg]);
    exit(EXIT_FAILURE);
  }
  if(*val < lower){
    (void) fprintf(stderr,"Illegal (too low) integer argument/parameter [%s %s] \n",label,argv[*numarg]);
    exit(EXIT_FAILURE);
  }
  if(*val > upper){
    (void) fprintf(stderr,"Illegal (too high) integer argument/parameter [%s %s] \n",label,argv[*numarg]);
    exit(EXIT_FAILURE);
  }
}

/*  Read all of the entries to match from gbl_sel.

   bigstring          a buffer
   header_name_list   pointer to an array of character pointers
   
   Returns the number of names to search for.

*/
int get_entries(char *bigstring, char ***header_name_list, char ***group_name_list){
   int end,size;
   char **newlist;
   char *newline;
   char *newstring;
   FILE *fin;
   int spanned;

   if(strcmp(gbl_sel,"-")){
      fin = fopen(gbl_sel,"r");
      if(fin==NULL)insane("fastaselecth: fatal error: could not read input file");
   }
   else {
      fin = stdin;
   }
   
   size=DEFENTRIES;
   end=0;

   // initial allocation 
   newlist=malloc(size*sizeof(char *));
   if(newlist==NULL)insane("fastaselecth: fatal error: could not allocate memory");
   *header_name_list=newlist;

   // initial allocation
   if(gbl_frag){
      newlist=malloc(size*sizeof(char *));
      if(newlist==NULL)insane("fastaselecth: fatal error: could not allocate memory");
      *group_name_list=newlist;
   }
   else {
      *group_name_list=NULL;
   }

   while (fgets(bigstring,MYMAXSTRING,fin) != NULL){
     newline=strstr(bigstring,"\n");
     if(newline != NULL){  
       *newline='\0';  /* replace the \n with a terminator */
       newline--;
     }
     else{ /* string truncated, record too long or EOF */
       if(!feof(fin)){
         (void) fprintf(stderr,"fastaselecth input record in fasta file exceeds %d characters\n",MYMAXSTRING); 
         exit(EXIT_FAILURE);
       }
       (void) fprintf(stderr,"fastaselecth warning: last line of file lacks a \\n \n"); 
       newline=&(bigstring[strlen(bigstring) - 1]);
     }
     if(newline>=bigstring && *newline=='\r')*newline='\0';
      
     
     spanned = strcspn(bigstring,gbl_hs);
     if(spanned>0){  /* ignore empty strings */
       bigstring[spanned] = '\0';
       newstring=malloc((spanned+1)*sizeof(char));
       if(newstring==NULL)insane("fastaselecth: fatal error: could not allocate memory");
       (*header_name_list)[end]=newstring;
       strcpy(newstring,bigstring);
       if(gbl_frag){
          char *rest = bigstring+spanned+1;
          spanned = strspn(rest,gbl_hs);   // consume all delimiters
          rest=rest+spanned;
          spanned = strcspn(rest,gbl_hs);  // find delimiter far side of "rest"
          
          if(spanned>0){  /* ignore empty strings */
            rest[spanned] = '\0';
            newstring=malloc((spanned+1)*sizeof(char));
            if(newstring==NULL)insane("fastaselecth: fatal error: could not allocate memory");
            (*group_name_list)[end]=newstring;
            strcpy(newstring,rest);
          }
       }
       end++;
       if(end >=size){
         size=size+DEFENTRIES;

         newlist=realloc(*header_name_list,size*sizeof(char *));
         if(newlist==NULL)insane("fastaselecth: fatal error: could not reallocate memory");
         *header_name_list=newlist;

         if(gbl_frag){
           newlist=realloc(*group_name_list,size*sizeof(char *));
           if(newlist==NULL)insane("fastaselecth: fatal error: could not reallocate memory");
           *group_name_list=newlist;
         }
       }
     }
   }
   if(fin!=stdin){
      fclose(fin);
   }
   return end;
}


/* modified combsort with restart capability - keeps sort from
degenerating into bubble sort on toxic data */
void sort_entries(char **header_name_list, char **group_name_list, int *order, int entrynum){
int done,gap,swaps;
int restart=0;
int i,j;
char *temp;
float shrink;
int itemp;

  shrink=1.3;
  for(done=0; !done; ){
    gap = (entrynum/2) - 1;
    if(gap<1)gap=1;
    for( ;gap>=1;gap/=shrink){
      if(gap<1)gap=1;
      swaps=0;
      for(i=0,j=gap;j<entrynum;i++,j++){
        if(strcmp(header_name_list[j],header_name_list[i])<0){
          swaps++;
          temp=header_name_list[i];
          header_name_list[i]=header_name_list[j];
          header_name_list[j]=temp;
          if(gbl_frag){
            temp=group_name_list[i];
            group_name_list[i]=group_name_list[j];
            group_name_list[j]=temp;
          }
          itemp=order[i];
          order[i]=order[j];
          order[j]=itemp;
        }
      }
      if(gap==1){
        if(swaps==0){
          done=1;
          break;
        }
        else {
          restart++;
          if(restart > 7)break;
        }
      }
    }
  }
}

void remove_dups(char **header_name_list, char **group_name_list, int *emit_order, int *entrynum){
int i;
   if(*entrynum==1)return;
   char *dst=header_name_list[0];
   int   didx=0;
   for(i=1;i<*entrynum;i++){
       if(strcmp(dst,header_name_list[i]) == 0){
          if(gbl_cod){
             fprintf(stderr,"fastaselecth: warning: duplicate entry name \"%s\" in -sel list, alternate header terminators may be needed\n",dst);
          }
          else {
             insane("fastaselecth: fatal error: duplicate entry names in list, alternate header terminators may be needed");
          }
          free(header_name_list[i]);
       }
       else {
          didx++;
          if(didx != i){
             header_name_list[didx] = header_name_list[i];
             if(gbl_frag){
                group_name_list[didx]  = group_name_list[i];
             }
             emit_order[didx]       = emit_order[i];
          }
          dst = header_name_list[didx];
       }
   }
   *entrynum=didx+1;
}

/* return position found or -1 if not found. */

int bin_search(char *find_me, char **list, int size ){
   int bot = 0;
   int mid;
   int top = size - 1;
   int ret = -1;

   while(bot <= top){
      mid = (bot + top)/2;
      if (strcmp(list[mid], find_me) == 0){
         ret = mid;
         break;
      } else if (strcmp(list[mid], find_me) > 0){
         top = mid - 1;
      } else if (strcmp(list[mid], find_me) < 0){
         bot = mid + 1;
      }
   }
   return(ret);
}


/* Convert text form for special characters to an unsigned character.  Handles:
    These C escape characters (ONLY) \\, \a,\b,\f,\t,\r, and \n.
    ASCII control characters like ^J (masks the 2nd character retaining only the lowest 6 bits)
      For a lone "^" use "/^".
    Numerically specified character values as \###,\o###,\x## (3,3, and 2 digits, as shown ONLY)
      Range is 0-255 only.
    Returns 1 on success, 0 on error
    */
int convert_escape(char *string){
unsigned char *parsed;  /* parsed string, writing is delayed and may overwrite ptr2 */
unsigned char *scan;    /* scanning pointer */
#define NORMAL   0
#define ESCAPE   1
#define CONTROL  2
#define DNUMERIC 3
#define ONUMERIC 4
#define XNUMERIC 5
int state  = NORMAL;
int sum    = 0;
int count  = 0;
int status = 1;
int ok     = 1;
  for(scan = parsed = (unsigned char *) string; ok; scan++){
     switch(state){
       case NORMAL:
         switch(*scan){
           case '\\':
             state=ESCAPE;
             break;
           case '^':
             state=CONTROL;
             break;
           case '\0':
             *parsed = *scan; ok = 0;
             break;
           default:
             state=NORMAL;   *parsed=*scan; parsed++;
         }
         break;
       case ESCAPE:
         switch(*scan){
           case '\\':
             state=NORMAL;   *parsed=*scan; parsed++;
             break;
           case 'a':
             state=NORMAL;   *parsed='\a';  parsed++;
             break;
           case 'b':
             state=NORMAL;   *parsed='\b';  parsed++;
             break;
           case 'f':
             state=NORMAL;   *parsed='\f';  parsed++;
             break;
           case 't':
             state=NORMAL;   *parsed='\t';  parsed++;
             break;
           case 'r':
             state=NORMAL;   *parsed='\r';  parsed++;
             break;
           case 'n':
             state=NORMAL;   *parsed='\n';  parsed++;
             break;
           case 'd':
             state=DNUMERIC; sum=0; count=0;
             break;
           case 'o':
             state=ONUMERIC; sum=0; count=0;
             break;
           case 'x':
             state=XNUMERIC; sum=0; count=0;
             break;
           case '0':
           case '1':
           case '2':
           case '3':
           case '4':
           case '5':
           case '6':
           case '7':
           case '8':
           case '9':
             state=DNUMERIC; sum = *scan - '0'; count=1;
             break;
           case '\0':
             ok = status = 0;
             break;
           default:
             state=NORMAL;   *parsed=*scan; parsed++;
         }
         break;
       case CONTROL:
         if(*scan=='\0'){
           ok = status = 0;
         }
         else {
           state=NORMAL;   *parsed = *scan & 31; parsed++;
         }
         break;
       case DNUMERIC:
         switch(*scan){
           case '0':
           case '1':
           case '2':
           case '3':
           case '4':
           case '5':
           case '6':
           case '7':
           case '8':
           case '9':
             sum=(10*sum) + (*scan - '0');
             if(++count == 3){  /* There must be exactly 3 digits */
               if(sum > UCHAR_MAX){ 
                 ok = status = 0;
               }
               else {
                 state = NORMAL; *parsed=sum; parsed++;
               }
             }
             break;
           default:
             ok = status = 0;
         }
         break;
       case ONUMERIC:
         switch(*scan){
           case '0':
           case '1':
           case '2':
           case '3':
           case '4':
           case '5':
           case '6':
           case '7':
             sum= (8*sum) + (*scan - '0');
             if(++count == 3){ /* There must be exactly 3 digits */
               if(sum > UCHAR_MAX){ 
                 ok = status = 0;
               }
               else {
                 state = NORMAL; *parsed=sum; parsed++;
               }
             }
             break;
           default:
             ok = status = 0;
         }
         break;
       case XNUMERIC:
         switch(*scan){
           case '0':
           case '1':
           case '2':
           case '3':
           case '4':
           case '5':
           case '6':
           case '7':
           case '8':
           case '9':
             sum=(16*sum) + (*scan - '0');
             if(++count == 2){ /* There must be exactly 2 digits */
               state = NORMAL; *parsed=sum; parsed++;
             }
             break;
           case 'A':
           case 'B':
           case 'C':
           case 'D':
           case 'E':
           case 'F':
             sum=(16*sum) + (10 + *scan - 'A');
             if(++count == 2){ /* There must be exactly 2 digits */
               state = NORMAL; *parsed=sum; parsed++;
             }
             break;
           case 'a':
           case 'b':
           case 'c':
           case 'd':
           case 'e':
           case 'f':
             sum=(16*sum) + (10 + *scan - 'a');
             if(++count == 2){ /* There must be exactly 2 digits */
               state = NORMAL; *parsed=sum; parsed++;
             }
             break;
           default:
             ok = status = 0;
         }
         break;
       default:  /* catch what should be an impossible state */
         ok = status = 0;
     }
   }
   return(status);
}

void emit_help(void){
   (void) fprintf(stderr,"Usage: fastaselecth [options]\n");
   (void) fprintf(stderr,"       select a subset of records in a fastafile by header values.\n\n");
   (void) fprintf(stderr,"Command line options:\n");
   (void) fprintf(stderr,"   -in FILE\n");
   (void) fprintf(stderr,"         Read fasta records from FILE.\n");
   (void) fprintf(stderr,"   -out FILE\n");
   (void) fprintf(stderr,"         Selected records go to FILE.  If omitted or FILE is \"-\" write to stdout instead..\n");
   (void) fprintf(stderr,"         If -frag[ca] is set FILE must be like \"template_%%s.fasta\"\n");
   (void) fprintf(stderr,"   -sel FILE\n");
   (void) fprintf(stderr,"         Name of a file containing record selection information.  Default or \"-\" is stdin.\n");
   (void) fprintf(stderr,"         If -frag[ca] is set every select string must have two fields: select and group.\n");
   (void) fprintf(stderr,"         The group field fills in the %%s in the output file name.  Multiple selections may be\n");
   (void) fprintf(stderr,"         directed to each output group file.\n");
   (void) fprintf(stderr,"   -com\n");
   (void) fprintf(stderr,"         Continue On Miss.  If a specified selector has no corresponding input record\n");
   (void) fprintf(stderr,"         a fatal error occurs.  If -com is specified a warning is issued\n");
   (void) fprintf(stderr,"         and processing continues.\n");
   (void) fprintf(stderr,"   -cod\n");
   (void) fprintf(stderr,"         Continue On Duplicates.  If a selector is a duplicate of another \n");
   (void) fprintf(stderr,"         normally a fatal error occurs.  If -cod is specified a warning is\n");
   (void) fprintf(stderr,"         issued, only one copy of the selector is retained, and processing continues.\n");
   (void) fprintf(stderr,"         A selector matching more than one fasta header is a fatal error unless.\n");
   (void) fprintf(stderr,"         -reject is also specified.  Such a selector may not trigger an error\n");
   (void) fprintf(stderr,"         without -reject if all other selectors have already matched, as the program\n");
   (void) fprintf(stderr,"         will exit normally at the first match and so never encounter the duplicates.\n");
   (void) fprintf(stderr,"   -fragc\n");
   (void) fprintf(stderr,"         Direct selections to multiple output files, none of which may exist.  Each group must be in\n");
   (void) fprintf(stderr,"         a contiguous series of select entries or a fatal error occurs.\n");
   (void) fprintf(stderr,"   -fraga\n");
   (void) fprintf(stderr,"         Direct selections to multiple output files, which may exist, and entries will be appended\n");
   (void) fprintf(stderr,"         to them.  Groups need not be clustered in the selection input.\n");
   (void) fprintf(stderr,"   -reject\n");
   (void) fprintf(stderr,"         Reject selected entries.  Default is to accept selected entries.  Not with -frag[ac]\n");
   (void) fprintf(stderr,"   -wl N\n");
   (void) fprintf(stderr,"         Width of Longest input line.  Default is %d.\n",MYMAXSTRING);
   (void) fprintf(stderr,"   -hs STRING\n");
   (void) fprintf(stderr,"   -ht STRING\n");
   (void) fprintf(stderr,"         Specify an alternate set of -sel FILE delimiters.  The first delimiter\n");
   (void) fprintf(stderr,"         encountered terminates the string.  Default delimiters are:\n");
   (void) fprintf(stderr,"         EOL, NULL, tab, space, vertical bar, and colon.  Possible values include:\n");
   (void) fprintf(stderr,"            C escape sequences: \\\\, \\a, \\b, \\f, \\t, \\r, and \\n;\n");
   (void) fprintf(stderr,"            Control characters like ^C;\n");
   (void) fprintf(stderr,"            Numeric character values: \\###, \\d###, \\o###, and \\x## (digital,digital,octal, and hex).\n");
   (void) fprintf(stderr,"   -hi STRING\n");
   (void) fprintf(stderr,"         Specify an alternate set of -input FILE delimiters.  Default is \"\\1 \\t\".\n");
   (void) fprintf(stderr,"         Syntax is the same as for -hs.\n");
   (void) fprintf(stderr,"   -h    Print this help message (also -help --h --help -? --?)\n");
   (void) fprintf(stderr,"   -hhead\n");
   (void) fprintf(stderr,"         Print explanation of header selection and header delimiters.\n");
   (void) fprintf(stderr,"   -i    Emit version, copyright, license and contact information\n\n");
}

void emit_hhead(void){
   (void) fprintf(stderr,"Fasta files contain one or more entries.\n");
   (void) fprintf(stderr,"Each entry starts with a header line which begins with \">NAME\".\n");
   (void) fprintf(stderr,"The rest of the entry contains any number of data lines.\n");
   (void) fprintf(stderr,"Data lines may hold any type of text other than another header line.\n");
   (void) fprintf(stderr,"The delimiter set applied is determined by the -hi parameter.\n");
   (void) fprintf(stderr,"\n");
   (void) fprintf(stderr,"Select files contain a series of entry names, one per line, terminated by a delimiter.\n");
   (void) fprintf(stderr,"If -fragc or -fraga is used the entry is followed by a group name which is .\n");
   (void) fprintf(stderr,"  used to construct the output file name.\n");
   (void) fprintf(stderr,"The delimiter set applied is determined by the -ht parameter.\n");
   (void) fprintf(stderr,"Entry names correspond to the \">NAME\" part of a fasta header line.\n");
   (void) fprintf(stderr,"Those fasta entries entries whose NAME matches exactly will be emitted\n");
   (void) fprintf(stderr,"to -out in the order in which they appear in the list file.\n");
   (void) fprintf(stderr,"\n");
   (void) fprintf(stderr,"Example: if the list file contains three lines \"YACL12\",\"SLACL2\", and \"LLEV12\"\n");
   (void) fprintf(stderr,"those entries from the fasta file will be emitted in the order YACL12,SLACL2,LLEV12.\n");
}

int  lcl_strcasecmp(const char *s1, const char *s2){
int c1;
int c2;
  for(; ;s1++,s2++){
    c1=toupper(*s1);
    c2=toupper(*s2);
    if(c1 < c2)return -1;
    if(c1 > c2)return  1;
    if(c1 == 0)return  0;  /*c2 also is 0 in this case */
  }
}

char *lcl_strdup(const char *string){
   int slen = 1 + strlen(string);
   char *retval = malloc(slen);
   if(!retval)insane("fastaselecth: fatal error: could not allocate memory");
   strcpy(retval,string);
   return(retval);
}

void process_command_line_args(int argc,char **argv){
   int numarg=0;
   gbl_hs  = lcl_strdup("|\t :");
   gbl_hi  = lcl_strdup("\1\t ");
   gbl_in  = NULL;
   gbl_sel = NULL;
   gbl_out = NULL;
   gbl_frag= FRAG_NONE;
   gbl_com = 0;
   gbl_cod = 0;
   gbl_wl  = MYMAXSTRING;
   gbl_reject = 0;

   while( ++numarg < argc){
      if( (lcl_strcasecmp(argv[numarg], "-h")==0)     ||
          (lcl_strcasecmp(argv[numarg], "--h")==0)    ||
          (lcl_strcasecmp(argv[numarg], "-?")==0)     ||
          (lcl_strcasecmp(argv[numarg], "--?")==0)    ||
          (lcl_strcasecmp(argv[numarg], "-help")==0)  ||
          (lcl_strcasecmp(argv[numarg], "--help")==0) ){
         emit_help();
         exit(EXIT_SUCCESS);
      }
      else if(lcl_strcasecmp(argv[numarg], "-i")==0){
         (void)fprintf(stderr,"Version:   %s\n",EXVERSTRING);
         (void)fprintf(stderr,"bugs to:   %s\n",BUGSTRING);
         (void)fprintf(stderr,"Copyright: %s\n",COPYSTRING);
         (void)fprintf(stderr,"License:   %s\n",LICSTRING);
         exit(EXIT_SUCCESS);
      }
      else if(lcl_strcasecmp(argv[numarg], "-hhead")==0){
         emit_hhead();
         exit(EXIT_SUCCESS);
      }
      else if(lcl_strcasecmp(argv[numarg], "-in")==0){
         gbl_in = argv[++numarg];
      }
      else if(lcl_strcasecmp(argv[numarg], "-out")==0){
         gbl_out = argv[++numarg];
      }
      else if(lcl_strcasecmp(argv[numarg], "-sel")==0){
         gbl_sel = argv[++numarg];
      }
      else if(lcl_strcasecmp(argv[numarg], "-fragc")==0){
         gbl_frag = FRAG_NEW;
      }
      else if(lcl_strcasecmp(argv[numarg], "-fraga")==0){
         gbl_frag = FRAG_APPEND;
      }
      else if(lcl_strcasecmp(argv[numarg], "-reject")==0){
         gbl_reject = 1;
      }
      else if(lcl_strcasecmp(argv[numarg], "-wl")==0){
         setirangenumeric(&gbl_wl,&numarg,1,INT_MAX,argc,argv,"-wl");
      }
      else if((lcl_strcasecmp(argv[numarg], "-ht")==0) || (lcl_strcasecmp(argv[numarg], "-hs")==0)){
         if(gbl_hs){ free(gbl_hs); }
         gbl_hs = malloc(sizeof(char) * (1 + strlen(argv[++numarg])));
         strcpy(gbl_hs, argv[numarg]);
         if(!convert_escape(gbl_hs)){
           insane("fastaselecth: fatal error: select header terminator string had syntax error");
         }
      }
      else if(lcl_strcasecmp(argv[numarg], "-hi")==0){
         if(gbl_hi){ free(gbl_hi); }
         gbl_hi = malloc(sizeof(char) * (1 + strlen(argv[++numarg])));
         strcpy(gbl_hi, argv[numarg]);
         if(!convert_escape(gbl_hi)){
           insane("fastaselecth: fatal error: file header terminator string had syntax error");
         }
      }
      else if(lcl_strcasecmp(argv[numarg], "-com")==0){
         gbl_com=1;
      }
      else if(lcl_strcasecmp(argv[numarg], "-cod")==0){
         gbl_cod=1;
      }
      else {
         (void) fprintf(stderr,"Unknown command line argument: %s\n",argv[numarg]);
         emit_help();
         exit(EXIT_FAILURE);
      }
   }

   /* sanity checking */
   if(!gbl_in)insane("fastaselecth: fatal error: no -in specified");
   if(!gbl_sel )insane("fastaselecth: fatal error: -sel must be specified");
   if(gbl_frag && !strstr(gbl_out,"%s"))insane("fastaselecth: fatal error: -frag set but -out does not contain %s");
   if(gbl_frag && gbl_reject)insane("fastaselecth: fatal error: -frag cannot be combined with -reject");
}

int main(int argc, char *argv[]){
   char *newline=NULL;
   char *tbuf=NULL;
   char **header_name_list=NULL;
   char **group_name_list=NULL;
   int  *emitlist=NULL;
   int  emitting;
   int  *emitorder=NULL;
   char **emitstrings=NULL;
   char **emitgroups=NULL;
   char *accumstring=NULL;
   int  entrynum,emit;
   int  DONE;
   int  i,tail,size;
   int  lastemitted;
   char *bptr=NULL;
   char *last_group;
   char empty_string[]="";
   char temp_name[1028];
   
   unsigned long long records;
   unsigned long long emitted;

   process_command_line_args(argc,argv);

   char *bigstring=malloc(gbl_wl + 1);
   if(!bigstring)insane("fastaselecth: fatal error: could not allocate memory");
   char *bigheader=malloc(gbl_wl + 1);
   if(!bigheader)insane("fastaselecth: fatal error: could not allocate memory");

   DONE        = 0;
   lastemitted = -1;

   
   entrynum    = get_entries(bigstring, &header_name_list, &group_name_list);
   if(!entrynum)insane("fastaselecth: fatal error: nothing was read from -sel");

   emitlist    = calloc(entrynum,sizeof(int));
   if(emitlist==NULL)insane("fastaselecth: fatal error: could not allocate memory");

   emitorder   =calloc(entrynum,sizeof(int));
   if(emitorder==NULL)insane("fastaselecth: fatal error: could not allocate memory");

   emitstrings =calloc(entrynum,sizeof(char *));
   if(emitstrings==NULL)insane("fastaselecth: fatal error: could not allocate memory");

   if(gbl_frag){
     emitgroups =calloc(entrynum,sizeof(char *));
     if(emitgroups==NULL)insane("fastaselecth: fatal error: could not allocate memory");
   }

   for(i=0;i<entrynum;i++){emitorder[i]=i;}
   sort_entries(header_name_list, group_name_list, emitorder, entrynum);
   remove_dups(header_name_list, group_name_list, emitorder, &entrynum);

   records=0;
   emit=0;
   emitted=0;
   emitting=0;
   accumstring=NULL;
   tail=size=0;
   lastemitted=-1;
   FILE *fin = fopen(gbl_in,"r");
   if(!fin)insane("fastaselecth: fatal error: could not open -in");
   FILE *fout=NULL;
   if(gbl_frag){
      fout = stdout;
      last_group=empty_string;
   }
   else {
      if(!gbl_out || !strcmp(gbl_out,"-")){
         fout = stdout;
      }
      else {
         fout = fopen(gbl_out,"w");
         if(!fout)insane("fastaselecth: fatal error: could not open -out");
      }
   }
   while( fgets(bigstring,MYMAXSTRING,fin) != NULL){
      newline=strstr(bigstring,"\n");
      if(newline != NULL){  
         *newline='\0';  /* replace the \n with a terminator */
         newline--;
      }
      else{ /* string truncated, record too long or EOF */
         if(!feof(stdin)){
            (void) fprintf(stderr,"fastaselecth: fatal error: input record in fasta file exceeds %d characters\n",MYMAXSTRING); 
            exit(EXIT_FAILURE);
         }
        (void) fprintf(stderr,"fastaselecth warning: last line of file lacks a \\n \n"); 
         newline=&(bigstring[strlen(bigstring) - 1]);
      }
      if(newline>=bigstring && *newline=='\r')*newline='\0';
      
      if(bigstring[0] == '>'){
         records++;
      
         /* A new entry.  If the preceding entry was in the emitting state store the pointer to it in emitstrings */
         if(!gbl_reject && accumstring!=NULL){
            if(emitstrings[emitorder[emitting]]!=NULL)insane("fastaselecth: fatal programming error: nonNULL storage");
            emitstrings[emitorder[emitting]]=accumstring;
            accumstring=NULL;
            emitting=0;
            size=0;
            tail=0;
            while(1){
               if(emitstrings[lastemitted+1]!=NULL){  
                  /* next one in order is available to emit.  The idea here is to emit the strings as soon as 
                     possible rather than waiting until all have been collected and then doing them all at once.
                     This is faster since it spreads the writes out over time, which can make a big difference if there
                     are many megabytes of writes. */
                  lastemitted++;
                  if(gbl_frag && strcmp(last_group,emitgroups[lastemitted])){
                     last_group = emitgroups[lastemitted];
                     if(fout){
                        fclose(fout);
                     }
                     sprintf(temp_name,gbl_out,last_group);
                     if(gbl_frag == FRAG_APPEND){
                         fout = fopen(temp_name,"a");
                     }
                     else if (gbl_frag == FRAG_NEW){
                         FILE *fprobe = fopen(temp_name,"r");
                         if(fprobe){
                            fprintf(stderr,"fastaselecth: fatal error: file name: %s\n",temp_name);
                            insane("fastaselecth: fatal error: -fragc mode output file already exists or noncontiguous group records");
                         }
                         else {
                            fout = fopen(temp_name,"w");
                         }
                     }
                     if(!fout){
                        fprintf(stderr,"fastaselecth: fatal error: file name: %s\n",temp_name);
                        insane("fastaselecth: fatal error: could not open output file in -frag mode");
                     }
                  }
                  (void) fprintf(fout,"%s",emitstrings[lastemitted]);
                  free(emitstrings[lastemitted]); /* release memory */
               }
               else {
                  break;
               }
               // There may be more data in the input file but all the selected entries have been found
               if(lastemitted == entrynum - 1){
                  goto bye;
               }
            }
         }

         /*does the name in bigstring match anything in the list?  Here "match" allows a space, tab, ^A or NULL to
           terminate the name in bigstring.  However, the name from the header_name_list must match exactly (end in null). 
           Replace the bigstring terminate with a \0 for the search, then put the original character back.
         */

         bptr=bigstring+1;
         strcpy(bigheader,bptr);
         size_t b_num_chars = 0;
         b_num_chars = strcspn(bptr,gbl_hi);
         char save_char='\0';
         if(b_num_chars){
            save_char = bigheader[b_num_chars];
            bigheader[b_num_chars]='\0';
         }
         emit = 0;
         int matched = bin_search(bigheader, header_name_list, entrynum);
         if((matched != -1) ^ gbl_reject){ // (matches and NOT reject) OR (NOT matches AND reject) == matches XOR reject
             if(!gbl_reject){
                emitting=matched;
                if(emitlist[matched]){
                   (void) fprintf(stderr,"fastaselecth: at fasta header: %s\n",bptr);
                   insane("fastaselecth: fatal error: duplicate entry name in FASTA file");
                }
                emitlist[matched]=1;
                if(gbl_frag){
                   if(!group_name_list[emitting] || !strlen(group_name_list[emitting]))insane("fastaselecth: fatal error: -frac[ac] used but one or more selectors lack second field");
                   emitgroups[emitorder[emitting]]=group_name_list[emitting];
                }
             }
             emit=1;
             accumstring=NULL;
             tail=0;
             size=0;
             emitted++;
         }
         if(b_num_chars){
            bigheader[b_num_chars] = save_char;
         }
      }

      if(emit){
        if(gbl_reject){ //write immediately
           (void) fprintf(fout,"%s\n",bigstring);
        }
        else {
           size=size + strlen(bigstring) + 2;
           tbuf=malloc(size*sizeof(char));
           if(tbuf==NULL)insane("fastaselecth: fatal error: ran out of memory during processing");
           if(accumstring != NULL){
             strcpy(tbuf,accumstring);
             free(accumstring);
           }
           accumstring=tbuf;
           (void) sprintf(&accumstring[tail],"%s\n",bigstring);
           size--;
           tail=size;
        }
      }
      if(DONE)break;
   } /* end of reading loop */
   
   /*if some were not found, now is the time to say so*/
   
   if(!gbl_reject && (emitted <= entrynum - 1)){
     for(i=0;i<entrynum;i++){
        if(! emitlist[i]){
           (void)fprintf(stderr,"fastaselecth: %s: did not find selector: %s\n",(gbl_com ? "warning" : "fatal error"), header_name_list[i]);
        }
     }
     if(!gbl_com){
        exit(EXIT_FAILURE);
     }
   }
   
   /* may still have been accumulating one.  In worst case this was the first
      to be emitted, so all the other strings are still in memory */
   
   if(accumstring!=NULL){
      emitstrings[emitorder[emitting]]=accumstring;
   }

   /* force out anything left in emitstrings.  If there was a miss it may have stalled
      very early...*/
   for(lastemitted++; lastemitted < entrynum; lastemitted++){
     if(emitstrings[lastemitted]!=NULL){  /* next one in order is available to emit */
        if(gbl_frag && strcmp(last_group,emitgroups[lastemitted])){
           last_group = emitgroups[lastemitted];
           fclose(fout);
           sprintf(temp_name,gbl_out,last_group);
           if(gbl_frag == FRAG_APPEND){
               fout = fopen(temp_name,"a");
           }
           else if (gbl_frag == FRAG_NEW){
               FILE *fprobe = fopen(temp_name,"r");
               if(fprobe){
                  fprintf(stderr,"fastaselecth: fatal error: file name: %s\n",temp_name);
                  insane("fastaselecth: fatal error: -fragc mode output file already exists or noncontiguous group records");
               }
               else {
                  fout = fopen(temp_name,"w");
               }
           }
           if(!fout){
              fprintf(stderr,"fastaselecth: fatal error: file name: %s\n",temp_name);
              insane("fastaselecth: fatal error: could not open output file in -frag mode");
           }
        }
        (void) fprintf(fout,"%s",emitstrings[lastemitted]);
        free(emitstrings[lastemitted]); /* release memory */
     }
   } 
bye:

   /* clean up */
   fclose(fin);
   if(fout!=stdout){
      fclose(fout);
   }
   free(bigheader);
   free(bigstring);
   free(emitlist);
   free(emitorder);
   free(emitstrings);
   for(i=0;i<entrynum;i++){
      free(header_name_list[i]);
      if(gbl_frag){
         free(group_name_list[i]);
      }
   }
   free(emitgroups);  // all entries in emitgroups were pointers to a group_name_list entry
   free(header_name_list);
   if(gbl_frag){
      free(group_name_list);
   }
   free(gbl_hs);
   free(gbl_hi);
   
   fprintf(stderr,"fastaselecth: status: selectors: %d, records read: %llu, emitted: %llu\n",entrynum, records,emitted);
   
   exit(EXIT_SUCCESS);
}