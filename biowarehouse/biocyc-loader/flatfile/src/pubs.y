/* *******************************************************************
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See
 * the License for the specific language governing rights and
 * limitations under the License.
 *
 * The Original Code is the BioWarehouse.
 *
 * The Initial Developer of the Original Code is SRI International.
 * Portions created by SRI International are Copyright (C) 2004.
 * All Rights Reserved.
 ******************************************************************* */
/*
 * Bison Parser for PUBS
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.pubs.c"
  #include "string-util.h"
  #include "pubs-parse.h"
  
  int pubs_records = 0;
  int pubs_errors  = 0;
  short first_pubs = 1; /* flag for first entry in file */
  struct pubs_entry current_pubs;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN
%token IGNORED_KEYWORD_T
%token COMMENT_T AUTHORS_T MEDLINE_UID_T PUBMED_ID_T REFERENT_FRAME_T SOURCE_T TITLE_T URL_T YEAR_T UNIQUE_ID_T 

%%

pubsfile         : pubsrecords
                 | // no entries
                 ;

pubsrecords      : pubsrecords pubsrecord
                 | pubsrecord
                 ;

pubsrecord       : pubsslots ENDOFRECORD
                   { pubs_next_entry(); first_pubs = 0; }
                 ;

pubsslots        : pubsslot
                 | pubsslots pubsslot
                 ;

pubsslot         : unique_id
                 | author
                 | comment
                 | medline_uid
                 | pubmed_id
                 | referent_frame
		 | source
                 | title
                 | url
                 | year
		 | ignored_keyword
		 | unknown_keyword
                 | error
                 ;

author           : AUTHORS_T SEPARATOR rhs NEWLINE
                   /* !!! Retain order of authors */
                   { current_pubs.authors = stringlist_cons($<string>3, current_pubs.authors); }
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_pubs.comments = stringlist_cons($<string>3, current_pubs.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

medline_uid      : MEDLINE_UID_T SEPARATOR strval NEWLINE
                   { current_pubs.medline_uid = $<string>3; }
                 ;

pubmed_id        : PUBMED_ID_T SEPARATOR strval NEWLINE
                   { current_pubs.pubmed_id = $<string>3; }
                 ;

referent_frame   : REFERENT_FRAME_T SEPARATOR strval NEWLINE
                   { current_pubs.referent_frame = $<string>3; }
                 ;

source           : SOURCE_T SEPARATOR rhs NEWLINE
                   { current_pubs.source = $<string>3; }
                 ;

title            : TITLE_T SEPARATOR rhs NEWLINE
                   { current_pubs.title = $<string>3; }
                 ;

url              : URL_T SEPARATOR rhs NEWLINE
                   { current_pubs.url = $<string>3; }
                 ;

year             : YEAR_T SEPARATOR strval NEWLINE
                   { current_pubs.year = $<string>3; }
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR rhs NEWLINE
                   { current_pubs.unique_id = $<string>3; }
                 ;

unknown_keyword  : id SEPARATOR rhs NEWLINE
                   { wh_add_undefined($<string>1, "Warning: undefined attribute %s\n"); }
		 | id SEPARATOR NEWLINE
                   { wh_add_undefined($<string>1, "Warning: undefined attribute %s\n");}
                 ;


/*** Ignored fields ***/

ignored_keyword  : IGNORED_KEYWORD_T SEPARATOR rhs NEWLINE
                 | IGNORED_KEYWORD_T SEPARATOR NEWLINE
                 | IGNORED_KEYWORD_T NEWLINE
                 ;


/* Tokens */

strval           : WORD
                   { $<string>$ = lexedword;  }
                 ;

id               : WORD
                   { $<string>$ = lexedword;  }
                 ;

rhs              : RHS
                   { $<string>$ = lexedword;  }
                 ;

%%

int
yyerror(char *s)
{
  printf("%s at `%s', line %d.\n", s, pubstext, pubs_lineno);
  pubs_parse_error();
}

void
free_pubs_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
pubs_clear_entry(short freeflag) {
  /* Deallocate current_pubs and clear it in preparation for next pubs
   */
  current_pubs.load_error = 0;
  free_pubs_value(&current_pubs.unique_id, freeflag);
  free_pubs_value(&current_pubs.referent_frame, freeflag);
  free_pubs_value(&current_pubs.medline_uid, freeflag);
  free_pubs_value(&current_pubs.pubmed_id, freeflag);
  free_pubs_value(&current_pubs.source, freeflag);
  free_pubs_value(&current_pubs.title, freeflag);
  free_pubs_value(&current_pubs.url, freeflag);
  free_pubs_value(&current_pubs.year, freeflag);

  /* Free all lists of items */
  if (freeflag && current_pubs.authors)
    free_stringlist(current_pubs.authors);
  current_pubs.authors = NULL;
  
  if (freeflag && current_pubs.comments) {
    free_stringlist(current_pubs.comments);
  }
  current_pubs.comments = NULL;
}
						   
void
pubs_next_entry(void) {

  /* Send this one to the database */
  pubs_load_entry(&current_pubs); 

  /* Then clear the structure for reuse */
  pubs_clear_entry(1); 
  pubs_records++;
}  

void
pubs_parse_error(void) {
  current_pubs.load_error = 1;
  pubs_errors++;
}
