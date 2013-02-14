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
 * Bison Parser for TRANSUNIT
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.transunit.c"
  #include "string-util.h"
  #include "transunit-parse.h"
  
  int transunit_records = 0;
  int transunit_errors  = 0;
  short first_transunit = 1; /* flag for first entry in file */
  struct transunit_entry current_transunit;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN
%token COMPONENTS_T 
%token CITATIONS_T COMMENT_T COMMON_NAME_T DBLINKS_T SYNONYMS_T UNIQUE_ID_T IGNORED_KEYWORD_T 

%%

transunitfile    : transunitrecords
                 | // no entries
                 ;

transunitrecords  : transunitrecords transunitrecord
                 | transunitrecord
                 ;

transunitrecord   : transunitslots ENDOFRECORD
                   { transunit_next_entry(); first_transunit = 0; }
                 ;

transunitslots    : transunitslot
                 | transunitslots transunitslot
                 ;

transunitslot     : unique_id
                 | citation
                 | comment
                 | common_name
                 | component
                 | dblinks
		 | synonyms
		 | ignored_keyword
		 | unknown_keyword
                 | error
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR rhs NEWLINE
                   { current_transunit.unique_id = $<string>3; }
                 ;

citation         : CITATIONS_T SEPARATOR rhs NEWLINE
                   { current_transunit.citations = stringlist_cons($<string>3, current_transunit.citations); }
                 | CITATIONS_T SEPARATOR NEWLINE  /* empty citation = a No-OP */
                 ;

component        : COMPONENTS_T SEPARATOR id NEWLINE  /* defines a subunit */
                   { current_transunit.components = stringlist_cons($<string>3, current_transunit.components); }
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_transunit.comments = stringlist_cons($<string>3, current_transunit.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

common_name      : COMMON_NAME_T SEPARATOR rhs NEWLINE
                   { current_transunit.common_name = $<string>3;  }
                 | COMMON_NAME_T SEPARATOR NEWLINE  /* empty common name = a No-OP */
                 ;

dblinks          : DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id ignored_ids CLOSEPAREN NEWLINE
                 | DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id  CLOSEPAREN NEWLINE
                 ;

dblinks_id       : id
                   { current_transunit.dblinks_ids = stringlist_cons($<string>1, current_transunit.dblinks_ids); }
                 ;

dblinks_db       : id
                   { current_transunit.dblinks_dbs = stringlist_cons($<string>1, current_transunit.dblinks_dbs); }
                 ;

synonyms         : SYNONYMS_T SEPARATOR rhs NEWLINE
                   { current_transunit.synonyms = stringlist_cons($<string>3, current_transunit.synonyms); }
                 | SYNONYMS_T SEPARATOR NEWLINE  /* empty synonyms = a No-OP */
                 ;

unknown_keyword  : id SEPARATOR rhs NEWLINE
                   { wh_add_undefined($<string>1, "Warning: undefined attribute %s\n"); }
		 | id SEPARATOR NEWLINE
                   { wh_add_undefined($<string>1, "Warning: undefined attribute %s\n");}
                 ;


/*** Ignored fields ***/

ignored_ids      : ignored_ids id
                 | id
                 ;

ignored_keyword  : IGNORED_KEYWORD_T SEPARATOR rhs NEWLINE
                 | IGNORED_KEYWORD_T SEPARATOR NEWLINE
                 | IGNORED_KEYWORD_T NEWLINE
                 ;


/* Tokens */

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
  printf("%s at `%s', line %d.\n", s, transunittext, transunit_lineno);
  transunit_parse_error();
}

void
free_transunit_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
transunit_clear_entry(short freeflag) {
  /* Deallocate current_transunit and clear it in preparation for next transunit
   */
  current_transunit.load_error = 0;
  
  free_transunit_value(&current_transunit.unique_id, freeflag);
  free_transunit_value(&current_transunit.common_name, freeflag);
  
  if (freeflag && current_transunit.components)
    free_stringlist(current_transunit.components);
  current_transunit.components = NULL;
  
  if (freeflag && current_transunit.synonyms)
    free_stringlist(current_transunit.synonyms);
  current_transunit.synonyms = NULL;
  
  if (freeflag && current_transunit.citations) {
    free_stringlist(current_transunit.citations);
  }
  current_transunit.citations = NULL;
  
  if (freeflag && current_transunit.comments) {
    free_stringlist(current_transunit.comments);
  }
  current_transunit.comments = NULL;

  if (freeflag && current_transunit.dblinks_dbs) {
    free_stringlist(current_transunit.dblinks_dbs);
  }
  current_transunit.dblinks_dbs = NULL;  

  if (freeflag && current_transunit.dblinks_ids) {
    free_stringlist(current_transunit.dblinks_ids);
  }
  current_transunit.dblinks_ids = NULL;  
}
						   
void
transunit_next_entry(void) {

  /* Send this one to the database */
  transunit_load_entry(&current_transunit); 

  /* Then clear the structure for reuse */
  transunit_clear_entry(1); 
  transunit_records++;
}  

void
transunit_parse_error(void) {
  current_transunit.load_error = 1;
  transunit_errors++;
}
