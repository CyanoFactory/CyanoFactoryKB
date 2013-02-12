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
 * Bison Parser for TERMINATOR
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.terminator.c"
  #include "string-util.h"
  #include "terminator-parse.h"
  
  int terminator_records = 0;
  int terminator_errors  = 0;
  short first_terminator = 1; /* flag for first entry in file */
  struct terminator_entry current_terminator;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN
%token COMPONENT_OF_T RIGHT_END_POSITION_T LEFT_END_POSITION_T
%token CITATIONS_T COMMENT_T COMMON_NAME_T DBLINKS_T SYNONYMS_T UNIQUE_ID_T IGNORED_KEYWORD_T 

%%

terminatorfile   : terminatorrecords
                 | // no entries
                 ;

terminatorrecords : terminatorrecords terminatorrecord
                 | terminatorrecord
                 ;

terminatorrecord : terminatorslots ENDOFRECORD
                   { terminator_next_entry(); first_terminator = 0; }
                 ;

terminatorslots : terminatorslot
                 | terminatorslots terminatorslot
                 ;

terminatorslot  : unique_id
                 | citation
                 | comment
                 | common_name
                 | component_of
                 | dblinks
                 | right_end_position
                 | left_end_position
		 | synonyms
		 | ignored_keyword
		 | unknown_keyword
                 | error
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR rhs NEWLINE
                   { current_terminator.unique_id = $<string>3; }
                 ;

citation         : CITATIONS_T SEPARATOR rhs NEWLINE
                   { current_terminator.citations = stringlist_cons($<string>3, current_terminator.citations); }
                 | CITATIONS_T SEPARATOR NEWLINE  /* empty citation = a No-OP */
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_terminator.comments = stringlist_cons($<string>3, current_terminator.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

common_name      : COMMON_NAME_T SEPARATOR rhs NEWLINE
                   { current_terminator.common_name = $<string>3;  }
                 | COMMON_NAME_T SEPARATOR NEWLINE  /* empty common name = a No-OP */
                 ;

right_end_position : RIGHT_END_POSITION_T SEPARATOR id NEWLINE
                   { current_terminator.right_end_position = $<string>3;  }
                 | RIGHT_END_POSITION_T SEPARATOR NEWLINE  /* empty = a No-OP */
                 ;

left_end_position : LEFT_END_POSITION_T SEPARATOR id NEWLINE
                   { current_terminator.left_end_position = $<string>3;  }
                 | LEFT_END_POSITION_T SEPARATOR NEWLINE  /* empty = a No-OP */
                 ;

component_of     : COMPONENT_OF_T SEPARATOR rhs NEWLINE
                   { current_terminator.component_of = stringlist_cons($<string>3, current_terminator.component_of); }
                 | COMPONENT_OF_T SEPARATOR NEWLINE  /* empty component_of = a No-OP */
                 ;

dblinks          : DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id ignored_ids CLOSEPAREN NEWLINE
                 | DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id  CLOSEPAREN NEWLINE
                 ;

dblinks_id       : id
                   { current_terminator.dblinks_ids = stringlist_cons($<string>1, current_terminator.dblinks_ids); }
                 ;

dblinks_db       : id
                   { current_terminator.dblinks_dbs = stringlist_cons($<string>1, current_terminator.dblinks_dbs); }
                 ;

synonyms         : SYNONYMS_T SEPARATOR rhs NEWLINE
                   { current_terminator.synonyms = stringlist_cons($<string>3, current_terminator.synonyms); }
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
  printf("%s at `%s', line %d.\n", s, terminatortext, terminator_lineno);
  terminator_parse_error();
}

void
free_terminator_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
terminator_clear_entry(short freeflag) {
  /* Deallocate current_terminator and clear it in preparation for next terminator
   */
  current_terminator.load_error = 0;
  
  free_terminator_value(&current_terminator.unique_id, freeflag);
  free_terminator_value(&current_terminator.common_name, freeflag);
  free_terminator_value(&current_terminator.right_end_position, freeflag);
  free_terminator_value(&current_terminator.left_end_position, freeflag);
  
  if (freeflag && current_terminator.component_of)
    free_stringlist(current_terminator.component_of);
  current_terminator.component_of = NULL;
  
  if (freeflag && current_terminator.synonyms)
    free_stringlist(current_terminator.synonyms);
  current_terminator.synonyms = NULL;
  
  if (freeflag && current_terminator.citations) {
    free_stringlist(current_terminator.citations);
  }
  current_terminator.citations = NULL;
  
  if (freeflag && current_terminator.comments) {
    free_stringlist(current_terminator.comments);
  }
  current_terminator.comments = NULL;

  if (freeflag && current_terminator.dblinks_dbs) {
    free_stringlist(current_terminator.dblinks_dbs);
  }
  current_terminator.dblinks_dbs = NULL;  

  if (freeflag && current_terminator.dblinks_ids) {
    free_stringlist(current_terminator.dblinks_ids);
  }
  current_terminator.dblinks_ids = NULL;  
}
						   
void
terminator_next_entry(void) {

  /* Send this one to the database */
  terminator_load_entry(&current_terminator); 

  /* Then clear the structure for reuse */
  terminator_clear_entry(1); 
  terminator_records++;
}  

void
terminator_parse_error(void) {
  current_terminator.load_error = 1;
  terminator_errors++;
}
