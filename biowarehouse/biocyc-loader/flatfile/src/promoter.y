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
 * Bison Parser for PROMOTER
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.promoter.c"
  #include "string-util.h"
  #include "promoter-parse.h"
  
  int promoter_records = 0;
  int promoter_errors  = 0;
  short first_promoter = 1; /* flag for first entry in file */
  struct promoter_entry current_promoter;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN
%token COMPONENT_OF_T ABSOLUTE_PLUS1_POS_T
%token CITATIONS_T COMMENT_T COMMON_NAME_T DBLINKS_T SYNONYMS_T UNIQUE_ID_T IGNORED_KEYWORD_T 

%%

promoterfile     : promoterrecords
                 | // no entries
                 ;

promoterrecords  : promoterrecords promoterrecord
                 | promoterrecord
                 ;

promoterrecord   : promoterslots ENDOFRECORD
                   { promoter_next_entry(); first_promoter = 0; }
                 ;

promoterslots    : promoterslot
                 | promoterslots promoterslot
                 ;

promoterslot     : unique_id
                 | absolute_plus1_pos
                 | citation
                 | comment
                 | common_name
                 | component_of
                 | dblinks
		 | synonyms
		 | ignored_keyword
		 | unknown_keyword
                 | error
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR rhs NEWLINE
                   { current_promoter.unique_id = $<string>3; }
                 ;

citation         : CITATIONS_T SEPARATOR rhs NEWLINE
                   { current_promoter.citations = stringlist_cons($<string>3, current_promoter.citations); }
                 | CITATIONS_T SEPARATOR NEWLINE  /* empty citation = a No-OP */
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_promoter.comments = stringlist_cons($<string>3, current_promoter.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

common_name      : COMMON_NAME_T SEPARATOR rhs NEWLINE
                   { current_promoter.common_name = $<string>3;  }
                 | COMMON_NAME_T SEPARATOR NEWLINE  /* empty common name = a No-OP */
                 ;

absolute_plus1_pos : ABSOLUTE_PLUS1_POS_T SEPARATOR rhs NEWLINE
                   { current_promoter.absolute_plus1_pos = $<string>3;  }
                 | ABSOLUTE_PLUS1_POS_T SEPARATOR NEWLINE  /* empty common name = a No-OP */
                 ;

component_of     : COMPONENT_OF_T SEPARATOR rhs NEWLINE
                   { current_promoter.component_of = stringlist_cons($<string>3, current_promoter.component_of); }
                 | COMPONENT_OF_T SEPARATOR NEWLINE  /* empty component_of = a No-OP */
                 ;

dblinks          : DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id ignored_ids CLOSEPAREN NEWLINE
                 | DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id  CLOSEPAREN NEWLINE
                 ;

dblinks_id       : id
                   { current_promoter.dblinks_ids = stringlist_cons($<string>1, current_promoter.dblinks_ids); }
                 ;

dblinks_db       : id
                   { current_promoter.dblinks_dbs = stringlist_cons($<string>1, current_promoter.dblinks_dbs); }
                 ;

synonyms         : SYNONYMS_T SEPARATOR rhs NEWLINE
                   { current_promoter.synonyms = stringlist_cons($<string>3, current_promoter.synonyms); }
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
  printf("%s at `%s', line %d.\n", s, promotertext, promoter_lineno);
  promoter_parse_error();
}

void
free_promoter_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
promoter_clear_entry(short freeflag) {
  /* Deallocate current_promoter and clear it in preparation for next promoter
   */
  current_promoter.load_error = 0;
  
  free_promoter_value(&current_promoter.unique_id, freeflag);
  free_promoter_value(&current_promoter.common_name, freeflag);
  free_promoter_value(&current_promoter.absolute_plus1_pos, freeflag);
  
  if (freeflag && current_promoter.component_of)
    free_stringlist(current_promoter.component_of);
  current_promoter.component_of = NULL;
  
  if (freeflag && current_promoter.synonyms)
    free_stringlist(current_promoter.synonyms);
  current_promoter.synonyms = NULL;
  
  if (freeflag && current_promoter.citations) {
    free_stringlist(current_promoter.citations);
  }
  current_promoter.citations = NULL;
  
  if (freeflag && current_promoter.comments) {
    free_stringlist(current_promoter.comments);
  }
  current_promoter.comments = NULL;

  if (freeflag && current_promoter.dblinks_dbs) {
    free_stringlist(current_promoter.dblinks_dbs);
  }
  current_promoter.dblinks_dbs = NULL;  

  if (freeflag && current_promoter.dblinks_ids) {
    free_stringlist(current_promoter.dblinks_ids);
  }
  current_promoter.dblinks_ids = NULL;  
}
						   
void
promoter_next_entry(void) {

  /* Send this one to the database */
  promoter_load_entry(&current_promoter); 

  /* Then clear the structure for reuse */
  promoter_clear_entry(1); 
  promoter_records++;
}  

void
promoter_parse_error(void) {
  current_promoter.load_error = 1;
  promoter_errors++;
}
