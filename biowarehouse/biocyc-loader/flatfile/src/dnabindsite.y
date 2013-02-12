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
 * Bison Parser for DNABINDSITE
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.dnabindsite.c"
  #include "string-util.h"
  #include "dnabindsite-parse.h"
  
  int dnabindsite_records = 0;
  int dnabindsite_errors  = 0;
  short first_dnabindsite = 1; /* flag for first entry in file */
  struct dnabindsite_entry current_dnabindsite;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN
%token COMPONENT_OF_T ABS_CENTER_POS_T REGULATED_PROMOTER_T
%token CITATIONS_T COMMENT_T COMMON_NAME_T DBLINKS_T SYNONYMS_T UNIQUE_ID_T IGNORED_KEYWORD_T 

%%

dnabindsitefile  : dnabindsiterecords
                 | // no entries
                 ;

dnabindsiterecords : dnabindsiterecords dnabindsiterecord
                 | dnabindsiterecord
                 ;

dnabindsiterecord : dnabindsiteslots ENDOFRECORD
                   { dnabindsite_next_entry(); first_dnabindsite = 0; }
                 ;

dnabindsiteslots : dnabindsiteslot
                 | dnabindsiteslots dnabindsiteslot
                 ;

dnabindsiteslot  : unique_id
                 | citation
                 | comment
                 | common_name
                 | component_of
                 | dblinks
                 | abs_center_pos
                 | regulated_promoter
		 | synonyms
		 | ignored_keyword
		 | unknown_keyword
                 | error
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR rhs NEWLINE
                   { current_dnabindsite.unique_id = $<string>3; }
                 ;

citation         : CITATIONS_T SEPARATOR rhs NEWLINE
                   { current_dnabindsite.citations = stringlist_cons($<string>3, current_dnabindsite.citations); }
                 | CITATIONS_T SEPARATOR NEWLINE  /* empty citation = a No-OP */
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_dnabindsite.comments = stringlist_cons($<string>3, current_dnabindsite.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

common_name      : COMMON_NAME_T SEPARATOR rhs NEWLINE
                   { current_dnabindsite.common_name = $<string>3;  }
                 | COMMON_NAME_T SEPARATOR NEWLINE  /* empty common name = a No-OP */
                 ;

abs_center_pos : ABS_CENTER_POS_T SEPARATOR id NEWLINE
                   { current_dnabindsite.abs_center_pos = $<string>3;  }
                 | ABS_CENTER_POS_T SEPARATOR NEWLINE  /* empty = a No-OP */
                 ;

regulated_promoter : REGULATED_PROMOTER_T SEPARATOR rhs NEWLINE
                   { current_dnabindsite.regulated_promoter = $<string>3;  }
                 | REGULATED_PROMOTER_T SEPARATOR NEWLINE  /* empty = a No-OP */
                 ;

component_of     : COMPONENT_OF_T SEPARATOR rhs NEWLINE
                   { current_dnabindsite.component_of = stringlist_cons($<string>3, current_dnabindsite.component_of); }
                 | COMPONENT_OF_T SEPARATOR NEWLINE  /* empty component_of = a No-OP */
                 ;

dblinks          : DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id ignored_ids CLOSEPAREN NEWLINE
                 | DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id  CLOSEPAREN NEWLINE
                 ;

dblinks_id       : id
                   { current_dnabindsite.dblinks_ids = stringlist_cons($<string>1, current_dnabindsite.dblinks_ids); }
                 ;

dblinks_db       : id
                   { current_dnabindsite.dblinks_dbs = stringlist_cons($<string>1, current_dnabindsite.dblinks_dbs); }
                 ;

synonyms         : SYNONYMS_T SEPARATOR rhs NEWLINE
                   { current_dnabindsite.synonyms = stringlist_cons($<string>3, current_dnabindsite.synonyms); }
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
  printf("%s at `%s', line %d.\n", s, dnabindsitetext, dnabindsite_lineno);
  dnabindsite_parse_error();
}

void
free_dnabindsite_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
dnabindsite_clear_entry(short freeflag) {
  /* Deallocate current_dnabindsite and clear it in preparation for next dnabindsite
   */
  current_dnabindsite.load_error = 0;
  
  free_dnabindsite_value(&current_dnabindsite.unique_id, freeflag);
  free_dnabindsite_value(&current_dnabindsite.common_name, freeflag);
  free_dnabindsite_value(&current_dnabindsite.abs_center_pos, freeflag);
  free_dnabindsite_value(&current_dnabindsite.regulated_promoter, freeflag);
  
  if (freeflag && current_dnabindsite.component_of)
    free_stringlist(current_dnabindsite.component_of);
  current_dnabindsite.component_of = NULL;
  
  if (freeflag && current_dnabindsite.synonyms)
    free_stringlist(current_dnabindsite.synonyms);
  current_dnabindsite.synonyms = NULL;
  
  if (freeflag && current_dnabindsite.citations) {
    free_stringlist(current_dnabindsite.citations);
  }
  current_dnabindsite.citations = NULL;
  
  if (freeflag && current_dnabindsite.comments) {
    free_stringlist(current_dnabindsite.comments);
  }
  current_dnabindsite.comments = NULL;

  if (freeflag && current_dnabindsite.dblinks_dbs) {
    free_stringlist(current_dnabindsite.dblinks_dbs);
  }
  current_dnabindsite.dblinks_dbs = NULL;  

  if (freeflag && current_dnabindsite.dblinks_ids) {
    free_stringlist(current_dnabindsite.dblinks_ids);
  }
  current_dnabindsite.dblinks_ids = NULL;  
}
						   
void
dnabindsite_next_entry(void) {

  /* Send this one to the database */
  dnabindsite_load_entry(&current_dnabindsite); 

  /* Then clear the structure for reuse */
  dnabindsite_clear_entry(1); 
  dnabindsite_records++;
}  

void
dnabindsite_parse_error(void) {
  current_dnabindsite.load_error = 1;
  dnabindsite_errors++;
}
