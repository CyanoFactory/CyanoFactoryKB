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
 * Bison Parser for REGULATION
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.regulation.c"
  #include "string-util.h"
  #include "regulation-parse.h"
  
  int regulation_records = 0;
  int regulation_errors  = 0;
  short first_regulation = 1; /* flag for first entry in file */
  struct regulation_entry current_regulation;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN
%token MODE_T MECHANISM_T PHYSIOLOGICALLY_RELEVANT_T REGULATED_ENTITY_T REGULATOR_T TYPES_T
%token CITATIONS_T COMMENT_T DBLINKS_T SYNONYMS_T UNIQUE_ID_T IGNORED_KEYWORD_T 

%%

regulationfile   : regulationrecords
                 | // no entries
                 ;

regulationrecords  : regulationrecords regulationrecord
                 | regulationrecord
                 ;

regulationrecord   : regulationslots ENDOFRECORD
                   { regulation_next_entry(); first_regulation = 0; }
                 ;

regulationslots    : regulationslot
                 | regulationslots regulationslot
                 ;

regulationslot     : unique_id
                 | citation
                 | comment
                 | dblinks
                 | mechanism
                 | mode
                 | physiologically_relevant
                 | regulated_entity
                 | regulator
                 | types
		 | synonyms
		 | ignored_keyword
		 | unknown_keyword
                 | error
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR rhs NEWLINE
                   { current_regulation.unique_id = $<string>3; }
                 ;

citation         : CITATIONS_T SEPARATOR rhs NEWLINE
                   { current_regulation.citations = stringlist_cons($<string>3, current_regulation.citations); }
                 | CITATIONS_T SEPARATOR NEWLINE  /* empty citation = a No-OP */
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_regulation.comments = stringlist_cons($<string>3, current_regulation.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

mechanism        : MECHANISM_T SEPARATOR rhs NEWLINE
                   { current_regulation.mechanism = $<string>3;  }
                 | MECHANISM_T SEPARATOR NEWLINE  /* empty = a No-OP */
                 ;

mode             : MODE_T SEPARATOR rhs NEWLINE
                   { current_regulation.mode = $<string>3;  }
                 | MODE_T SEPARATOR NEWLINE  /* empty = a No-OP */
                 ;

physiologically_relevant : PHYSIOLOGICALLY_RELEVANT_T SEPARATOR rhs NEWLINE
                   { current_regulation.physiologically_relevant = $<string>3;  }
                 | PHYSIOLOGICALLY_RELEVANT_T SEPARATOR NEWLINE  /* empty = a No-OP */
                 ;

regulated_entity : REGULATED_ENTITY_T SEPARATOR rhs NEWLINE
                   { current_regulation.regulated_entity = $<string>3;  }
                 | REGULATED_ENTITY_T SEPARATOR NEWLINE  /* empty = a No-OP */
                 ;

regulator        : REGULATOR_T SEPARATOR rhs NEWLINE
                   { current_regulation.regulator = $<string>3;  }
                 | REGULATOR_T SEPARATOR NEWLINE  /* empty = a No-OP */
                 ;

types            : TYPES_T SEPARATOR rhs NEWLINE
                   { current_regulation.type = $<string>3;  }
                 | TYPES_T SEPARATOR NEWLINE  /* empty = a No-OP */
                 ;

dblinks          : DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id ignored_ids CLOSEPAREN NEWLINE
                 | DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id  CLOSEPAREN NEWLINE
                 ;

dblinks_id       : id
                   { current_regulation.dblinks_ids = stringlist_cons($<string>1, current_regulation.dblinks_ids); }
                 ;

dblinks_db       : id
                   { current_regulation.dblinks_dbs = stringlist_cons($<string>1, current_regulation.dblinks_dbs); }
                 ;

synonyms         : SYNONYMS_T SEPARATOR rhs NEWLINE
                   { current_regulation.synonyms = stringlist_cons($<string>3, current_regulation.synonyms); }
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
  printf("%s at `%s', line %d.\n", s, regulationtext, regulation_lineno);
  regulation_parse_error();
}

void
free_regulation_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
regulation_clear_entry(short freeflag) {
  /* Deallocate current_regulation and clear it in preparation for next regulation
   */
  current_regulation.load_error = 0;
  
  free_regulation_value(&current_regulation.unique_id, freeflag);
  free_regulation_value(&current_regulation.mechanism, freeflag);
  free_regulation_value(&current_regulation.mode, freeflag);
  free_regulation_value(&current_regulation.physiologically_relevant, freeflag);
  free_regulation_value(&current_regulation.regulated_entity, freeflag);
  free_regulation_value(&current_regulation.regulator, freeflag);
  free_regulation_value(&current_regulation.type, freeflag);
  
  if (freeflag && current_regulation.synonyms)
    free_stringlist(current_regulation.synonyms);
  current_regulation.synonyms = NULL;
  
  if (freeflag && current_regulation.citations) {
    free_stringlist(current_regulation.citations);
  }
  current_regulation.citations = NULL;
  
  if (freeflag && current_regulation.comments) {
    free_stringlist(current_regulation.comments);
  }
  current_regulation.comments = NULL;

  if (freeflag && current_regulation.dblinks_dbs) {
    free_stringlist(current_regulation.dblinks_dbs);
  }
  current_regulation.dblinks_dbs = NULL;  

  if (freeflag && current_regulation.dblinks_ids) {
    free_stringlist(current_regulation.dblinks_ids);
  }
  current_regulation.dblinks_ids = NULL;  
}
						   
void
regulation_next_entry(void) {

  /* Send this one to the database */
  regulation_load_entry(&current_regulation); 

  /* Then clear the structure for reuse */
  regulation_clear_entry(1); 
  regulation_records++;
}  

void
regulation_parse_error(void) {
  current_regulation.load_error = 1;
  regulation_errors++;
}
