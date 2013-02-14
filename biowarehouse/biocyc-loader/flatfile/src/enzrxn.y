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
 * Bison Parser for ENZRXN
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.enzrxn.c"
  #include "string-util.h"
  #include "enzrxn-parse.h"
  
  int enzrxn_records = 0;
  int enzrxn_errors  = 0;
  short first_enzrxn = 1; /* flag for first entry in file */
  struct enzrxn_entry current_enzrxn;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN
%token CITATIONS_T COMMON_NAME_T COMMENT_T DBLINKS_T SYNONYMS_T UNIQUE_ID_T IGNORED_KEYWORD_T
%token COMPLEX_T ENZYME_T REACTION_T REACTION_DIRECTION_T
%token COFACTOR_BINDING_COMMENT_T COFACTORS_T COFACTORS_OR_PROSTHETIC_GROUPS_T PROSTHETIC_GROUPS_T
%token ALTERNATIVE_COFACTORS_T ALTERNATIVE_SUBSTRATES_T

%%

enzrxnfile       : enzrxnrecords
                 | // no entries
                 ;

enzrxnrecords    : enzrxnrecords enzrxnrecord
                 | enzrxnrecord
                 ;

enzrxnrecord     : enzrxnslots ENDOFRECORD
                   { enzrxn_next_entry(); first_enzrxn = 0; }
                 ;

enzrxnslots      : enzrxnslot
                 | enzrxnslots enzrxnslot
                 ;

enzrxnslot       : unique_id
                 | citation
                 | common_name
                 | comment
                 | dblinks
		 | synonyms
                 | complex
                 | enzyme
                 | reaction
                 | reaction_direction
                 | cofactor_binding_comment
                 | cofactors
                 | cofactors_or_prosthetic_groups
                 | prosthetic_groups
                 | alt_cofactors 
                 | alt_substrates 
		 | ignored_keyword
		 | unknown_keyword
                 | error
                   { yyerror("Unrecognized attribute"); }
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR id NEWLINE
                   { current_enzrxn.unique_id = $<string>3; }
                 ;

cofactors        : COFACTORS_T SEPARATOR rhs NEWLINE
                   { current_enzrxn.cofactors = stringlist_cons($<string>3, current_enzrxn.cofactors); }
                 ;

prosthetic_groups : PROSTHETIC_GROUPS_T SEPARATOR rhs NEWLINE
                   { current_enzrxn.cofactor_prosthetics = stringlist_cons($<string>3, current_enzrxn.cofactor_prosthetics); }
                 ;

cofactors_or_prosthetic_groups : COFACTORS_OR_PROSTHETIC_GROUPS_T SEPARATOR rhs NEWLINE
                   { current_enzrxn.cofactor_unknowns = stringlist_cons($<string>3, current_enzrxn.cofactor_unknowns); }
                 ;

cofactor_binding_comment : COFACTOR_BINDING_COMMENT_T SEPARATOR rhs NEWLINE
                   { current_enzrxn.cofactor_comments = stringlist_cons($<string>3, current_enzrxn.cofactor_comments); }
                 | COFACTOR_BINDING_COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

alt_cofactors    : ALTERNATIVE_COFACTORS_T SEPARATOR OPENPAREN alt_cofactor_ids CLOSEPAREN NEWLINE
                 ;

alt_cofactor_ids : alt_cofactor_ids alt_cofactor_id
                 | alt_cofactor_id
                 ;

alt_cofactor_id  : id
                   { current_enzrxn.alt_cofactors = stringlist_cons($<string>1, current_enzrxn.alt_cofactors); }
                 ;

alt_substrates   : ALTERNATIVE_SUBSTRATES_T SEPARATOR OPENPAREN alt_substrate_ids CLOSEPAREN NEWLINE
                 ;

alt_substrate_ids : alt_substrate_ids alt_substrate_id
                 | alt_substrate_id
                 ;

alt_substrate_id : id
                   { current_enzrxn.alt_substrates = stringlist_cons($<string>1, current_enzrxn.alt_substrates); }
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_enzrxn.comments = stringlist_cons($<string>3, current_enzrxn.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

common_name      : COMMON_NAME_T SEPARATOR rhs NEWLINE
                   { current_enzrxn.common_name = $<string>3; }
                 ;

complex          : COMPLEX_T SEPARATOR id NEWLINE
                   { current_enzrxn.complex = $<string>3; }
                 ;

dblinks          : DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id ignored_ids CLOSEPAREN NEWLINE
                 | DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id  CLOSEPAREN NEWLINE
                 ;

dblinks_id       : id
                   { current_enzrxn.dblinks_ids = stringlist_cons($<string>1, current_enzrxn.dblinks_ids); }
                 ;

dblinks_db       : id
                   { current_enzrxn.dblinks_dbs = stringlist_cons($<string>1, current_enzrxn.dblinks_dbs); }
                 ;

enzyme           : ENZYME_T SEPARATOR id NEWLINE
                   { current_enzrxn.enzyme = $<string>3; }
                 ;

reaction         : REACTION_T SEPARATOR id NEWLINE
                   { current_enzrxn.reaction = $<string>3; }
                 ;

reaction_direction  : REACTION_DIRECTION_T SEPARATOR id NEWLINE
                   { current_enzrxn.reaction_direction = $<string>3; }
                 ;

citation         : CITATIONS_T SEPARATOR rhs NEWLINE
                   { current_enzrxn.citations = stringlist_cons($<string>3, current_enzrxn.citations); }
                 | CITATIONS_T SEPARATOR NEWLINE  /* empty citation = a No-OP */
                 ;

synonyms         : SYNONYMS_T SEPARATOR rhs NEWLINE
                   { current_enzrxn.synonyms = stringlist_cons($<string>3, current_enzrxn.synonyms); }
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

/* strval           : WORD */
/*                    { $<string>$ = lexedword;  } */
/*                  ; */

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
  printf("%s at `%s', line %d.\n", s, enzrxntext, enzrxn_lineno);
  enzrxn_parse_error();
}

void
free_enzrxn_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
enzrxn_clear_entry(short freeflag) {
  /* Deallocate current_enzrxn and clear it in preparation for next enzrxn
   */
  current_enzrxn.load_error = 0;
  
  free_enzrxn_value(&current_enzrxn.unique_id, freeflag);
  free_enzrxn_value(&current_enzrxn.common_name, freeflag);
  free_enzrxn_value(&current_enzrxn.complex, freeflag);
  free_enzrxn_value(&current_enzrxn.enzyme, freeflag);
  free_enzrxn_value(&current_enzrxn.reaction, freeflag);
  free_enzrxn_value(&current_enzrxn.reaction_direction, freeflag);
  
  if (freeflag && current_enzrxn.cofactors) {
    free_stringlist(current_enzrxn.cofactors);
  }
  current_enzrxn.cofactors = NULL;
  
  if (freeflag && current_enzrxn.cofactor_prosthetics) {
    free_stringlist(current_enzrxn.cofactor_prosthetics);
  }
  current_enzrxn.cofactor_prosthetics = NULL;
  
  if (freeflag && current_enzrxn.cofactor_unknowns) {
    free_stringlist(current_enzrxn.cofactor_unknowns);
  }
  current_enzrxn.cofactor_unknowns = NULL;
  
  if (freeflag && current_enzrxn.cofactor_comments) {
    free_stringlist(current_enzrxn.cofactor_comments);
  }
  current_enzrxn.cofactor_comments = NULL;
  
  if (freeflag && current_enzrxn.alt_cofactors) {
    free_stringlist(current_enzrxn.alt_cofactors);
  }
  current_enzrxn.alt_cofactors = NULL;  

  if (freeflag && current_enzrxn.alt_substrates) {
    free_stringlist(current_enzrxn.alt_substrates);
  }
  current_enzrxn.alt_substrates = NULL;  

  if (freeflag && current_enzrxn.dblinks_dbs) {
    free_stringlist(current_enzrxn.dblinks_dbs);
  }
  current_enzrxn.dblinks_dbs = NULL;  

  if (freeflag && current_enzrxn.dblinks_ids) {
    free_stringlist(current_enzrxn.dblinks_ids);
  }
  current_enzrxn.dblinks_ids = NULL;  

  if (freeflag && current_enzrxn.synonyms)
    free_stringlist(current_enzrxn.synonyms);
  current_enzrxn.synonyms = NULL;
  
  if (freeflag && current_enzrxn.citations) {
    free_stringlist(current_enzrxn.citations);
  }
  current_enzrxn.citations = NULL;

  if (freeflag && current_enzrxn.comments) {
    free_stringlist(current_enzrxn.comments);
  }
  current_enzrxn.comments = NULL;
}
						   
void
enzrxn_next_entry(void) {

  /* Send this one to the database */
  enzrxn_load_entry(&current_enzrxn); 

  /* Then clear the structure for reuse */
  enzrxn_clear_entry(1); 
  enzrxn_records++;
}  

void
enzrxn_parse_error(void) {
  current_enzrxn.load_error = 1;
  enzrxn_errors++;
}
