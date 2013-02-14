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
 * Bison Parser for REACTION
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.reaction.c"
  #include "string-util.h"
  #include "reaction-parse.h"
  
  int reaction_records = 0;
  int reaction_errors  = 0;
  short first_reaction = 1; /* flag for first entry in file */
  struct reaction_entry current_reaction;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN
%token CITATIONS_T COMMENT_T COMMON_NAME_T COEFFICIENT_T DBLINKS_T SYNONYMS_T SPONTANEOUS_FLAG_T OFFICIAL_EC_FLAG_T EC_NUMBER_T RIGHT_T LEFT_T DELTA_G0_T UNIQUE_ID_T IGNORED_KEYWORD_T

%%

reactionfile     : reactionrecords
                 | // no entries
                 ;

reactionrecords  : reactionrecords reactionrecord
                 | reactionrecord
                 ;

reactionrecord   : reactionslots ENDOFRECORD
                   { reaction_next_entry(); first_reaction = 0; }
                 ;

reactionslots    : reactionslot
                 | reactionslots reactionslot
                 ;

reactionslot :     unique_id
                 | citation
                 | coefficient
                 | comment
                 | common_name
                 | delta_g0
                 | dblinks
                 | ec_flag
                 | ec_number
                 | left
                 | right
		 | spontaneous
		 | synonyms
                 | ignored_keyword
                 | unknown_keyword
                 | error
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR id NEWLINE
                   { current_reaction.unique_id = $<string>3; }
                 ;

citation         : CITATIONS_T SEPARATOR rhs NEWLINE
                   { current_reaction.citations = stringlist_cons($<string>3, current_reaction.citations); }
                 | CITATIONS_T SEPARATOR NEWLINE  /* empty citation = a No-OP */
                 ;

coefficient      : COEFFICIENT_T SEPARATOR strval NEWLINE
                   { if (current_reaction.side_last_parsed == LEFT)
                         current_reaction.left_coeffs->string = $<string>3;
                     else if (current_reaction.side_last_parsed == RIGHT)
                         current_reaction.right_coeffs->string = $<string>3;
                     else fatal_error("Internal error"); }
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_reaction.comments = stringlist_cons($<string>3, current_reaction.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

common_name      : COMMON_NAME_T SEPARATOR rhs NEWLINE
                   { current_reaction.synonyms = stringlist_cons($<string>3, current_reaction.synonyms); }
                 ;

dblinks          : DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id ignored_ids CLOSEPAREN NEWLINE
                 | DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id  CLOSEPAREN NEWLINE
                 ;

dblinks_id       : id
                   { current_reaction.dblinks_ids = stringlist_cons($<string>1, current_reaction.dblinks_ids); }
                 ;

dblinks_db       : id
                   { current_reaction.dblinks_dbs = stringlist_cons($<string>1, current_reaction.dblinks_dbs); }
                 ;

ec_number        : EC_NUMBER_T SEPARATOR id NEWLINE
                   { current_reaction.ec_number = $<string>3; }
                 ;

ec_flag          : OFFICIAL_EC_FLAG_T SEPARATOR id NEWLINE
                   { current_reaction.official_ec_flag = $<string>3; }
                 ;

delta_g0         : DELTA_G0_T SEPARATOR strval NEWLINE
                   { current_reaction.delta_g0 = $<string>3; }
                 ;

left             : LEFT_T SEPARATOR rhs NEWLINE
                   { current_reaction.side_last_parsed = LEFT;
		     current_reaction.left = stringlist_cons($<string>3, current_reaction.left);
		     current_reaction.left_coeffs = stringlist_cons(strdup("1"), current_reaction.left_coeffs); }
                 ;

right            : RIGHT_T SEPARATOR rhs NEWLINE
                   { current_reaction.side_last_parsed = RIGHT;
                     current_reaction.right = stringlist_cons($<string>3, current_reaction.right); 
		     current_reaction.right_coeffs = stringlist_cons(strdup("1"), current_reaction.right_coeffs); }
                 ;

spontaneous      : SPONTANEOUS_FLAG_T SEPARATOR strval NEWLINE
                   { current_reaction.spontaneous = $<string>3; }
                 ;

synonyms         : SYNONYMS_T SEPARATOR rhs NEWLINE
                   { current_reaction.synonyms = stringlist_cons($<string>3, current_reaction.synonyms); }
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
  printf("%s at `%s', line %d.\n", s, reactiontext, reaction_lineno);
  reaction_parse_error();
}

void
free_reaction_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
reaction_clear_entry(short freeflag) {
  /* Deallocate current_reaction and clear it in preparation for next reaction
   */
  current_reaction.load_error = 0;
  
  free_reaction_value(&current_reaction.unique_id, freeflag);
  free_reaction_value(&current_reaction.ec_number, freeflag);
  free_reaction_value(&current_reaction.delta_g0, freeflag);
  free_reaction_value(&current_reaction.spontaneous, freeflag);
  
  if (freeflag && current_reaction.synonyms)
    free_stringlist(current_reaction.synonyms);
  current_reaction.synonyms = NULL;
  if (freeflag && current_reaction.left)
    free_stringlist(current_reaction.left);
  current_reaction.left = NULL;
  if (freeflag && current_reaction.left_coeffs)
    free_stringlist(current_reaction.left_coeffs);
  current_reaction.left_coeffs = NULL;
  if (freeflag && current_reaction.right)
    free_stringlist(current_reaction.right);
  current_reaction.right = NULL;
  if (freeflag && current_reaction.right_coeffs)
    free_stringlist(current_reaction.right_coeffs);
  current_reaction.right_coeffs = NULL;
  
  if (freeflag && current_reaction.citations) {
    free_stringlist(current_reaction.citations);
  }
  current_reaction.citations = NULL;
  
  if (freeflag && current_reaction.comments) {
    free_stringlist(current_reaction.comments);
  }
  current_reaction.comments = NULL;

  if (freeflag && current_reaction.dblinks_dbs) {
    free_stringlist(current_reaction.dblinks_dbs);
  }
  current_reaction.dblinks_dbs = NULL;  

  if (freeflag && current_reaction.dblinks_ids) {
    free_stringlist(current_reaction.dblinks_ids);
  }
  current_reaction.dblinks_ids = NULL;  
}
						   
void
reaction_next_entry(void) {

  /* Send this one to the database */
  reaction_load_entry(&current_reaction); 

  /* Then clear the structure for reuse */
  reaction_clear_entry(1); 
  reaction_records++;
}  

void
reaction_parse_error(void) {
  current_reaction.load_error = 1;
  reaction_errors++;
}
