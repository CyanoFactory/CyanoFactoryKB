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
 * Bison Parser for GENE
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.gene.c"
  #include "string-util.h"
  #include "gene-parse.h"
  
  int gene_records = 0;
  int gene_errors  = 0;
  short first_gene = 1; /* flag for first entry in file */
  struct gene_entry current_gene;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN
%token CITATIONS_T COMMENT_T COMMON_NAME_T DBLINKS_T SYNONYMS_T UNIQUE_ID_T TYPES_T IGNORED_KEYWORD_T
%token EVIDENCE_T INTERRUPTED_FLAG_T LEFT_END_POSITION_T PRODUCT_T RIGHT_END_POSITION_T TRANSCRIPTION_DIRECTION_T COMPONENT_OF_T

%%

genefile         : generecords
                 | // no entries
                 ;

generecords      : generecords generecord
                 | generecord
                 ;

generecord       : geneslots ENDOFRECORD
                   { gene_next_entry(); first_gene = 0; }
                 ;

geneslots        : geneslot
                 | geneslots geneslot
                 ;

geneslot         : unique_id
                 | citation
                 | comment
                 | common_name
                 | component_of
                 | dblinks
		 | synonyms
                 | evidence
                 | interrupted
                 | left_end_position
                 | right_end_position
                 | product
                 | transcription_direction
                 | types
		 | ignored_keyword
		 | unknown_keyword
                 | error
                 ;

citation         : CITATIONS_T SEPARATOR rhs NEWLINE
                   { current_gene.citations = stringlist_cons($<string>3, current_gene.citations); }
                 | CITATIONS_T SEPARATOR NEWLINE  /* empty citation = a No-OP */
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_gene.comments = stringlist_cons($<string>3, current_gene.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

dblinks          : DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id ignored_ids CLOSEPAREN NEWLINE
                 | DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id  CLOSEPAREN NEWLINE
                 ;

dblinks_id       : id
                   { current_gene.dblinks_ids = stringlist_cons($<string>1, current_gene.dblinks_ids); }
                 ;

dblinks_db       : id
                   { current_gene.dblinks_dbs = stringlist_cons($<string>1, current_gene.dblinks_dbs); }
                 ;

evidence         : EVIDENCE_T SEPARATOR strval NEWLINE
                   { current_gene.evidence = $<string>3; }
                 ;

interrupted      : INTERRUPTED_FLAG_T SEPARATOR strval NEWLINE
                   { current_gene.interrupted = $<string>3; }
                 ;

left_end_position : LEFT_END_POSITION_T SEPARATOR strval NEWLINE
                   { current_gene.left_end_position = $<string>3; }
                 ;

right_end_position : RIGHT_END_POSITION_T SEPARATOR strval NEWLINE
                   { current_gene.right_end_position = $<string>3; }
                 ;

product          : PRODUCT_T SEPARATOR id NEWLINE
                   { current_gene.products = stringlist_cons($<string>3, current_gene.products); }
                 ;

transcription_direction : TRANSCRIPTION_DIRECTION_T SEPARATOR strval NEWLINE
                   { current_gene.transcription_direction = $<string>3; }
                 ;

types            : TYPES_T SEPARATOR id NEWLINE
                   { current_gene.types = stringlist_cons($<string>3, current_gene.types); }
                 ;

component_of     : COMPONENT_OF_T SEPARATOR id NEWLINE
                   { current_gene.component_of = stringlist_cons($<string>3, current_gene.component_of); }
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR id NEWLINE
                   { current_gene.unique_id = $<string>3; }
                 ;

common_name      : COMMON_NAME_T SEPARATOR rhs NEWLINE
                   { current_gene.common_name = $<string>3; }
                 | COMMON_NAME_T SEPARATOR NEWLINE  /* empty common name = a No-OP */
                 ;

synonyms         : SYNONYMS_T SEPARATOR rhs NEWLINE
                   { current_gene.synonyms = stringlist_cons($<string>3, current_gene.synonyms); }
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
  printf("%s at `%s', line %d.\n", s, genetext, gene_lineno);
  gene_parse_error();
}

void
free_gene_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
gene_clear_entry(short freeflag) {
  /* Deallocate current_gene and clear it in preparation for next gene
   */
  current_gene.load_error = 0;
  
  free_gene_value(&current_gene.unique_id, freeflag);
  free_gene_value(&current_gene.common_name, freeflag);
  free_gene_value(&current_gene.evidence, freeflag);
  free_gene_value(&current_gene.left_end_position, freeflag);
  free_gene_value(&current_gene.right_end_position, freeflag);
  free_gene_value(&current_gene.interrupted, freeflag);
  free_gene_value(&current_gene.transcription_direction, freeflag);

  /* Free all lists of items */
  
  if (freeflag && current_gene.synonyms)
    free_stringlist(current_gene.synonyms);
  current_gene.synonyms = NULL;
  
  if (freeflag && current_gene.products)
    free_stringlist(current_gene.products);
  current_gene.products = NULL;
  
  if (freeflag && current_gene.types)
    free_stringlist(current_gene.types);
  current_gene.types = NULL;
  
  if (freeflag && current_gene.component_of)
    free_stringlist(current_gene.component_of);
  current_gene.component_of = NULL;
  
  if (freeflag && current_gene.comments) {
    free_stringlist(current_gene.comments);
  }
  current_gene.comments = NULL;

  if (freeflag && current_gene.citations) {
    free_stringlist(current_gene.citations);
  }
  current_gene.citations = NULL;

  if (freeflag && current_gene.dblinks_dbs) {
    free_stringlist(current_gene.dblinks_dbs);
  }
  current_gene.dblinks_dbs = NULL;  

  if (freeflag && current_gene.dblinks_ids) {
    free_stringlist(current_gene.dblinks_ids);
  }
  current_gene.dblinks_ids = NULL;  
}
						   
void
gene_next_entry(void) {

  /* Send this one to the database */
  gene_load_entry(&current_gene); 

  /* Then clear the structure for reuse */
  gene_clear_entry(1); 
  gene_records++;
}  

void
gene_parse_error(void) {
  current_gene.load_error = 1;
  gene_errors++;
}
