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
 * Bison Parser for PROTEIN
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.protein.c"
  #include "string-util.h"
  #include "protein-parse.h"
  
  int protein_records = 0;
  int protein_errors  = 0;
  short first_protein = 1; /* flag for first entry in file */
  struct protein_entry current_protein;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN
%token COMPONENTS_T COEFFICIENT_T
%token CITATIONS_T COMMENT_T COMMON_NAME_T DBLINKS_T SYNONYMS_T UNIQUE_ID_T IGNORED_KEYWORD_T GENE_T MOLECULAR_WEIGHT_EXP_T MOLECULAR_WEIGHT_T PI_T LOCATIONS_T GO_TERMS_T

%%

proteinfile     : proteinrecords
                 | // no entries
                 ;

proteinrecords  : proteinrecords proteinrecord
                 | proteinrecord
                 ;

proteinrecord   : proteinslots ENDOFRECORD
                   { protein_next_entry(); first_protein = 0; }
                 ;

proteinslots    : proteinslot
                 | proteinslots proteinslot
                 ;

proteinslot     : unique_id
                 | citation
                 | coefficient
                 | comment
                 | common_name
                 | component
                 | dblinks
                 | gene
                 | go_terms
                 | locations
                 | molecular_weight
                 | molecular_weight_exp
                 | pi
		 | synonyms
		 | ignored_keyword
		 | unknown_keyword
                 | error
                 ;

citation         : CITATIONS_T SEPARATOR rhs NEWLINE
                   { current_protein.citations = stringlist_cons($<string>3, current_protein.citations); }
                 | CITATIONS_T SEPARATOR NEWLINE  /* empty citation = a No-OP */
                 ;

coefficient      : COEFFICIENT_T SEPARATOR strval NEWLINE
                   { free(current_protein.coeffs->string);
		     current_protein.coeffs->string = $<string>3; }  /* overwrite default */
                 ;

component        : COMPONENTS_T SEPARATOR id NEWLINE  /* defines a subunit */
                   { current_protein.components = stringlist_cons($<string>3, current_protein.components);
		     current_protein.coeffs = stringlist_cons(strdup("1"), current_protein.coeffs); } /* default */
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR rhs NEWLINE
                   { current_protein.unique_id = $<string>3; }
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_protein.comments = stringlist_cons($<string>3, current_protein.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

common_name      : COMMON_NAME_T SEPARATOR rhs NEWLINE
                   { current_protein.common_name = $<string>3;  }
                 | COMMON_NAME_T SEPARATOR NEWLINE  /* empty common name = a No-OP */
                 ;

dblinks          : DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id ignored_ids CLOSEPAREN NEWLINE
                 | DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id  CLOSEPAREN NEWLINE
                 ;

dblinks_id       : id
                   { current_protein.dblinks_ids = stringlist_cons($<string>1, current_protein.dblinks_ids); }
                 ;

dblinks_db       : id
                   { current_protein.dblinks_dbs = stringlist_cons($<string>1, current_protein.dblinks_dbs); }
                 ;

gene             : GENE_T SEPARATOR id NEWLINE
                   { current_protein.gene = $<string>3; }
                 ;

go_terms         : GO_TERMS_T SEPARATOR id NEWLINE
                 { current_protein.go_terms = stringlist_cons($<string>3, current_protein.go_terms); }
                 | GO_TERMS_T SEPARATOR NEWLINE /* empty go_terms = a No-OP */
                 ;

locations        : LOCATIONS_T SEPARATOR rhs NEWLINE
                 { current_protein.locations = stringlist_cons($<string>3, current_protein.locations); }
                 | LOCATIONS_T SEPARATOR NEWLINE /* empty locations = a No-OP */
                 ;

molecular_weight : MOLECULAR_WEIGHT_T SEPARATOR strval NEWLINE
                   { current_protein.molecular_weight = $<string>3; }
                 ;

molecular_weight_exp : MOLECULAR_WEIGHT_EXP_T SEPARATOR strval NEWLINE
                   { current_protein.molecular_weight_exp = $<string>3; }
                 ;

pi               : PI_T SEPARATOR strval NEWLINE
                   { current_protein.pi = $<string>3; }
                 ;

synonyms         : SYNONYMS_T SEPARATOR rhs NEWLINE
                   { current_protein.synonyms = stringlist_cons($<string>3, current_protein.synonyms); }
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
  printf("%s at `%s', line %d.\n", s, proteintext, protein_lineno);
  protein_parse_error();
}

void
free_protein_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
protein_clear_entry(short freeflag) {
  /* Deallocate current_protein and clear it in preparation for next protein
   */
  current_protein.load_error = 0;
  
  // Kludge to disable silent numeric conversion errors in number_column
  current_protein.pi_ind = INDICATE_NULL;
  current_protein.molecular_weight_ind = INDICATE_NULL;
  current_protein.molecular_weight_exp_ind = INDICATE_NULL;
  
  free_protein_value(&current_protein.unique_id, freeflag);
  free_protein_value(&current_protein.common_name, freeflag);
  free_protein_value(&current_protein.gene, freeflag);
  free_protein_value(&current_protein.molecular_weight, freeflag);
  free_protein_value(&current_protein.molecular_weight_exp, freeflag);
  free_protein_value(&current_protein.pi, freeflag);
  
  if (freeflag && current_protein.coeffs)
    free_stringlist(current_protein.coeffs);
  current_protein.coeffs = NULL;
  
  if (freeflag && current_protein.components)
    free_stringlist(current_protein.components);
  current_protein.components = NULL;
  
  if (freeflag && current_protein.synonyms)
    free_stringlist(current_protein.synonyms);
  current_protein.synonyms = NULL;
  
  if (freeflag && current_protein.locations)
    free_stringlist(current_protein.locations);
  current_protein.locations = NULL;
  
  if (freeflag && current_protein.go_terms)
    free_stringlist(current_protein.go_terms);
  current_protein.go_terms = NULL;
  
  if (freeflag && current_protein.citations) {
    free_stringlist(current_protein.citations);
  }
  current_protein.citations = NULL;
  
  if (freeflag && current_protein.comments) {
    free_stringlist(current_protein.comments);
  }
  current_protein.comments = NULL;

  if (freeflag && current_protein.dblinks_dbs) {
    free_stringlist(current_protein.dblinks_dbs);
  }
  current_protein.dblinks_dbs = NULL;  

  if (freeflag && current_protein.dblinks_ids) {
    free_stringlist(current_protein.dblinks_ids);
  }
  current_protein.dblinks_ids = NULL;  
}
						   
void
protein_next_entry(void) {

  /* Send this one to the database */
  protein_load_entry(&current_protein); 

  /* Then clear the structure for reuse */
  protein_clear_entry(1); 
  protein_records++;
}  

void
protein_parse_error(void) {
  current_protein.load_error = 1;
  protein_errors++;
}
