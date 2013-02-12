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
 * Bison Parser for COMPOUND
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.compound.c"
  #include "string-util.h"
  #include "compound-parse.h"
  
  void skip_compound_line(void);
  int compound_records = 0;
  int compound_errors  = 0;
  int last_compound_line_skipped = -1;
  short first_compound = 1;
  struct compound_entry current_compound;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN
%token ATOM_CHARGES_T CAS_REGISTRY_NUMBERS_T CHARGE_T CHEMICAL_FORMULA_T CITATIONS_T
%token COMMON_NAME_T COMMENT_T MOLECULAR_WEIGHT_T PKA1_T PKA2_T PKA3_T
%token SMILES_T SYNONYMS_T SYSTEMATIC_NAME_T UNIQUE_ID_T DBLINKS_T IGNORED_KEYWORD_T

%%

compoundfile     : compoundrecords
                 | // no entries
                 ;

compoundrecords  : compoundrecords compoundrecord
                 | compoundrecord
                 ;

compoundrecord   : compoundslots ENDOFRECORD
                   { compound_next_entry(); first_compound = 0; }
                 ;

compoundslots    : compoundslot
                 | compoundslots compoundslot
                 ;

compoundslot :     unique_id
		 | cas_registry_numbers
                 | charge
                 | chemical_formula
                 | citation
                 | comment
                 | common_name
                 | dblinks
                 | molecular_weight
                 | pka1
                 | pka2
                 | pka3
                 | systematic_name
                 | smiles
		 | synonyms
                 | atom_charges
                 | ignored_keyword
                 | unknown_keyword
                 | error
                 ;

cas_registry_numbers : CAS_REGISTRY_NUMBERS_T SEPARATOR strval NEWLINE
                   { current_compound.cas_registry_numbers = $<string>3; }
                 ;

charge           : CHARGE_T SEPARATOR strval NEWLINE
                   { current_compound.charge = $<string>3; }
                 ;

chemical_formula : CHEMICAL_FORMULA_T SEPARATOR OPENPAREN id number CLOSEPAREN NEWLINE
                   { add_element_to_formula($<string>4, $<string>5); }
                 ;

citation         : CITATIONS_T SEPARATOR rhs NEWLINE
                   { current_compound.citations = stringlist_cons($<string>3, current_compound.citations); }
                 | CITATIONS_T SEPARATOR NEWLINE  /* empty citation = a No-OP */
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_compound.comments = stringlist_cons($<string>3, current_compound.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

common_name      : COMMON_NAME_T SEPARATOR rhs NEWLINE
                   { current_compound.common_name = $<string>3; }
                 ;

dblinks          : DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id ignored_ids CLOSEPAREN NEWLINE
                 | DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id CLOSEPAREN NEWLINE
                 ;

dblinks_id       : id
                   { current_compound.dblinks_ids = stringlist_cons($<string>1, current_compound.dblinks_ids); }
                 ;

dblinks_db       : id
                   { current_compound.dblinks_dbs = stringlist_cons($<string>1, current_compound.dblinks_dbs); }
                 ;

smiles           : SMILES_T SEPARATOR strval NEWLINE
                   { current_compound.smiles = $<string>3; }
                 ;

molecular_weight : MOLECULAR_WEIGHT_T SEPARATOR strval NEWLINE
                   { current_compound.molecular_weight = $<string>3; }
                 ;
pka1             : PKA1_T SEPARATOR strval NEWLINE
                   { current_compound.pka1 = $<string>3; }
                 ;

pka2             : PKA2_T SEPARATOR strval NEWLINE
                   { current_compound.pka2 = $<string>3; }
                 ;

pka3             : PKA3_T SEPARATOR strval NEWLINE
                   { current_compound.pka3 = $<string>3; }
                 ;

systematic_name  : SYSTEMATIC_NAME_T SEPARATOR rhs NEWLINE
                  { current_compound.systematic_name = $<string>3;  }
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR id NEWLINE
                   { current_compound.unique_id = $<string>3; }
                 ;

synonyms        : SYNONYMS_T SEPARATOR rhs NEWLINE
                   { current_compound.synonyms = stringlist_cons($<string>3, current_compound.synonyms); }
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
atom_charges     : ATOM_CHARGES_T SEPARATOR rhs NEWLINE
                 ;

/* Tokens */

strval           : WORD
                   { $<string>$ = lexedword;  }
                 ;

number           : WORD
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
  printf("%s at `%s', line %d.\n", s, compoundtext, compound_lineno);
  compound_parse_error();
}

void
free_compound_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
compound_clear_entry(short freeflag) {
  /* Deallocate current_compound and clear it in preparation for next compound
   */
  current_compound.load_error = 0;

  // Kludge to disable silent numeric conversion errors in number_column
  current_compound.molecular_weight_ind = INDICATE_NULL;
  
  free_compound_value(&current_compound.unique_id, freeflag);
  free_compound_value(&current_compound.common_name, freeflag);
  free_compound_value(&current_compound.systematic_name, freeflag);
  free_compound_value(&current_compound.cas_registry_numbers, freeflag);
  free_compound_value(&current_compound.smiles, freeflag);
  free_compound_value(&current_compound.formula, freeflag);
  free_compound_value(&current_compound.charge, freeflag);
  free_compound_value(&current_compound.molecular_weight , freeflag);
  free_compound_value(&current_compound.pka1, freeflag);
  free_compound_value(&current_compound.pka2, freeflag);
  free_compound_value(&current_compound.pka3, freeflag);
  
  if (freeflag && current_compound.synonyms) {
    free_stringlist(current_compound.synonyms);
  }
  current_compound.synonyms = NULL;
  
  if (freeflag && current_compound.comments) {
    free_stringlist(current_compound.comments);
  }
  current_compound.comments = NULL;

  if (freeflag && current_compound.citations) {
    free_stringlist(current_compound.citations);
  }
  current_compound.citations = NULL;

  if (freeflag && current_compound.dblinks_dbs) {
    free_stringlist(current_compound.dblinks_dbs);
  }
  current_compound.dblinks_dbs = NULL;  

  if (freeflag && current_compound.dblinks_ids) {
    free_stringlist(current_compound.dblinks_ids);
  }
  current_compound.dblinks_ids = NULL;  
}
						   
void
compound_next_entry(void) {

  /* Send this one to the database */
  compound_load_entry(&current_compound); 

  /* Then clear the structure for reuse */
  compound_clear_entry(1); 
  compound_records++;
}  

void
compound_parse_error(void) {
  current_compound.load_error = 1;
  compound_errors++;
}

void
add_element_to_formula(char *element, char *number) {
  /* Builds up the formula by concatenating element and number from an
     CHEMICAL-FORMULA - (element number) input
  */
  if (!current_compound.formula)
    current_compound.formula = element;
  else
    current_compound.formula = string_append(current_compound.formula, element);
  current_compound.formula = string_append(current_compound.formula, number);
}

void
skip_compound_line(void) {
  if (last_compound_line_skipped != compound_records) {
    last_compound_line_skipped = compound_records; /* remember last line that was skipped due to error */
    compound_parse_error();
  }
}
