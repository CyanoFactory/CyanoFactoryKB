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
 * Bison Parser for PROTSEQ
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.protseq.c"
  #include "string-util.h"
  #include "protseq-parse.h"
  
  #define YYDEBUG 1
  #define DEBUG 0
   
  int protseq_records = 0;
  int protseq_errors  = 0;
  short first_protseq = 1; /* flag for first entry in file */
  struct protseq_entry current_protseq;

 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token NEWLINE WORD RHS START

%%

protseqfile     : protseqrecords
                 | // no entries
                ;

protseqrecords  : protseqrecords protseqrecord
                | protseqrecord
                ;

protseqrecord   : proteininfo sequence 
                  { protseq_next_entry(); first_protseq = 0; }
                ;

proteininfo     : START id rhs NEWLINE
                  { current_protseq.protein = $<string>2; }
                ;

sequence        : fragments NEWLINE
                  { current_protseq.aa_sequence = string_append(current_protseq.aa_sequence, $<string>1 ); }
                | /* empty */
                ;

fragments       : fragments fragment
                   { $<string>$ = string_concat($<string>1, $<string>2); }
                | fragment
                ;

fragment        : strval
		;


/* Tokens */

strval           : WORD
                   { $<string>$ = lexedword; }
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
  printf("%s at `%s', line %d.\n", s, protseqtext, protseq_lineno);
  protseq_parse_error();
}

void
free_protseq_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
protseq_clear_entry(short freeflag) {
  /* Deallocate current_protseq and clear it in preparation for next protseq
   */
  current_protseq.load_error = 0;
  
  free_protseq_value(&current_protseq.protein, freeflag);
  free_protseq_value(&current_protseq.aa_sequence, freeflag);
}
						   
void
protseq_next_entry(void) {
  yydebug = DEBUG;

  /* Send this one to the database */
  protseq_load_entry(&current_protseq); 

  /* Then clear the structure for reuse */
  protseq_clear_entry(1); 
  protseq_records++;
}  

void
protseq_parse_error(void) {
  current_protseq.load_error = 1;
  protseq_errors++;
}
