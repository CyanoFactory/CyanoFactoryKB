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
 * flex scanner for protseq.fasta file of the Pathway/Genome Databases
 * See http://bioinformatics.ai.sri.com/ptools/flatfile-format.html
 */

%{
 #include <string.h>
 #include "string-util.h"
 int protseq_lineno = 1;
 extern char *lexedword;
 extern char lexedchar;

 /* BUMP increments the line counter for the source file */
 #define BUMP { protseq_lineno++; }
 
 #define SAVE { lexedword = (char *) strdup(protseqtext); }
%}

/*** Definitions ***/

newline     [\n]
identifier  [[:alnum:][:punct:]]+

text        .*

whitespace  [[:blank:]]

aafragment  [[:alpha:]*]+

%option stack noyywrap

%x grabID grabRHS AASequence


%% /* Rules */

<INITIAL>{
^"> "                    BEGIN(grabID); return START; /* skip record starter */
^">"                     BEGIN(grabID); return START; /* skip record starter */
">"                      BEGIN(grabID); return START; /* skip record starter */
{newline}                BUMP; return NEWLINE;
}

<grabID>{
{whitespace}             ; // Ignore 
{identifier}             BEGIN(grabRHS); SAVE; return WORD;
{newline}                BUMP; return NEWLINE;
}

<grabRHS>{
{text}                   SAVE; return RHS;
{newline}/">"            BEGIN(INITIAL); BUMP; return NEWLINE;  /* Allow for no sequence lines */
{newline}                BEGIN(AASequence); BUMP; return NEWLINE;
} 

<AASequence>{
{aafragment}/"\n>"       BEGIN(INITIAL); SAVE; return WORD;
{aafragment}             SAVE; return WORD;
{newline}/.              BUMP; /* newline followed by non-EOF */
{newline}                BUMP; return NEWLINE; /* newline followed by EOF */
}

.                        SAVE; return WORD;          /* Error catch all */

%% /*** User Code ***/
