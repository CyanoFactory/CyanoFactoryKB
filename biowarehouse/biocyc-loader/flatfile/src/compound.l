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
 * flex scanner for compound.dat file of the Pathway/Genome Databases
 * See http://bioinformatics.ai.sri.com/ptools/flatfile-format.html
 */

%{
 #include <string.h>
 #include "string-util.h"
  
 #define DEBUG 0

 int compound_lineno = 1;
 extern char *lexedword;
 extern char lexedchar;

 /* BUMP increments the line counter for the source file */
 #define BUMP { compound_lineno++; }
 
 /* SAVEUNQUOTED saves the text except for the first and last chars. */
 #define SAVEUNQUOTED { lexedword = (char *) string_truncate(compoundtext+1, strlen(compoundtext)-2); if (DEBUG) printf("UnquotedSave %s.\n", lexedword); }
  
 /* SAVE copies scanned text into a global */
 #define SAVE { lexedword = (char *) strdup(compoundtext) ; if (DEBUG) printf("Saved %s.\n", lexedword); }
 
 /* RESET sets scanned text to empty, to handle some special cases of empty tokens */
 #define RESET { lexedword = (char *) strdup("") ; if (DEBUG) printf("Reset.\n"); }
 
 /* CONTINUE appends scanned text onto previous text, chopping off
    the lead "\n/" and replacing with a blank */
 #define CONTINUE { lexedword = string_concat(lexedword, strdup(compoundtext+2)); if (DEBUG) printf("CONTINUE %s.\n", lexedword);}
  
 /* APPEND appends scanned text onto previous text, addinf a blank in between */
 #define APPEND { lexedword = string_concat(lexedword, strdup(compoundtext)); if (DEBUG) printf("APPEND %s.\n", lexedword);}
%}

/*** Definitions ***/

whitespace  [[:blank:]]
newline     [\n]
openparen   [(]
closeparen  [)]
doublequote ["]
integer     ("+"|"-")?[[:digit:]]+
identifier  [[:alnum:][:punct:]]+
quotedid    [|][^|\n]*[|]
noparenid   [[:alnum:]'\-\+/,i\._\[\]{}<>=\*@$]+
contin      [\n][/]

separator     " - " 

linecomment ^[#].*[\n]
text        .+


%option stack noyywrap

%x fullRHS grabRHS continueRHS formulaRHS simpleRHS listRHS quotedString


%% /* Rules */

<INITIAL>{ /* State for start of a text line */
^"//\n"                  if (DEBUG) printf("INITIAL 1\n"); BUMP; return ENDOFRECORD;
{newline}                if (DEBUG) printf("INITIAL 2\n"); BUMP; return NEWLINE;
{whitespace}             ;  /* skip whitespace in INITIAL state */
{linecomment}            BUMP;

^"ATOM-CHARGES"          BEGIN(fullRHS); return ATOM_CHARGES_T;
^"CAS-REGISTRY-NUMBERS"  return CAS_REGISTRY_NUMBERS_T;
^"CHARGE"                return CHARGE_T;
^"CHEMICAL-FORMULA"      BEGIN(formulaRHS); return CHEMICAL_FORMULA_T;
^"CITATIONS"             BEGIN(fullRHS); return CITATIONS_T;
^"COMMENT"               BEGIN(fullRHS); return COMMENT_T;
^"COMMON-NAME"           BEGIN(fullRHS); return COMMON_NAME_T;
^"DBLINKS"               BEGIN(listRHS); return DBLINKS_T;
^"MOLECULAR-WEIGHT"      return MOLECULAR_WEIGHT_T;
^"PKA1"                  return PKA1_T;
^"PKA2"                  return PKA2_T;
^"PKA3"                  return PKA3_T;
^"SMILES"                return SMILES_T;
^"SYNONYMS"              BEGIN(fullRHS); return SYNONYMS_T;
^"SYSTEMATIC-NAME"       BEGIN(fullRHS); return SYSTEMATIC_NAME_T;
^"UNIQUE-ID"             return UNIQUE_ID_T;

^"ABBREV-NAME" |
^"ACTIVATORS-ALLOSTERIC-OF" |
^"ACTIVATORS-NONALLOSTERIC-OF" |
^"ACTIVATORS-UNKMECH-OF" |
^"COFACTORS-OF" |
^"COFACTORS-OR-PROSTHETIC-GROUPS-OF" |
^"COMPONENT-OF" |
^"HAS-NO-STRUCTURE?" |
^"INCHI" |
^"INHIBITORS-ALLOSTERIC-OF" |
^"INHIBITORS-COMPETITIVE-OF" |
^"INHIBITORS-IRREVERSIBLE-OF" |
^"INHIBITORS-NONCOMPETITIVE-OF" |
^"INHIBITORS-OTHER-OF" |
^"INHIBITORS-UNCOMPETITIVE-OF" |
^"INHIBITORS-UNKMECH-OF" |
^"PROSTHETIC-GROUPS-OF" |
^"MODIFIED-FORM" |
^"UNMODIFIED-FORM" |
^"N-NAME" |
^"N-1-NAME" |
^"N+1-NAME" |
^"ORIGINAL-FRAME-TYPE" |
^"REGULATES" |
^"SUPERATOMS" |
^"TYPES" |
^"ACID-BASE-VARIANTS"    BEGIN(fullRHS); return IGNORED_KEYWORD_T;

{separator}              BEGIN(simpleRHS); return SEPARATOR;
^{identifier}            if (DEBUG) printf("INITIAL ID\n"); SAVE; BEGIN(fullRHS); return WORD; /* Unrecognized attribute */
}

<simpleRHS>{
^"//\n"                  BUMP; return ENDOFRECORD;
{newline}                BUMP; return NEWLINE;
{whitespace}             ;  /* skip whitespace in simpleRHS state */
{identifier}             SAVE; BEGIN(INITIAL); return WORD; 
}

<formulaRHS>{
{whitespace}             ;  /* skip whitespace */
{separator}              return SEPARATOR;
[[:alpha:]][[:alnum:]]*  SAVE; return WORD;
{quotedid}               SAVE; return WORD;
{integer}                SAVE; return WORD;
{openparen}              return OPENPAREN;
{closeparen}             return CLOSEPAREN;
{newline}                BEGIN(INITIAL); BUMP; return NEWLINE;
{linecomment}            BEGIN(INITIAL); BUMP; return NEWLINE;
} 

<fullRHS>{
{separator}              BEGIN(grabRHS); return SEPARATOR;
} 

<grabRHS>{ /* rhs must generally be a RHS token but can be a NEWLINE */
{contin}/[^/]            if (DEBUG) printf("grabRHS 0\n"); RESET; BUMP; BEGIN(continueRHS);
{contin}/{newline}("//"|[^/]) 	if (DEBUG) printf("grabRHS 1\n"); RESET; BEGIN(INITIAL); BUMP;  return RHS; /* Hack for special case of empty continuation */
{text}/([\n]("//"|[^/])) if (DEBUG) printf("grabRHS 2\n"); SAVE;  return RHS;
{text}/{contin}[^/]      if (DEBUG) printf("grabRHS 3\n"); SAVE;  BEGIN(continueRHS);
{newline}                if (DEBUG) printf("grabRHS 4\n"); BEGIN(INITIAL); BUMP; return NEWLINE;
{linecomment}            if (DEBUG) printf("grabRHS 5\n"); BEGIN(INITIAL); BUMP; return NEWLINE;
} 

<continueRHS>{ /* Appends text from continuation line(s) onto existing text. */
{contin}/[^/]                 if (DEBUG) printf("continueRHS 0\n"); BUMP;
{text}                        if (DEBUG) printf("continueRHS 1\n"); APPEND;
{text}/{newline}("//"|[^/])   if (DEBUG) printf("continueRHS 2\n"); APPEND; return RHS;
{contin}/{newline}("//"|[^/]) if (DEBUG) printf("continueRHS 3\n"); BUMP;  return RHS;
{newline}/("//"|[^/])         if (DEBUG) printf("continueRHS 4\n"); BUMP;  BEGIN(INITIAL); return NEWLINE;
}

<listRHS>{
{whitespace}             ;  /* skip whitespace */
{separator}              if (DEBUG) printf("listRHS 0\n"); return SEPARATOR;
{quotedid}               if (DEBUG) printf("listRHS 1\n"); SAVEUNQUOTED; return WORD;
{noparenid}              if (DEBUG) printf("listRHS 2\n"); SAVE; return WORD;
{doublequote}            if (DEBUG) printf("listRHS 3\n"); RESET; BEGIN(quotedString);
{openparen}              if (DEBUG) printf("listRHS 4\n"); return OPENPAREN;
{closeparen}             if (DEBUG) printf("listRHS 5\n"); return CLOSEPAREN;
{newline}                if (DEBUG) printf("listRHS 6\n"); BEGIN(INITIAL); BUMP; return NEWLINE;
{linecomment}            if (DEBUG) printf("listRHS 7\n"); BEGIN(INITIAL); BUMP; return NEWLINE;
} 

<quotedString>{          /*  string beginning with double quote embedded in list */
[^"]*                    if (DEBUG) printf("quotedString 0\n"); SAVE;
{doublequote}            if (DEBUG) printf("quotedString 1\n"); BEGIN(listRHS); return WORD;
}

.                        SAVE; return WORD;          /* Error catch all */

%% /*** User Code ***/


