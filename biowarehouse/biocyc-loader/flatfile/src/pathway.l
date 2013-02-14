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
 * flex scanner for pathway.dat file of the Pathway/Genome Databases
 * See http://bioinformatics.ai.sri.com/ptools/flatfile-format.html
 */

%{
 #include <string.h>
 #include "string-util.h"

 #define DEBUG 0

 int pathway_lineno = 1;
 extern char *lexedword;
 extern char lexedchar;
 
 /* BUMP increments the line counter for the source file */
 #define BUMP { pathway_lineno++; }
 
 /* SAVEUNQUOTED saves the text except for the first and last chars. */
 #define SAVEUNQUOTED { lexedword = (char *) string_truncate(pathwaytext+1, strlen(pathwaytext)-2); if (DEBUG) printf("UnquotedSave %s.\n", lexedword); }
  
 /* SAVE copies scanned text into a global */
 #define SAVE { lexedword = (char *) strdup(pathwaytext); if (DEBUG) printf("Saved %s.\n", lexedword); }
 
 /* RESET sets scanned text to empty, to handle some special cases of empty tokens */
 #define RESET { lexedword = (char *) strdup(""); if (DEBUG) printf("Reset.\n"); }
 
 /* CONTINUE appends scanned text onto previous text, chopping off
    the lead "\n/" and replacing with a blank */
 #define CONTINUE { lexedword = string_concat(lexedword, strdup(pathwaytext+2)); if (DEBUG) printf("CONTINUE %s.\n", lexedword);}
  
 /* APPEND appends scanned text onto previous text, adding a blank in between */
 #define APPEND { lexedword = string_concat(lexedword, strdup(pathwaytext)); if (DEBUG) printf("APPEND %s.\n", lexedword);}
%}

/*** Definitions ***/

whitespace  [[:blank:]]
newline     [\n]
openparen   [(]
closeparen  [)]
element     [[:alpha:]]|[[:alpha:]][[:alpha:]]
integer     ("+"|"-")?[[:digit:]]+
identifier  [[:alnum:][:punct:]]+
noparenid   [[:alnum:]'\-\+/,i\._"\[\]{}<>=\*@$&]+
quotedid    [|][^|\n]*[|]
quotedstring ["][^"\n]*["]
contin      [\n][/]

separator   " - " 

linecomment ^[#].*[\n]
text        .+

%option stack noyywrap

%x fullRHS grabRHS continueRHS listRHS


%% /* Rules */

<INITIAL>{
^"//\n"                  BUMP; return ENDOFRECORD;
{newline}                BUMP; return NEWLINE;
{whitespace}             ;  /* skip whitespace in INITIAL state */
{linecomment}            BUMP;

^"CITATIONS"             BEGIN(fullRHS); return CITATIONS_T;
^"COMMENT"               BEGIN(fullRHS); return COMMENT_T;
^"COMMON-NAME"           BEGIN(fullRHS); return COMMON_NAME_T;
^"DBLINKS"               BEGIN(listRHS); return DBLINKS_T;
^"HYPOTHETICAL-REACTIONS"                return HYPOTHETICAL_REACTIONS_T;
^"NET-REACTION-EQUATION" BEGIN(fullRHS); return COMMENT_T;
^"PATHWAY-INTERACTIONS"  BEGIN(fullRHS); return COMMENT_T;
^"PATHWAY-LINKS"         BEGIN(listRHS); return PATHWAY_LINKS_T;
^"PREDECESSORS"          BEGIN(listRHS); return PREDECESSORS_T;
^"REACTION-LIST"         return REACTION_LIST_T;
^"SPECIES"               BEGIN(fullRHS); return COMMENT_T;
^"SUB-PATHWAYS"          return SUB_PATHWAYS_T;
^"SUPER-PATHWAYS"        return SUPER_PATHWAYS_T;
^"SYNONYMS"              BEGIN(fullRHS); return SYNONYMS_T;
^"UNIQUE-ID"             BEGIN(fullRHS); return UNIQUE_ID_T;

^"CLASS-INSTANCE-LINKS" |
^"CREDITS" |
^"DELTAGO" |
^"DELTAG0" |
^"ENZYME-USE" |
^"IN-PATHWAY" |
^"LAYOUT-ADVICE" |
^"PATHWAY-LINKS" |
^"POLYMERIZATION-LINKS" |
^"PRIMARIES" |
^"REACTION-LAYOUT" |
^"TAXONOMIC-RANGE" |
^"TYPES"                 BEGIN(fullRHS); return IGNORED_KEYWORD_T;

{separator}              return SEPARATOR;
^{identifier}            SAVE; BEGIN(fullRHS); return WORD; /* Unrecognized attribute */
{quotedid}               SAVEUNQUOTED; return WORD;
{identifier}             SAVE; return WORD;
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
{contin}/{newline}("//"|[^/]) if (DEBUG) printf("continueRHS 3\n"); BUMP; return RHS;
{newline}/("//"|[^/])         if (DEBUG) printf("continueRHS 4\n");  BEGIN(INITIAL); BUMP; return NEWLINE;
}

<listRHS>{
{whitespace}             ;  /* skip whitespace */
{separator}              if (DEBUG) printf("listRHS 1\n"); return SEPARATOR;
"."                      if (DEBUG) printf("listRHS 2\n"); return DOT;
{quotedid}               if (DEBUG) printf("listRHS 3\n"); SAVEUNQUOTED; return WORD;
{noparenid}              if (DEBUG) printf("listRHS 4\n"); SAVE; return WORD;
{openparen}              if (DEBUG) printf("listRHS 5\n"); return OPENPAREN;
{closeparen}             if (DEBUG) printf("listRHS 6\n"); return CLOSEPAREN;
{newline}                if (DEBUG) printf("listRHS 7\n"); BEGIN(INITIAL); BUMP; return NEWLINE;
{linecomment}            if (DEBUG) printf("listRHS 8\n"); BEGIN(INITIAL); BUMP; return NEWLINE;
":"{noparenid}           if (DEBUG) printf("listRHS 9\n"); SAVE; return WORD;
{quotedstring}           if (DEBUG) printf("listRHS 10\n"); SAVE; return WORD;
}

.                        SAVE; return WORD;          /* Error catch all */

%% /*** User Code ***/
