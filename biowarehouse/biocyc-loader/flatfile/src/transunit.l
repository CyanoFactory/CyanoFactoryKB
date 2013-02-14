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
 * flex scanner for transunits.dat file of the Pathway/Genome Databases
 * See http://bioinformatics.ai.sri.com/ptools/flatfile-format.html
 */

%{
 #include <string.h>
 #include "string-util.h"
  
 #define DEBUG 0

 int transunit_lineno = 1;
 extern char *lexedword;
 extern char lexedchar;
 
 /* BUMP increments the line counter for the source file */
 #define BUMP { transunit_lineno++; }
 
 /* SAVEUNQUOTED saves the text except for the first and last chars. */
 #define SAVEUNQUOTED { lexedword = (char *) string_truncate(transunittext+1, strlen(transunittext)-2); if (DEBUG) printf("UnquotedSave %s.\n", lexedword); }
  
 /* SAVE copies scanned text into a global */
 #define SAVE { lexedword = (char *) strdup(transunittext) ; if (DEBUG) printf("Saved %s.\n", lexedword); }
 
 /* RESET sets scanned text to empty, to handle some special cases of empty tokens */
 #define RESET { lexedword = (char *) strdup("") ; if (DEBUG) printf("Reset.\n"); }
 
 /* CONTINUE appends scanned text onto previous text, chopping off
    the lead "\n/" and replacing with a blank */
 #define CONTINUE { lexedword = string_concat(lexedword, strdup(transunittext+2)); if (DEBUG) printf("CONTINUE %s.\n", lexedword);}
  
 /* APPEND appends scanned text onto previous text, addinf a blank in between */
 #define APPEND { lexedword = string_concat(lexedword, strdup(transunittext)); if (DEBUG) printf("APPEND %s.\n", lexedword);}
%}

/*** Definitions ***/

whitespace  [[:blank:]]
newline     [\n]
openparen   [(]
closeparen  [)]
doublequote ["]
element     [[:alpha:]]|[[:alpha:]][[:alpha:]]
integer     ("+"|"-")?[[:digit:]]+
identifier  [[:alnum:][:punct:]]+
quotedid    [|][^|\n]*[|]
noparenid   [[:alnum:]'\-\+/,i\._\[\]{}<>=\*@$]+
contin      [\n][/]

separator   " - " 

linecomment ^[#].*
text        .+

%option stack noyywrap

%x fullRHS grabRHS continueRHS listRHS quotedString


%% /* Rules */

<INITIAL>{
^"//\n"                  BUMP; return ENDOFRECORD;
{newline}                BUMP; return NEWLINE;
{whitespace}             ;  /* skip whitespace in INITIAL state */
{linecomment}            ;
{linecomment}[\n]        BUMP;

^"CITATIONS"             BEGIN(fullRHS); return CITATIONS_T;
^"COMMENT"               BEGIN(fullRHS); return COMMENT_T;
^"COMMON-NAME"           BEGIN(fullRHS); return COMMON_NAME_T; 
^"COMPONENTS"            return COMPONENTS_T;
^"DBLINKS"               BEGIN(listRHS); return DBLINKS_T; 
^"SYNONYMS"              BEGIN(fullRHS); return SYNONYMS_T;
^"UNIQUE-ID"             BEGIN(fullRHS); return UNIQUE_ID_T;

^"EXTENT-UNKNOWN?",
^"TYPES"                 BEGIN(fullRHS); return IGNORED_KEYWORD_T;

{separator}              return SEPARATOR;
^{identifier}            SAVE; BEGIN(fullRHS); return WORD; /* Unrecognized attribute */
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
{linecomment}[\n]        if (DEBUG) printf("grabRHS 5\n"); BEGIN(INITIAL); BUMP; return NEWLINE;
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
{linecomment}[\n]        if (DEBUG) printf("listRHS 7\n"); BEGIN(INITIAL); BUMP; return NEWLINE;
} 

<quotedString>{ /*  string beginning with " embedded in list */
[^"]*                    if (DEBUG) printf("quotedString 0\n"); SAVE;
{doublequote}            if (DEBUG) printf("quotedString 1\n"); BEGIN(listRHS); return WORD;
}

.                        SAVE; return WORD;          /* Error catch all */

%% /*** User Code ***/
