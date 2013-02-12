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
 * Bison Parser for PATHWAY
 */

%{
  #include <stdio.h>
  #include "main.h"
  #include "lex.pathway.c"
  #include "string-util.h"
  #include "widtable.h"
  #include "pathway-parse.h"
  
  int pathway_records = 0;
  int pathway_errors  = 0;
  short first_pathway = 1; /* flag for first entry in file */
  struct pathway_entry current_pathway;

  void free_pathway_link_list(struct pathway_link *list);
  struct pathway_link * pathway_link_cons(char *chemical, struct stringlist *pathways, struct pathway_link *list);
 %}

%union {
  char              *string;
  struct stringlist *list;
};

%token ENDOFRECORD NEWLINE WORD RHS
%token SEPARATOR OPENPAREN CLOSEPAREN DOT
%token CITATIONS_T COMMENT_T COMMON_NAME_T DBLINKS_T SYNONYMS_T UNIQUE_ID_T IGNORED_KEYWORD_T
%token HYPOTHETICAL_REACTIONS_T PREDECESSORS_T SUB_PATHWAYS_T SUPER_PATHWAYS_T PATHWAY_LINKS_T REACTION_LIST_T

%%

pathwayfile      : pathwayrecords
                 | // no entries
                 ;

pathwayrecords  : pathwayrecords pathwayrecord
                 | pathwayrecord
                 ;

pathwayrecord   : pathwayslots ENDOFRECORD
                   { pathway_next_entry(); first_pathway = 0; }
                 ;

pathwayslots    : pathwayslot
                 | pathwayslots pathwayslot
                 ;

pathwayslot     : unique_id
                 | citation
                 | comment
                 | common_name
                 | dblinks
		 | synonyms
                 | pathway_links
                 | reaction_list
                 | predecessor
                 | super_pathway
                 | sub_pathway
                 | hypothetical_reaction
		 | ignored_keyword
		 | unknown_keyword
                 | error
                 ;

unique_id        : UNIQUE_ID_T SEPARATOR rhs NEWLINE
                   { current_pathway.unique_id = $<string>3; }
                 ;

citation         : CITATIONS_T SEPARATOR rhs NEWLINE
                   { current_pathway.citations = stringlist_cons($<string>3, current_pathway.citations); }
                 | CITATIONS_T SEPARATOR NEWLINE  /* empty citation = a No-OP */
                 ;

comment          : COMMENT_T SEPARATOR rhs NEWLINE
                   { current_pathway.comments = stringlist_cons($<string>3, current_pathway.comments); }
                 | COMMENT_T SEPARATOR NEWLINE  /* empty comment = a No-OP */
                 ;

common_name      : COMMON_NAME_T SEPARATOR rhs NEWLINE
                   { current_pathway.common_name = $<string>3; }
                 ;

dblinks          : DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id ignored_ids CLOSEPAREN NEWLINE
                 | DBLINKS_T SEPARATOR OPENPAREN dblinks_db dblinks_id  CLOSEPAREN NEWLINE
                 ;

dblinks_id       : id
                   { current_pathway.dblinks_ids = stringlist_cons($<string>1, current_pathway.dblinks_ids); }
                 ;

dblinks_db       : id
                   { current_pathway.dblinks_dbs = stringlist_cons($<string>1, current_pathway.dblinks_dbs); }
                 ;

sub_pathway      : SUB_PATHWAYS_T SEPARATOR id NEWLINE
                   { current_pathway.sub_pathways = stringlist_cons($<string>3, current_pathway.sub_pathways); } 
                 ;

super_pathway    : SUPER_PATHWAYS_T SEPARATOR id NEWLINE
                   { current_pathway.super_pathways = stringlist_cons($<string>3, current_pathway.super_pathways); } 
                 ;

reaction_list    : REACTION_LIST_T SEPARATOR id NEWLINE
                   { current_pathway.reactions = stringlist_cons($<string>3, current_pathway.reactions); } 
                 ;

pathway_links    : PATHWAY_LINKS_T SEPARATOR OPENPAREN id pathway_ids CLOSEPAREN NEWLINE
                   { current_pathway.links = pathway_link_cons($<string>4, current_pathway.templist, current_pathway.links);
		     current_pathway.templist = NULL; }
                 |  PATHWAY_LINKS_T SEPARATOR OPENPAREN id CLOSEPAREN NEWLINE  // ignore these, probably a lone pathway ID
                 ;

predecessor      : PREDECESSORS_T SEPARATOR OPENPAREN id CLOSEPAREN NEWLINE  // no predecessors for this successor
                   { current_pathway.successors   = stringlist_cons($<string>4, current_pathway.successors);
                     current_pathway.predecessors = stringlist_cons(strdup(" "), current_pathway.predecessors); }
                 | PREDECESSORS_T SEPARATOR OPENPAREN id predecessor_ids CLOSEPAREN NEWLINE  // 1+ predecessors
                   {current_pathway.successors   = stringlist_rep($<string>4, stringlist_length(current_pathway.predecessors)); // N copies of successor
		    current_pathway.predecessors = stringlist_cons($<string>5, current_pathway.predecessors);
                    current_pathway.temppreds = NULL; }  // clear list of preds for this one successor
/*                  | PREDECESSORS_T SEPARATOR OPENPAREN id predecessor_ids CLOSEPAREN NEWLINE  // several predecessors */
/*                    {current_pathway.successors    = stringlist_cons($<string>4, current_pathway.successors); */
/* 		    current_pathway.predecessors  = stringlist_cons($<string>5, current_pathway.predecessors); } */
/*                  | PREDECESSORS_T SEPARATOR OPENPAREN id id CLOSEPAREN NEWLINE  /\* one predecessor *\/ */
/*                    {current_pathway.successors    = stringlist_cons($<string>4, current_pathway.successors); */
/* 		    current_pathway.predecessors  = stringlist_cons($<string>5, current_pathway.predecessors); } */
/* 		 | PREDECESSORS_T SEPARATOR OPENPAREN id id id CLOSEPAREN NEWLINE  /\* two predecessors *\/ */
/*                    {current_pathway.successors    = stringlist_cons($<string>4, current_pathway.successors); */
/* 		    current_pathway.predecessors  = stringlist_cons($<string>5, current_pathway.predecessors); */
/* 		    current_pathway.successors    = stringlist_cons(strdup($<string>4), current_pathway.successors); */
/* 		    current_pathway.predecessors  = stringlist_cons($<string>6, current_pathway.predecessors); } */
                 | PREDECESSORS_T SEPARATOR id NEWLINE  // non-list, subpathway predecessor
                   // ignore, info is in SUBAPTHWAYS attr
                 ;

hypothetical_reaction : HYPOTHETICAL_REACTIONS_T SEPARATOR id NEWLINE
                        { current_pathway.hypothetical_rxns = stringlist_cons($<string>3, current_pathway.hypothetical_rxns); }
                 ;

synonyms         : SYNONYMS_T SEPARATOR rhs NEWLINE
                   { current_pathway.synonyms = stringlist_cons($<string>3, current_pathway.synonyms); }
                 | SYNONYMS_T SEPARATOR NEWLINE  /* empty synonyms = a No-OP */
                 ;

predecessor_ids  : predecessor_ids predecessor_id
                 | predecessor_id
//------- looks wrong { current_pathway.temppreds = stringlist_cons($<string>1, current_pathway.temppreds); }
                 ;

predecessor_id   : id
                   { current_pathway.temppreds = stringlist_cons($<string>1, current_pathway.temppreds); }
                 ;

pathway_ids      : pathway_ids pathway_id
                 | pathway_id
//------- looks wrong { current_pathway.templist = stringlist_cons($<string>1, current_pathway.templist); }
                 ;

pathway_id       : id
                   { current_pathway.templist = stringlist_cons($<string>1, current_pathway.templist); }
                 | OPENPAREN id DOT id CLOSEPAREN   // 2nd ID indicates direction
                   { current_pathway.templist = stringlist_cons($<string>2, current_pathway.templist); }
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
  printf("%s at `%s', line %d.\n", s, pathwaytext, pathway_lineno);
  pathway_parse_error();
}

void
free_pathway_value(char **value, short freeflag) {
  if (freeflag) free(*value);
  *value = NULL;
}

void
pathway_clear_entry(short freeflag) {
  /* Deallocate current_pathway and clear it in preparation for next pathway
   */
  current_pathway.load_error = 0;
  
  free_pathway_value(&current_pathway.unique_id, freeflag);
  free_pathway_value(&current_pathway.common_name, freeflag);
  
  if (freeflag && current_pathway.sub_pathways)
    free_stringlist(current_pathway.sub_pathways);
  current_pathway.sub_pathways = NULL;
  
  if (freeflag && current_pathway.successors)
    free_stringlist(current_pathway.successors);
  current_pathway.successors = NULL;
  if (freeflag && current_pathway.predecessors)
    free_stringlist(current_pathway.predecessors);
  current_pathway.predecessors = NULL;
  
  if (freeflag && current_pathway.links)
    free_pathway_link_list(current_pathway.links);
  current_pathway.links = NULL;

  current_pathway.templist = NULL;  /* Don't free this list - it's copied elsewhere and freed separately */
  current_pathway.temppreds = NULL;  /* Don't free this list - it's copied elsewhere and freed separately */
  
  if (freeflag && current_pathway.super_pathways)
    free_stringlist(current_pathway.super_pathways);
  current_pathway.super_pathways = NULL;
  
  if (freeflag && current_pathway.reactions)
    free_stringlist(current_pathway.reactions);
  current_pathway.reactions = NULL;
  
  if (freeflag && current_pathway.hypothetical_rxns)
    free_stringlist(current_pathway.hypothetical_rxns);
  current_pathway.hypothetical_rxns = NULL;
  
  if (freeflag && current_pathway.synonyms)
    free_stringlist(current_pathway.synonyms);
  current_pathway.synonyms = NULL;
  
  if (freeflag && current_pathway.citations) {
    free_stringlist(current_pathway.citations);
  }
  current_pathway.citations = NULL;

  if (freeflag && current_pathway.comments) {
    free_stringlist(current_pathway.comments);
  }
  current_pathway.comments = NULL;

  if (freeflag && current_pathway.dblinks_dbs) {
    free_stringlist(current_pathway.dblinks_dbs);
  }
  current_pathway.dblinks_dbs = NULL;  

  if (freeflag && current_pathway.dblinks_ids) {
    free_stringlist(current_pathway.dblinks_ids);
  }
  current_pathway.dblinks_ids = NULL;  
}
						   
void
pathway_next_entry(void) {

  /* Send this one to the database */
  pathway_load_entry(&current_pathway); 

  /* Then clear the structure for reuse */
  pathway_clear_entry(1); 
  pathway_records++;
}  

void
pathway_parse_error(void) {
  current_pathway.load_error = 1;
  pathway_errors++;
}


/*** pathway_link utilities ***/

void
free_pathway_link_list(struct pathway_link *list) {
  if (list) {
    free_pathway_link_list(list->next);
    free_stringlist(list->pathways);
    free(list->chemical);
    free(list);
  }
}

struct pathway_link *
pathway_link_cons(char *chemical, struct stringlist *pathways, struct pathway_link *list) {

  struct pathway_link *newlist;

  newlist = (struct pathway_link *) malloc (sizeof(struct pathway_link));
  if (!newlist)
    fatal_error("malloc failed in pathway_link_cons");

  newlist->next     = list;
  newlist->chemical = chemical;
  newlist->pathways = pathways;
  return(newlist);
}

