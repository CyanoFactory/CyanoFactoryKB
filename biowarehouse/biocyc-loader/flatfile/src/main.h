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
#ifndef _MAIN_H
#define _MAIN_H 1

// System includes
#include <stdio.h>
#include <stddef.h>
#include <assert.h>
#include <string.h>

// Control declaration of global vars and externs for them
#ifdef MAIN_PGM
#define EXTERN
#else
#define EXTERN extern
#endif

// Warehouse includes
#include "wh.h"

/*
 *  Definitions used throughout the loader.
 */

extern char *yytext; /* Added 3.6.1, apparently required by some platforms */

/*** Used by utilities in utils/src/c ***/

EXTERN int dataset_wid;    
EXTERN int multifun_dataset_wid;
EXTERN int go_dataset_wid;

EXTERN int commit_freq;  /* commit every this many WIDs, or sometimes input records;
			    nonzero also causes commit after load of each table */
EXTERN char *userid; /* for Oracle: USERID/PASSWORD
			for MYSQL : USERID */



EXTERN int organism_wid;

EXTERN WIDTABLE citation_wids;
EXTERN WIDTABLE protein_wids;
EXTERN WIDTABLE pathway_wids;
EXTERN WIDTABLE promoter_positions;
EXTERN WIDTABLE transunit_wids;

EXTERN WIDNAMELIST predecessor_subpathways;  // Associates a superpathway WID with a named subpathway that is a xpredecessor

/*
 * List of external references
 */
struct dblinks {
  char *database;
  char *xid;
  char *description;
  struct dblinks *next;
};

#define NULL_DBLINKS (struct dblinks *) NULL

/* theses are holdovers from KEGG and are not used in BioCyc */
struct dblinks * mk_dblink(char *, char *, struct stringlist *);
struct dblinks * mk_dblinks(char *, struct stringlist *);
struct dblinks * prepend_dblink(struct dblinks *, struct dblinks *);
void free_dblinks(struct dblinks *);

void
insert_dblinks(int object_wid, struct stringlist * dbs, struct stringlist * ids);

int
find_transunit(char *name);

int
find_pathway(char *name);

int
find_protein(char *name);

int
find_citation(char *citation);

short
insert_citation(int other_wid, char * citation);

short
insert_support(int other_wid, char *citation);

int
obtain_enzrxn_compound(char *name, int allow_protein);

                                                       
#define REFERENCE_ORGANISM (struct organisms **) NULL


/*
 * Supporting Functions
 */

void fatal_error (char *);

#define advise_error(args) fprintf(stderr, args)

#endif
