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
 * Function prototypes, structs, global variables, etc...
 */

#ifndef _STRING_UTIL_H
#define _STRING_UTIL_H 1

#include <stddef.h>

void
string_uppercase (char *str);

char *
string_append (char *, char *);

char *
string_concat (char *, char *);

char *
string_colon_concat (char *, char *);

char *
string_concatenate(char *s1, char *s2, int space);

char *
string_truncate(char *str, int maxlen);

char * 
string_number(char * str);

/* Linked list of strings */
struct stringlist {
  char *string;
  struct stringlist *next;
};

#define NULL_LIST (struct stringlist *) NULL

void
stringlist_print(struct stringlist *list);

struct stringlist *
error_stringlist (void);

int
stringlist_find(char *string, struct stringlist *list);
							  
int
stringlist_length(struct stringlist *list);

struct stringlist *
stringlist_cons (char *string, struct stringlist *list);

// Return a list with N copies of the argument string
struct stringlist *
stringlist_rep(char *string, int num_copies);

struct stringlist *
stringlist_append (char *string, struct stringlist *list);

char *
flatten_list_old (struct stringlist *);

char *
flatten_list (struct stringlist *list, int space);

void
free_stringlist (struct stringlist *);

#define list_flatten_list(list) stringlist_cons(flatten_list(list), NULL_LIST)

#endif
