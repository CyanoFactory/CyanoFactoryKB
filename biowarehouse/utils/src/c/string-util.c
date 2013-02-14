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
 * String Handling Utility functions for the warehouse loader
 */

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "string-util.h"

/////////////////////////////////////////////////////////////////////////////
//// Basic utilities
/////////////////////////////////////////////////////////////////////////////

void
string_uppercase (char *str) {
  /* Convert str in place to upper case */
  char *strp = str;
  while (*strp = toupper(*strp))
    strp++;
}

char *
string_truncate(char *str, int maxlen) {
  /* Return a string str or a substring of it truncated to length maxlen. */
  char *newstr;
  if (!str)
    return(str);
  if (strlen(str) <= maxlen)
    return(str);
  /* Allocate a new string, copy initial substring of old to it */
  if (maxlen < 0) maxlen = 0;  /* To be safe */
  newstr = (char *) malloc(maxlen+1);
  if (!newstr)
    fatal_error("malloc failed in string_truncate");
  strncpy(newstr, str, maxlen);
  newstr[maxlen] = '\0';
  /* free(str); */ /****!!! memory leak */
  return(newstr);
}

char * 
string_number (char * str) {
  char * num_str; 
  int index =0;
  if(!str)
    return (str);
  num_str = (char *)malloc(strlen(str) +1);
  if(!num_str)
	fatal_error("malloc failed in string_num");
  while((str[index] != '\0') && (str[index] >= '0' && str[index] <= '9'))
  {
	num_str[index] = str[index];
	index++;
  } 
  num_str[index] = '\0';
  return num_str;
}


/////////////////////////////////////////////////////////////////////////////
//// Concatenate and append
////
//// NOTE: these fns free its string arguments,
////  so make sure they are freeable (eg. via strdup() )
/////////////////////////////////////////////////////////////////////////////

char *
string_concatenate(char *s1, char *s2, int space) {
/* 
   space = 0 -> no separator
   space = 1 -> separate by " ",
   space = 2 -> separate by ", "
   space = 3 -> separate by " : "
   space = 4 -> separate by tab
*/
  
  char *new_string, *t;

  if (!s1 || *s1 == '\0') return s2;
  if (!s2 || *s2 == '\0') return s1;
  
  new_string = (char *) malloc (strlen(s1) + strlen(s2) + 1 + space); /* tab wastes 3 bytes */
  if (!new_string)
    fatal_error("malloc failed in string_concatenate");

  (void) strcpy(new_string, s1);

  if (space > 0) {
    t = new_string + strlen(s1); /* position t at end of str1 */
    if (space == 1) {
      *t++ = ' '; 
    }
    else if (space == 2) {
      *t++ = ','; *t++ = ' ';
    }
    else if (space == 3) {
      *t++ = ' '; *t++ = ':'; *t++ = ' ';
    }
    else if (space == 4) {
      *t++ = '\t'; 
    }
    *t = '\0';
  }
  (void) strcat(new_string, s2);
  free(s1); free(s2);
  return(new_string);
}

char *
string_append(char *s1, char *s2) {
  return(string_concatenate(s1, s2, 0));
}

char *
string_concat(char *s1, char *s2) {
  return(string_concatenate(s1, s2, 1));
}


char *
string_colon_concat(char *s1, char *s2) {
  return(string_concatenate(s1, s2, 3));
}


/////////////////////////////////////////////////////////////////////////////
//// Stringlist fns
/////////////////////////////////////////////////////////////////////////////

void
stringlist_print(struct stringlist *list) {
  // For debugging
  struct stringlist *slp = list;
  assert(list);

  printf(">> stringlist=\n");
  for (; slp; slp=slp->next)
    printf("   %s\n", slp->string);
}

char *
flatten_list_old (struct stringlist *list) {
  /* Old version, unused(?) in any loaders */
  
  struct stringlist *ptr;
  char *newstring;

  if (!list)
    return NULL;

  newstring = list->string;
  
  for (ptr=list->next;ptr;ptr=ptr->next) 
    newstring = string_concat(ptr->string, newstring);
  
  return(newstring);
}

void
free_stringlist_shallow(struct stringlist *list) {
  // Frees stringlist structs but not the strings referred to from it
  if (list) {
    free_stringlist_shallow(list->next);
    free(list);
  }
}

void
free_stringlist(struct stringlist *list) {
  if (list) {
    free_stringlist(list->next);
    free(list->string);
    free(list);
  }
}

char *
flatten_list (struct stringlist *list, int space) {
  /* KEGG version */
  
  struct stringlist *ptr;
  char *newstring;

  if (!list)
    return NULL;

  newstring = strdup(list->string); // tomlee added strdup() 2009-05-27
  
  for (ptr=list->next;ptr;ptr=ptr->next) 
    newstring = string_concatenate(ptr->string, newstring, space);

  free_stringlist_shallow(list);
  return(newstring);
}

struct stringlist *
error_stringlist(void) {

  struct stringlist *newlist;

  newlist = (struct stringlist *) malloc (sizeof(struct stringlist));
  if (!newlist)
    fatal_error("malloc failed in error_stringlist");

  newlist->next   = NULL_LIST;
  newlist->string = (char *) strdup("ERROR");
  return(newlist);
}

int
stringlist_find(char *string, struct stringlist *list) {
  /* Returns 0 unless string is in list. Comparison is case-INsensitive. */
  struct stringlist *slp = list;
  while (slp) {
    if (0 == strcasecmp(string, slp->string)) return 1;
    slp = slp->next;
  }
  return 0;
}

int
stringlist_length(struct stringlist *list) {
  int len = 0;
  struct stringlist *slp = list;
  while (slp) {
    len++;
    slp = slp->next;
  }
  return len;
}

struct stringlist *
stringlist_cons(char *string, struct stringlist *list) {

  struct stringlist *newlist;

  if (!string || *string == '\0') return list;

  newlist = (struct stringlist *) malloc (sizeof(struct stringlist));
  if (!newlist)
    fatal_error("malloc failed in stringlist_cons");

  newlist->next   = list;
  newlist->string = string;
  return(newlist);
}

// Return a list with N strdup'ed copies of the argument string
struct stringlist *
stringlist_rep(char *string, int num_copies) {

  struct stringlist *newlist = NULL_LIST;
  int j;
  
  if (!string || *string == '\0') return newlist;
  
  for (j=0; j<num_copies; j++)
    newlist = stringlist_cons(strdup(string), newlist);
  return(newlist);
}

struct stringlist *
stringlist_append(char *string, struct stringlist *list) {
  
  struct stringlist *t = list;
  assert(list);

  for (;t->next;t=t->next);
  t->next = stringlist_cons(string,NULL_LIST);
  
  return(list);
}

