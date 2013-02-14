This directory contains files that are common to two or more
C-based loaders for the BioSPICE Warehouse.

SOURCE FILES

wh.{h,c}::
	Specify functions used by all loaders for all DBMSes.
	All functions are named wh_NAME of the form: 
	  #ifdef ORACLE
	    wh_oracle_NAME(...);
	  #elif DEF_MYSQL
	    wh_mysql_NAME(...);
	  #endif
	  
wh_mysql_util.{h,c},
wh_oracle_util.{h,pc,c}:
	Implement the DBMS-specific functions that are called from wh.c.
	All fns are named wh_mysql... or wh_oracle... respectively.

widtable.{h,c}:
	Hash table to map a string name to a Warehouse ID.

string-util.{h,c}:
	Various utilities for string manipulation.
