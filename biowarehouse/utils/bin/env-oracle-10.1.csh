#!/bin/csh
setenv ORACLE_HOME /u01/app/oracle/product/10.1.0/db_1 
setenv ORACLE_SID biospice
setenv PATH $ORACLE_HOME/bin:$PATH
if ($?LD_LIBRARY_PATH == 1) then
  setenv LD_LIBRARY_PATH	${ORACLE_HOME}/lib:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH	${ORACLE_HOME}/lib
endif
