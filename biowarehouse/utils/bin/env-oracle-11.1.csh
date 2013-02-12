#!/bin/csh
# Use cilantro's Oracle 11.1.0 installation 
setenv ORACLE_HOME /opt/oracle/ClientFDev
setenv ORACLE_SID biospice
# Do this until an 11.1 tnsnames.ora is in place, if ever
setenv TNS_ADMIN /obase/app/oracle/product/10.1.0/db_1/network/admin
setenv ORACLE_SID biospice
setenv PATH $ORACLE_HOME/bin:$PATH
if ($?LD_LIBRARY_PATH == 1) then
  setenv LD_LIBRARY_PATH	${ORACLE_HOME}/lib:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH	${ORACLE_HOME}/lib
endif
