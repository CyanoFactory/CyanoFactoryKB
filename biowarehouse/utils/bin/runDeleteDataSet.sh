#!/bin/sh

if [ "${JAVA_HOME}" = "" ]; then
    echo 1>&2 Error: JAVA_HOME must be defined
    exit 127
fi

if [ "$WINDIR" = "" ]; then
    PATHSEP=':'
else
    PATHSEP=';'
fi

LOADER_DIR=.
LIB_DIR=../src/java/lib
CLIB_DIR=${LIB_DIR}
WH_DIR=../src/java/dist

CPJ=${CPJ}${PATHSEP}${WH_DIR}/biowarehouse.jar
CPJ=${CPJ}${PATHSEP}${CLIB_DIR}/db-schema.jar
CPJ=${CPJ}${PATHSEP}${CLIB_DIR}/element-map-beans.jar
CPJ=${CPJ}${PATHSEP}${CLIB_DIR}/xbean.jar
CPJ=${CPJ}${PATHSEP}${CLIB_DIR}/commons-cli-1.0.jar
CPJ=${CPJ}${PATHSEP}${CLIB_DIR}/log4j-1.2.8.jar
CPJ=${CPJ}${PATHSEP}${CLIB_DIR}/ojdbc14.jar
CPJ=${CPJ}${PATHSEP}${CLIB_DIR}/mysql-connector-java-3.0.9-stable-bin.jar

"${JAVA_HOME}/bin/java" -ms100M -mx512M -cp "${CPJ}" com.sri.biospice.warehouse.util.DeleteDataSet -f ../../schema/all-schema.xml $*

