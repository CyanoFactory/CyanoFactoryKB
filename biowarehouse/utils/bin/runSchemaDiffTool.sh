#!/bin/sh

echo ${JAVA_HOME}

if [ "$JAVA_HOME" = "" ]; then
    echo 1>&2 Error: JAVA_HOME must be defined
    exit 127
fi

if [ "$WINDIR" = "" ]; then
    PATHSEP=':'
else
    PATHSEP=';'
fi

LOADER_DIR=.
LIB_DIR=../src/java//lib
CLIB_DIR=${LIB_DIR}
WH_DIR=../src/java/dist

CPJ=${CPJ}${PATHSEP}${WH_DIR}/biowarehouse.jar
CPJ=${CPJ}${PATHSEP}${CLIB_DIR}/jaxb-api.jar
CPJ=${CPJ}${PATHSEP}${CLIB_DIR}/jaxb-impl.jar
CPJ=${CPJ}${PATHSEP}${CLIB_DIR}/jsr173_api.jar


"${JAVA_HOME}/bin/java" -ms100M -mx512M -cp "${CPJ}" com.sri.biospice.warehouse.schema.tools.SchemaDiffTool $*
