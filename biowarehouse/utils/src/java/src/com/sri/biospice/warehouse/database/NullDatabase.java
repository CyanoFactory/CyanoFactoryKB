/* $Id: NullDatabase.java,v 1.2 2008/10/10 22:02:19 valerie Exp $ */
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
package com.sri.biospice.warehouse.database;

import java.sql.*;
import java.math.BigDecimal;
import java.io.InputStream;
import java.io.Reader;
import java.util.Calendar;
import java.net.URL;

/**
 * @author Valerie Wagner
 *         Date: Jul 19, 2004
 */
public class NullDatabase implements Database {

    public void connectToDatabase(String username, String password) throws SQLException {
    }


    public ResultSet executeQuery(String query) throws SQLException {
        return null;
    }


    public int executeUpdate(String command) throws SQLException {
        return 0;
    }


    public void close() throws SQLException {
    }


    public void commit() {
    }


    public void rollback() {
    }


    public int getDBMSType() {
        return 0;
    }


    public String getJDBCDate() {
        return null;
    }


    public String getJDBCTimestamp() {
        return null;
    }


    public void setTimeoutLimit(int seconds) throws SQLException {
    }


    public String getVersion() {
        return null;
    }


    public ResultSet getTableMetaData(String tableName) {
        return null;
    }


    public PreparedStatement createPreparedStatement(String sqlLine) throws SQLException {
        PreparedStatement ps = new PreparedStatement() {

            public ResultSet executeQuery() throws SQLException {
                return null;
            }

            public int executeUpdate() throws SQLException {
                return 0;
            }

            public void setNull(int i, int i1) throws SQLException {

            }

            public void setBoolean(int i, boolean b) throws SQLException {

            }

            public void setByte(int i, byte b) throws SQLException {

            }

            public void setShort(int j, short i) throws SQLException {

            }

            public void setInt(int i, int i1) throws SQLException {

            }

            public void setLong(int i, long l) throws SQLException {

            }

            public void setFloat(int i, float v) throws SQLException {

            }

            public void setDouble(int i, double v) throws SQLException {

            }

            public void setBigDecimal(int i, BigDecimal bigDecimal) throws SQLException {

            }

            public void setString(int i, String s) throws SQLException {

            }

            public void setBytes(int i, byte[] bytes) throws SQLException {

            }

            public void setDate(int i, Date date) throws SQLException {

            }

            public void setTime(int i, Time time) throws SQLException {

            }

            public void setTimestamp(int i, Timestamp timestamp) throws SQLException {

            }

            public void setAsciiStream(int i, InputStream inputStream, int i1) throws SQLException {

            }

            public void setUnicodeStream(int i, InputStream inputStream, int i1) throws SQLException {

            }

            public void setBinaryStream(int i, InputStream inputStream, int i1) throws SQLException {

            }

            public void clearParameters() throws SQLException {

            }

            public void setObject(int i, Object o, int i1) throws SQLException {

            }

            public void setObject(int i, Object o) throws SQLException {

            }

            public boolean execute() throws SQLException {
                return false;
            }

            public void addBatch() throws SQLException {

            }

            public void setCharacterStream(int i, Reader reader, int i1) throws SQLException {

            }

            public void setRef(int i, Ref ref) throws SQLException {

            }

            public void setBlob(int i, Blob blob) throws SQLException {

            }

            public void setClob(int i, Clob clob) throws SQLException {

            }

            public void setArray(int i, Array array) throws SQLException {

            }

            public ResultSetMetaData getMetaData() throws SQLException {
                return null;
            }

            public void setDate(int i, Date date, Calendar calendar) throws SQLException {

            }

            public void setTime(int i, Time time, Calendar calendar) throws SQLException {

            }

            public void setTimestamp(int i, Timestamp timestamp, Calendar calendar) throws SQLException {

            }

            public void setNull(int i, int i1, String s) throws SQLException {

            }

            public void setURL(int i, URL url) throws SQLException {

            }

            public ParameterMetaData getParameterMetaData() throws SQLException {
                return null;
            }

            public void setRowId(int i, RowId rowId) throws SQLException {

            }

            public void setNString(int i, String s) throws SQLException {

            }

            public void setNCharacterStream(int i, Reader reader, long l) throws SQLException {

            }

            public void setNClob(int i, NClob nClob) throws SQLException {

            }

            public void setClob(int i, Reader reader, long l) throws SQLException {

            }

            public void setBlob(int i, InputStream inputStream, long l) throws SQLException {

            }

            public void setNClob(int i, Reader reader, long l) throws SQLException {

            }

            public void setSQLXML(int i, SQLXML sqlxml) throws SQLException {

            }

            public void setObject(int i, Object o, int i1, int i2) throws SQLException {

            }

            public void setAsciiStream(int i, InputStream inputStream, long l) throws SQLException {

            }

            public void setBinaryStream(int i, InputStream inputStream, long l) throws SQLException {

            }

            public void setCharacterStream(int i, Reader reader, long l) throws SQLException {

            }

            public void setAsciiStream(int i, InputStream inputStream) throws SQLException {

            }

            public void setBinaryStream(int i, InputStream inputStream) throws SQLException {

            }

            public void setCharacterStream(int i, Reader reader) throws SQLException {

            }

            public void setNCharacterStream(int i, Reader reader) throws SQLException {

            }

            public void setClob(int i, Reader reader) throws SQLException {

            }

            public void setBlob(int i, InputStream inputStream) throws SQLException {

            }

            public void setNClob(int i, Reader reader) throws SQLException {

            }

            public ResultSet executeQuery(String s) throws SQLException {
                return null;
            }

            public int executeUpdate(String s) throws SQLException {
                return 0;
            }

            public void close() throws SQLException {

            }

            public int getMaxFieldSize() throws SQLException {
                return 0;
            }

            public void setMaxFieldSize(int i) throws SQLException {

            }

            public int getMaxRows() throws SQLException {
                return 0;
            }

            public void setMaxRows(int i) throws SQLException {

            }

            public void setEscapeProcessing(boolean b) throws SQLException {

            }

            public int getQueryTimeout() throws SQLException {
                return 0;
            }

            public void setQueryTimeout(int i) throws SQLException {

            }

            public void cancel() throws SQLException {

            }

            public SQLWarning getWarnings() throws SQLException {
                return null;
            }

            public void clearWarnings() throws SQLException {

            }

            public void setCursorName(String s) throws SQLException {

            }

            public boolean execute(String s) throws SQLException {
                return false;
            }

            public ResultSet getResultSet() throws SQLException {
                return null;
            }

            public int getUpdateCount() throws SQLException {
                return 0;
            }

            public boolean getMoreResults() throws SQLException {
                return false;
            }

            public void setFetchDirection(int i) throws SQLException {

            }

            public int getFetchDirection() throws SQLException {
                return 0;
            }

            public void setFetchSize(int i) throws SQLException {

            }

            public int getFetchSize() throws SQLException {
                return 0;
            }

            public int getResultSetConcurrency() throws SQLException {
                return 0;
            }

            public int getResultSetType() throws SQLException {
                return 0;
            }

            public void addBatch(String s) throws SQLException {

            }

            public void clearBatch() throws SQLException {

            }

            public int[] executeBatch() throws SQLException {
                return new int[0];
            }

            public Connection getConnection() throws SQLException {
                return null;
            }

            public boolean getMoreResults(int i) throws SQLException {
                return false;
            }

            public ResultSet getGeneratedKeys() throws SQLException {
                return null;
            }

            public int executeUpdate(String s, int i) throws SQLException {
                return 0;
            }

            public int executeUpdate(String s, int[] ints) throws SQLException {
                return 0;
            }

            public int executeUpdate(String s, String[] strings) throws SQLException {
                return 0;
            }

            public boolean execute(String s, int i) throws SQLException {
                return false;
            }

            public boolean execute(String s, int[] ints) throws SQLException {
                return false;
            }

            public boolean execute(String s, String[] strings) throws SQLException {
                return false;
            }

            public int getResultSetHoldability() throws SQLException {
                return 0;
            }

            public boolean isClosed() throws SQLException {
                return false;
            }

            public void setPoolable(boolean b) throws SQLException {

            }

            public boolean isPoolable() throws SQLException {
                return false;
            }

            public <T> T unwrap(Class<T> tClass) throws SQLException {
                return null;
            }

            public boolean isWrapperFor(Class<?> aClass) throws SQLException {
                return false;
            }
        }        ;
        return ps;
    }


    public DatabaseMetaData getMetaData() throws SQLException {
        return null;
    }


    public String getSchemaName() {
        return null;
    }

    /**
     * Put database connection into read-only mode
     *
     * @param readOnly true to enable read-only mode; false disables read-only mode
     */
    public void setReadOnly(boolean readOnly) {
    }
}
