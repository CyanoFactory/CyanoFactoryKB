/* $Id: Database.java,v 1.1 2006/07/07 15:03:35 valerie Exp $ */
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

import java.sql.DatabaseMetaData;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Defines the interaction with a Database
 * @author Valerie Wagner
 *         Date: Jun 14, 2004
 */
public interface Database
{
    public static final int MYSQL = 0;
    public static final int ORACLE = 1;
    public static final int NULL = 2;

    public static final String MYSQL_STR = "mysql";
    public static final String ORACLE_STR = "oracle";
    public static final String NULL_STR = "null";

    /**
     * The limit in seconds for how long a JDBC statement can
     * attempt to execute a query on the database.
     * Currently equals 600 seconds (10 minutes).
     */
    public static final int TIMEOUT_LIMIT = 600;


    /**
     * Connect to a database using the username and password
     * @param username Username for database
     * @param password Password for user
     */
    void connectToDatabase( String username, String password ) throws SQLException;


    /**
     * Sends a query to the database and returns the ResultSet
     * @param query An SQL String.  Note: omit semicolon from query string
     * @return The result set from the query
     */
    ResultSet executeQuery( String query ) throws SQLException;


    /**
     * Sends an update (INSERT, MODIFY, or DELETE) command.
     * @param command The command to execute
     * @return The number of rows modified by the command
     * @throws java.sql.SQLException
     */
    int executeUpdate( String command ) throws SQLException;


    /**
     * Closes the connection to the database.
     * @throws java.sql.SQLException
     */
    void close() throws SQLException;


    /**
     * Commit all changes to the database.
     */
    void commit();


    /**
     * Rollback the changes since the last commit.
     */
    void rollback();


    /**
     * @return The type of DBMS (<code>ORACLE</code> or <code>MYSQL</code>)
     */
    int getDBMSType();


    /**
     * Sets the timeout limit for an SQL Statement while performing a query.
     * The timeout starts at <code>Database.TIMEOUT_LIMIT</code>.  Set
     * the limit to zero seconds to make there be no time limit.
     * @param seconds The number of seconds to wait for a query to execute
     * @throws SQLException
     */
    void setTimeoutLimit( int seconds ) throws SQLException;


    /**
     * @return the version number and build number of this version
     *         of the database utilities.
     */
    String getVersion();


    /**
     * Returns the metadata about a single table.
     * This operation incurs a large overhead and should be done
     * only once per table per program.  For more information,
     * see the <code>java.sql.DatabaseMetaData.getColumns()</code> javadoc in the
     * Java API.
     * @param tableName The name of the table in the Warehouse
     * @return The metadata associated with this table
     */
    ResultSet getTableMetaData( String tableName );


    /**
     * Creates a prepared (pre-compiled) SQL statement.
     * @param sqlLine The line from which the statement will be created.  Should be of the form
     *                "insert into tablename (col1, col2, ...) values (?, ?, ...)"
     * @return A prepared statement that may be used repeatedly to execute the sqlLine
     * @throws SQLException
     */
    PreparedStatement createPreparedStatement( String sqlLine ) throws SQLException;


    /**
     * @return The meta data associated with this database connection.
     * @throws SQLException
     */
    public DatabaseMetaData getMetaData() throws SQLException;


    /**
     * @return The name of the schema for this instance of the database
     */
    public String getSchemaName();

    /**
     * Put database connection into read-only mode
     * @param readOnly true to enable read-only mode; false disables read-only mode
     */
    public void setReadOnly(boolean readOnly);

}
