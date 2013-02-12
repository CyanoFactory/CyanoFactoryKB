/* $Id: DatabaseProxy.java,v 1.4 2007/08/29 23:44:02 nguo Exp $ */
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

import org.apache.log4j.Logger;
import java.sql.*;
import java.util.Vector;

/**
 * Implements basic Database functionality not specific
 * to any DBMS type.
 * @author Valerie Wagner
 *         Date: Apr 15, 2004
 */
public abstract class DatabaseProxy implements Database
{
    private static final Logger log = Logger.getLogger( DatabaseProxy.class );
    protected Connection connection;  // Connection to the database
    protected Statement stmt;         // Statement used for executing queries
    protected String driverName;
    protected String databaseURL;
    protected String databaseHost;
    protected String databaseName;
    protected String databasePort;
    protected String schemaName;
    private int timeoutLimit = Database.TIMEOUT_LIMIT;


    /**
     * Creates a class to act as a proxy for a remote database
     * @param databaseHost Machine name or IP address of database server
     * @param databaseName Name of database
     * @throws ClassNotFoundException if the driver for the database cannot be loaded
     */
    DatabaseProxy( String databaseHost, String databasePort, String databaseName ) throws ClassNotFoundException
    {
        this.databaseHost = databaseHost;
        this.databaseName = databaseName;
        this.databasePort = databasePort;

        initialize();
        loadDriver();
    }


    /**
     * Extending classes will specify the driver name and connection string
     * for this database connection.
     */
    abstract void initialize();


    /**
     * Loads the appropriate driver for database.
     */
    public void loadDriver() throws ClassNotFoundException
    {
        log.debug( "Attempting to load driver..." );
        Class tmpClass = Class.forName( driverName );
        log.debug( "driver class: " + tmpClass );

        if( tmpClass == null )
        {
            log.error( "Could not loader driver: " + driverName );
        }
        else
        {
            log.debug( "Driver successfully loaded" );
        }
    }


    public void connectToDatabase( String username, String password ) throws SQLException
    {
        // todo: if the user has registered the wrong driver for the system
        // you are trying to connect to, you get the very unhelpful message
        // about not being able to connect because connection refused
        // and the stack trace goes via the "unregistered driver" classes

        log.info( "Attempting connection to database with url=" + databaseURL + " and username=" + username );

        connection = DriverManager.getConnection( databaseURL, username, password );
        schemaName = username;
        log.debug( "Connection " + connection );
        stmt = connection.createStatement( ResultSet.TYPE_SCROLL_INSENSITIVE,
                                           ResultSet.CONCUR_UPDATABLE );
        stmt.setQueryTimeout( timeoutLimit );
        connection.setAutoCommit( false );
        schemaName = username.toUpperCase();

        log.info( "Connection successful!" );
        log.debug( "Statement query timeout = " + stmt.getQueryTimeout() + " seconds" );
        log.info("Autocommit = " + connection.getAutoCommit());
    }


    public ResultSet executeQuery( String query ) throws SQLException
    {
        log.debug( "Executing query: " + query );
        ResultSet result = stmt.executeQuery( query );
        return result;
    }


    public int executeUpdate( String command ) throws SQLException
    {
        log.debug( "Executing update: " + command );
        return stmt.executeUpdate( command );
    }


    public void close() throws SQLException
    {
        log.debug( "Closing the connection to the database" );
        connection.close();
    }


    public static String escapeSingleQuotes( String input )
    {
        if( input == null )
        {
            return null;
        }

        String output = input.replaceAll( "'", "''" );
        return output;
    }


    public void commit()
    {
        try
        {
            connection.commit();
        }
        catch( SQLException e )
        {
            log.error( "Error trying to perform commit on database", e );
        }
    }


    public void rollback()
    {
        try
        {
            connection.rollback();
        }
        catch( SQLException e )
        {
            log.error( "Error trying to perform rollback on database", e );
        }
    }


    public void setTimeoutLimit( int seconds ) throws SQLException
    {
        this.timeoutLimit = seconds;
        stmt.setQueryTimeout( timeoutLimit );
    }


    public String getVersion()
    {
        return "BioWarehouse database module.  Version number: @VERSION_NUMBER@, " +
               "Build number: @BUILD_NUMBER@";
    }



    public PreparedStatement createPreparedStatement( String sqlLine ) throws SQLException
    {
        return connection.prepareStatement( sqlLine );
    }


    public DatabaseMetaData getMetaData() throws SQLException
    {
        DatabaseMetaData metadata = connection.getMetaData();
        if( metadata == null )
        {
            throw new SQLException( "Database metadata could not be obtained" );
        }
        return metadata;
    }


    public String getSchemaName()
    {
        return schemaName;
    }

    public void setReadOnly(boolean readOnly) {
        try {
            connection.setReadOnly(readOnly);
        } catch (SQLException e) {
            log.error("Unable to create read-only connection.  Program can continue however.");
        }
    }


}
