/* $Id: MySQLDatabase.java,v 1.1 2006/07/07 15:03:35 valerie Exp $ */
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
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Implements Database functionality specific to MySQL (driver name and JDBC connection string)
 * @author Valerie Wagner
 *         Date: Apr 15, 2004
 */
public class MySQLDatabase extends DatabaseProxy
{
    private static final Logger log = Logger.getLogger( MySQLDatabase.class );


    MySQLDatabase( String databaseHost, String databasePort, String databaseName )
            throws ClassNotFoundException
    {
        super( databaseHost, databasePort, databaseName );
    }


    /**
     * Initializes the driver name for MySQL.
     * Creates the database connection string for MySQL.
     */
    void initialize()
    {
        driverName = "com.mysql.jdbc.Driver";

        /*
        The JDBC URL format for MySQL Connector/J is as follows, with items in square brackets ([, ]) being optional:
        jdbc:mysql://[host][,failoverhost...][:port]/[database][?propertyName1][=propertyValue1][&propertyName2][=propertyValue2]...
        */

        databaseURL = "jdbc:mysql://" + databaseHost +":" + databasePort + "/" + databaseName;
        log.debug( "Mysql connection url = " + databaseURL );
    }


    /*
    public void updateClob( String tableName, String columnName, long WID, String value ) throws SQLException
    {
        String sqlLine = "SELECT " + columnName + " FROM " + tableName +
            " WHERE WID = " + WID + " for update nowait";
log.debug( "Inserting a Clob.  Selection statement is: " +sqlLine );

        try
        {
            ResultSet rs = executeQuery( sqlLine );

            if( rs.next() )
            {

                Clob clob = rs.getClob( columnName );
                String updateSQLLine = "UPDATE " + tableName + " SET " + columnName + " + ?";
log.debug( "Clob update line: " + updateSQLLine );
                PreparedStatement pstmt = createPreparedStatement( updateSQLLine );
                pstmt.setClob( 1, clob );
                pstmt.executeUpdate();
            }
            else
            {
                log.error( "There was no result set returned by the query: " + sqlLine );
            }

        }
        catch( SQLException e )
        {
            log.error( "Error when updating CLOB", e );
        }
    }
    */


    public int getDBMSType()
    {
        return Database.MYSQL;
    }


    public ResultSet getTableMetaData( String tableName )
    {
        ResultSet results = null;

        try
        {
            // The "%" character matches all column names
            results = connection.getMetaData().getColumns( null, null, tableName, "%" );
        }
        catch( SQLException e )
        {
            log.error( "" );
        }

        return results;
    }


}
