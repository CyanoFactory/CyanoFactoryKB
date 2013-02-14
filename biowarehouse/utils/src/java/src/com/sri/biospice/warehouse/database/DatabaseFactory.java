/* $Id: DatabaseFactory.java,v 1.1 2006/07/07 15:03:35 valerie Exp $ */
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

/**
 * Factory class to create general Database instances.
 * @author Valerie Wagner
 *         Date: Apr 15, 2004
 */
public class DatabaseFactory
{
    private static final Logger log = Logger.getLogger( DatabaseFactory.class );

    private static Database currentProxy = null;


    /**
     * Calls <code>createDatabase( int databaseType, String databaseHost, String databaseName )</code>
     * with named constant
     * @param databaseType Name of database type ("mysql" or "oracle")
     * @param databaseHost Hostname or IP address of database server
     * @param databaseName Database name (SID in the case of Oracle)
     * @return A proxy connection to the database or null if databaseType is unsupported or the JDBC driver could not be loaded
     */
    public static Database createDatabase( String databaseType, String databaseHost, String databasePort, String databaseName )
    {
        if( databaseType.equalsIgnoreCase( Database.MYSQL_STR ) )
        {
            return DatabaseFactory.createDatabase( Database.MYSQL, databaseHost, databasePort, databaseName );
        }
        else if( databaseType.equalsIgnoreCase( Database.ORACLE_STR ) )
        {
            return DatabaseFactory.createDatabase( Database.ORACLE, databaseHost, databasePort, databaseName );
        }
        else if( databaseType.equalsIgnoreCase( Database.NULL_STR ) )
        {
            return DatabaseFactory.createDatabase( Database.NULL, databaseHost, databasePort, databaseName );
        }
        else
        {
            return null;
        }
    }


    /**
     * Use this method to get an instance of a MYSQL or ORACLE database proxy.
     * @param databaseType Currently supports MYSQL and ORACLE
     * @param databaseHost Host machine name or IP address
     * @param databaseName Name of the database
     * @return a <code>DatabaseProxy</code> class for a specific database type
     */
    public static Database createDatabase( int databaseType, String databaseHost, String databasePort, String databaseName )
    {

        try
        {
            switch( databaseType )
            {
                case Database.MYSQL:
                    currentProxy = new MySQLDatabase( databaseHost, databasePort, databaseName );
                    break;
                case Database.ORACLE:
                    currentProxy = new OracleDatabase( databaseHost, databasePort, databaseName );
                    break;
                case Database.NULL:
                    currentProxy = new NullDatabase();
                    break;
            }
        }
        catch( ClassNotFoundException e )
        {
            log.error( "Class not found exception results from not being able to load the database driver.", e );
        }

        return currentProxy;
    }


    /**
     * Use this method to get a handle on the currently open database
     * @return The current instance of the databaseProxy
     */
    public static Database getDatabase()
    {
        return currentProxy;
    }
}
