/* $Id: WarehouseManager.java,v 1.1 2006/07/07 15:03:36 valerie Exp $ */
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
 * Factory class for creating instances of a connection to the
 * BioSPICE Warehouse database.  It also enforces the Singleton policy of having only
 * one instance of the Warehouse at any given time.
 * @author Valerie Wagner
 *         Date: Jun 11, 2004
 */
public class WarehouseManager
{
    private static final Logger log = Logger.getLogger( WarehouseManager.class );

    private static Warehouse theWarehouse;


    /**
     * Calls <code>initWarehouse( int databaseType, String databaseHost, String databaseName )</code>
     * with named constant
     * @param databaseType Name of database type ("mysql" or "oracle")
     * @param databaseHost Hostname or IP address of database server
     * @param databaseName Database name (SID in the case of Oracle)
     * @return A proxy connection to the Warehouse or null if databaseType is unsupported or the JDBC driver could not be loaded
     */
    public static Warehouse initWarehouse( String databaseType, String databaseHost, String databaseName, String databasePort )
    {
        if( databaseType.equalsIgnoreCase( Database.MYSQL_STR ) )
        {
            return WarehouseManager.initWarehouse( Database.MYSQL, databaseHost, databaseName, databasePort );
        }
        else if( databaseType.equalsIgnoreCase( Database.ORACLE_STR ) )
        {
            return WarehouseManager.initWarehouse( Database.ORACLE, databaseHost, databaseName, databasePort );
        }
        else if( databaseType.equalsIgnoreCase( Database.NULL_STR ) )
        {
            return WarehouseManager.initWarehouse( Database.NULL, databaseHost, databaseName, databasePort );
        }
        else
        {
            return null;
        }
    }


    /**
     * Use this method to get an instance of a MYSQL or ORACLE database proxy to the Warehouse.
     * @param databaseType Currently supports MYSQL and ORACLE
     * @param databaseHost Host machine name or IP address
     * @param databaseName Name of the database
     * @return a <code>Warehouse</code> instnace for a specific database type
     */
    public static Warehouse initWarehouse( int databaseType, String databaseHost, String databaseName, String databasePort )
    {

        try
        {
            switch( databaseType )
            {
                case Database.MYSQL:
                    theWarehouse = new MySQLWarehouse( databaseHost, databaseName, databasePort );
                    break;
                case Database.ORACLE:
                    theWarehouse = new OracleWarehouse( databaseHost, databaseName, databasePort );
                    break;
                case Database.NULL:
                    theWarehouse = new NullWarehouse();
                    break;

            }
        }
        catch( ClassNotFoundException e )
        {
            log.error( "Class not found exception results from not being able to load the database driver.", e );
        }

        return theWarehouse;
    }


    /**
     * @return The current instance of the Warehouse
     */
    public static Warehouse getWarehouse()
    {
        return theWarehouse;
    }

}
