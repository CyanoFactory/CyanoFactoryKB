/* $Id: DeletionFactory.java,v 1.1 2006/07/07 15:03:35 valerie Exp $ */
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
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.Hashtable;
import java.util.Vector;
import com.sri.biospice.warehouse.schema.linking.BioSourceWIDBioSubtypeWID;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Oct 1, 2004
 */
public class DeletionFactory
{
    private static final Logger log = Logger.getLogger( DeletionFactory.class );

    private static Warehouse warehouse = WarehouseManager.getWarehouse();
    private static Hashtable widRelatedDeletes = new Hashtable();
    private static Vector linkingTableDeletes = new Vector();
    private static final boolean IS_MYSQL = (warehouse.getDBMSType() == Database.MYSQL);


    private static PreparedStatement createStatement( String statement )
    {
        try
        {
            return warehouse.createPreparedStatement( statement );
        }
        catch( SQLException e )
        {
            log.error( "Could not create deletion statment: " + statement, e );
        }
        return null;
    }


    public static void deletObjectTable( String tableName, long wid )
    {
        deleteWIDRelatedTable( "WID", tableName, wid );
    }


    public static void deletOtherWIDTable( String tableName, long wid )
    {
        deleteWIDRelatedTable( "OtherWID", tableName, wid );
    }


    private static void deleteWIDRelatedTable( String columnName, String tableName, long wid )
    {
        if( widRelatedDeletes.get( tableName ) == null )
        {
            widRelatedDeletes.put( tableName, createStatement( "delete from " + tableName + " where " + columnName + "=?" ) );
        }

        PreparedStatement stmt = (PreparedStatement)widRelatedDeletes.get( tableName );
        try
        {
            stmt.clearParameters();
            stmt.setLong( 1, wid );
            stmt.executeUpdate();
        }
        catch( SQLException e )
        {
            log.error( "Error deleting table " + tableName, e );
        }
    }


    public static void deleteLinkingTableEntries( long wid )
    {
        if( IS_MYSQL )
        {
            if( linkingTableDeletes.isEmpty() )
            {
                createLinkingTableDeletes();
            }

            for( int i = 0; i < linkingTableDeletes.size(); i++ )
            {
                PreparedStatement stmt = (PreparedStatement)linkingTableDeletes.elementAt( i );
                try
                {
                    stmt.clearParameters();
                    stmt.setLong( 1, wid );
                    stmt.setLong( 2, wid );
                    stmt.executeUpdate();
                }
                catch( SQLException e )
                {
                    log.error( "Error deleting linking table", e );
                }
            }
        }
    }


    private static void createLinkingTableDeletes()
    {
        linkingTableDeletes.add( createLinkingTableDelete( "BioSourceWID", "BioSubtypeWID" ) );
        linkingTableDeletes.add( createLinkingTableDelete( "BioSourceWID", "GeneWID" ) );
        linkingTableDeletes.add( createLinkingTableDelete( "BioSourceWID", "ProteinWID" ) );
        linkingTableDeletes.add( createLinkingTableDelete( "CitationWID", "OtherWID" ) );
        linkingTableDeletes.add( createLinkingTableDelete( "GeneWID", "NucleicAcidWID" ) );
        linkingTableDeletes.add( createLinkingTableDelete( "GeneWID", "ProteinWID" ) );
        linkingTableDeletes.add( createLinkingTableDelete( "ProteinWID", "FunctionWID" ) );
    }


    private static PreparedStatement createLinkingTableDelete( String column1, String column2 )
    {
        String deleteCommand = "delete from " + column1 + column2 + " where " +
                               column1 + "=? or " + column2 + "=?";

        return createStatement( deleteCommand );
    }
}
