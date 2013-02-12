/* $Id */
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

import com.sri.biospice.warehouse.schema.DataSet;
import com.sri.biospice.warehouse.schema.object.Protein;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Simple example of how to use the Warhouse to insert a new data set.
 * @author Valerie Wagner
 */
public class WarehouseExample
{
    private static final Logger log = Logger.getLogger( WarehouseExample.class );


    public static void main( String[] args )
    {
        String databaseDBMS = "oracle"; // or "mysql"
        String databaseHost = "localhost";
        String databasePort = "1523";
        String databaseName = "biospice";
        String databaseUsername = "user";
        String databasePassword = "pass";

        // Set up a simple log4j configuration that logs to the console.
        BasicConfigurator.configure();

        log.info( "Starting Warehouse Example..." );

        // Create an instance of the Warehouse, given the type (here, MySQL)
        Warehouse warehouse = WarehouseManager.initWarehouse( databaseDBMS,
                                                              databaseHost, databaseName, databasePort );

        // Try connecting to the Warehouse
        try
        {
            log.debug( "Trying to connect to database" );
            warehouse.connectToDatabase( databaseUsername, databasePassword );
        }
        catch( SQLException e )
        {
            log.error( "Error trying to connect to warehouse", e );
            System.exit( 1 );
        }

        // First create an entry in the DataSet table
        DataSet myDataSet = new DataSet( "My Dataset" );
        try
        {
            myDataSet.store();
        }
        catch( SQLException e )
        {
            log.error( "Error storing DataSet", e );
            System.exit( 1 );
        }

        // Now enter data into a table.  For example, Protein
        com.sri.biospice.warehouse.schema.object.Protein protein = new Protein( myDataSet.getWID() );
        protein.setName( "My protein's name" );
        protein.addSynonym( "Synonym for my protein" );
        try
        {
            protein.store();
        }
        catch( SQLException e )
        {
            log.error( "Error storing Protein", e );
        }

        // Load a protein already stored in the database
        // using the WID of the protein we stored above
        com.sri.biospice.warehouse.schema.object.Protein storedProtein = new com.sri.biospice.warehouse.schema.object.Protein( myDataSet.getWID(), protein.getWID() );
        try
        {
            storedProtein.load();
        }
        catch( SQLException e )
        {
            log.error( "Could not load stored Protein", e );
        }
        // Now you can access the values of Protein
        String name = storedProtein.getName();
        log.info( "Name of Protein: " + name );

        // To update the values of a row in a table, it is not necessary
        // to load it first (unless you want to).
        com.sri.biospice.warehouse.schema.object.Protein updatedProtein = new com.sri.biospice.warehouse.schema.object.Protein( myDataSet.getWID(), protein.getWID() );
        updatedProtein.setName( "New Protein Name" );
        try
        {
            updatedProtein.update();
        }
        catch( SQLException e )
        {
            log.error( "Error updating Protein table", e );
        }

        // Do a query on the Warehouse
        String query = "select * from DataSet";
        try
        {
            ResultSet results = warehouse.executeQuery( query );
            // Do whatever with your results
            // Close the ResultSet or there will be a resource leak!
            results.close();
        }
        catch( SQLException e )
        {
            log.error( "Error trying to execute query", e );
        }

        // Commit the changes and close the connection to the Warehouse
        log.info( "Committing changes" );
        warehouse.commit();
        try
        {
            warehouse.close();
        }
        catch( SQLException e )
        {
            log.error( "Error while trying to close connect to the database", e );
        }

        log.info( "Program finished" );
    }
}
