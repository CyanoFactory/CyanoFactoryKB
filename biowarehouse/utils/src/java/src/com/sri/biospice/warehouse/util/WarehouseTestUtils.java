/* $Id: WarehouseTestUtils.java,v 1.1 2006/07/07 15:03:38 valerie Exp $ */
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
package com.sri.biospice.warehouse.util;

import com.sri.biospice.warehouse.database.Warehouse;
import com.sri.biospice.warehouse.database.WarehouseManager;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import java.sql.SQLException;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Jul 7, 2004
 */
public class WarehouseTestUtils {
    private static final Logger log = Logger.getLogger(WarehouseTestUtils.class);

    public static String databaseHost;
    public static String databaseName;
    public static String databasePort;
    public static String username;
    public static String password;
    public static String databaseDBMSType;
    private static boolean initialized = false;


    /**
     * Sets up the connection to the warehouse, inits log4j,
     * gets database connection parameters from system properties,
     * which must be set in the Ant build file.
     *
     * @throws SQLException
     */
    public static void initWarehouseTest() throws SQLException {
        if (!initialized) {
            // Configure log4j
            BasicConfigurator.configure();

            // Get the properties for connecting to the database
            // These should be set by the Ant build file
            databaseHost = System.getProperty("warehouse.test.database.host");
            databaseName = System.getProperty("warehouse.test.database.name");
            databasePort = System.getProperty("warehouse.test.database.port");
            username = System.getProperty("warehouse.test.database.username");
            password = System.getProperty("warehouse.test.database.password");
            databaseDBMSType = System.getProperty("warehouse.test.database.dbms.type");
            log.debug("Connecting to " + username);
            // Connect to the Warehosue
            Warehouse warehouse = WarehouseManager.initWarehouse(databaseDBMSType, databaseHost, databaseName, databasePort);
            warehouse.connectToDatabase(username, password);
            initialized = true;
        }
    }


    /**
     * Closes the connection to the Warehouse.
     *
     * @throws SQLException
     */
    public static void endWarehouseTest() throws SQLException {
        WarehouseManager.getWarehouse().commit();
        WarehouseManager.getWarehouse().close();
    }


    /**
     * List of object tables.
     * And object table is a table having both a "WID" and "DataSetWID" column
     *
     * @return
     */
    public static Vector<String> getObjectTables() {

        // todo: automatically determine these
        Vector<String> tables = new Vector<String>();

        tables.add("Archive");
        tables.add("BioSource");
        tables.add("BioSubtype");
        tables.add("Chemical");
        tables.add("Citation");
        tables.add("Computation");
        tables.add("Division");
        tables.add("EnzymaticReaction");
        tables.add("Experiment");
        tables.add("ExperimentData");
        tables.add("Feature");
        tables.add("Function");
        tables.add("Gene");
        tables.add("GeneticCode");
        tables.add("NucleicAcid");
        tables.add("Pathway");
        tables.add("Protein");
        tables.add("Reaction");
        tables.add("Subsequence");
        tables.add("Support");
        tables.add("Taxon");
        tables.add("Term");

        return tables;
    }

    /**
     * A list of meta-table names
     * A meta-table is a table having an "OtherWID" column that is a foreign
     * key to an object table having a "WID" column.  These tables describe
     * objects in the object tables.   This list excludes the Entry table;
     *
     * @return
     */
    public static Vector<String> getMetaTables() {
        Vector<String> tables = new Vector<String>();

        tables.add("CommentTable");
        tables.add("CrossReference");
        tables.add("DBID");
        tables.add("Description");
        tables.add("Support");
        tables.add("SynonymTable");
        tables.add("ToolAdvice");

        return tables;
    }
}
