/* $Id: ForceDropTables.java,v 1.4 2008/10/10 22:02:20 valerie Exp $ */
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
import org.apache.log4j.Logger;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Jul 21, 2004
 */
public class ForceDropTables extends LoaderMain {
    private static final Logger log = Logger.getLogger(ForceDropTables.class);

    private int lastDroppedCount = 0;
    private Warehouse warehouse;
    private String dropCommandEnding = "";


    public static void main(String[] args) throws SQLException {
        new ForceDropTables(args);
    }


    public ForceDropTables(String[] args) throws SQLException {
        super(args, "ForceDropTables", "no-script");
        connectToWarehouse();
        warehouse = WarehouseManager.getWarehouse();
        if (warehouse.getDBMSType() == Warehouse.ORACLE) {
            dropTablesOracle();
        } else {
            dropTablesMySQL();
        }
    }


    public void dropTablesOracle() throws SQLException {
        String tableName;
        Vector<String> tableNames = new Vector<String>();

        dropCommandEnding = "cascade constraints purge";

        String findUndroppedTables = "select table_name from user_tables where Dropped='NO'";
        ResultSet undroppedTables = warehouse.executeQuery(findUndroppedTables);
        while (undroppedTables.next()) {
            tableName = undroppedTables.getString("table_name");
            if (tableName != null) {
                tableNames.add(tableName);
            }
        }
        undroppedTables.close();
        doTheDrop(tableNames);
        dropOracleSequences();
    }

    private void dropOracleSequences() {
        try {
            warehouse.executeUpdate("drop SEQUENCE WID_sequence");
        } catch (SQLException e) {
            log.warn(e.getMessage());
        }
        try {
            warehouse.executeUpdate("drop SEQUENCE SpecialWID_sequence");
        } catch (SQLException e) {
            log.warn(e.getMessage());
        }
    }

    public void dropTablesMySQL() throws SQLException {
        String tableName;
        Vector<String> tableNames = new Vector<String>();

        // todo: figure out how to not require this command
        String fkCmd = "SET FOREIGN_KEY_CHECKS = 0";
        warehouse.executeUpdate(fkCmd);

        String findUndroppedTables = "show tables";
        ResultSet undroppedTables = warehouse.executeQuery(findUndroppedTables);
        while (undroppedTables.next()) {
            tableName = undroppedTables.getString(1);
            if (tableName != null) {
                tableNames.add(tableName);
            }
        }
        undroppedTables.close();
        doTheDrop(tableNames);

        fkCmd = "SET FOREIGN_KEY_CHECKS = 1";
        warehouse.executeUpdate(fkCmd);
    }

    private void doTheDrop(Vector<String> tableNames) {
        dropTablesInList(tableNames);
        warehouse.commit();
    }


    private void dropTablesInList(Vector<String> tableNames) {
        String tableName, dropCommand;

        if (tableNames.size() == 0) {
            log.debug("Finished.  All tables dropped.");
            return;
        } else if (tableNames.size() == lastDroppedCount) {
            log.debug("Unable to drop any tables in the last round.  Giving up.");
            return;
        }

        lastDroppedCount = tableNames.size();
        log.debug("Current drop count: " + lastDroppedCount);

        for (int i = 0; i < tableNames.size(); i++) {
            tableName = tableNames.elementAt(i);
            dropCommand = "drop table " + tableName + " " + dropCommandEnding;
//            dropCommand = "drop table " + tableName + " cascade constraints purge";

            try {
                log.debug("Trying to drop " + tableName);
                warehouse.executeUpdate(dropCommand);
                tableNames.removeElementAt(i);
            } catch (SQLException e) {
                log.debug("Drop failed: " + e.getMessage().trim());
                System.exit(1);
            }
        }

        dropTablesInList(tableNames);
    }


    public String getLoaderVersionNumber() {
        return null;
    }


    public String getLoaderBuildNumber() {
        return null;
    }

    public float getEarliestSupportedDataVersion() {
        return 0;
    }

    public float getLatestSupportedDataVersion() {
        return 0;
    }

}
