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
import com.sri.bw.dbschema.SchemaDocument;
import com.sri.bw.dbschema.SchemaDocument.Schema;
import com.sri.bw.dbschema.TableDocument.Table;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.log4j.Logger;
import java.io.File;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Jul 13, 2006
 */
public class DeleteDataSet extends LoaderMain {

    private static final Logger log = Logger.getLogger(DeleteDataSet.class);

    private String filename;
    private String datasetWID;
    private Warehouse warehouse;

    public static void main(String[] args) throws Exception {
        DeleteDataSet delete = new DeleteDataSet(args);
    }


    public DeleteDataSet(String[] args) throws Exception {
        super("DeleteDataSet", null, null);
        addOption(HELP_OPTION, false);
        addOption(PROPERTIES_FILE_OPTION, false);
        addOption(DB_HOST_OPTION, true);
        addOption(DB_NAME_OPTION, true);
        addOption(DBMS_OPTION, true);
        addOption(DB_USERNAME_OPTION, true);
        addOption(DB_PASSWORD_OPTION, true);
        addOption(DB_PORT_OPTION, true);
        addOption(SINGLE_INPUT_FILE_OPTION, true);


        Option dataSetListOption = OptionBuilder.withArgName("datasetwids")
                .hasArg()
                .withLongOpt("datasetwids")
                .withDescription("Comma-separated list of DataSet.WIDs of datasets to be deleted")
                .create('x');
        addOption(dataSetListOption, true);

        parseCommandLine(args, "runDelteDataSet.sh");

        if (properties.validate()) {
            filename = properties.getFile();
            String widList = properties.getProperty(dataSetListOption.getLongOpt());

            connectToWarehouse();
            warehouse = WarehouseManager.getWarehouse();

            String[] wids = widList.split(",");
            for (String wid : wids) {
                datasetWID = wid;
                deleteDataSet();
            }
        }
    }


    private void deleteDataSet() throws Exception {
        log.info("Deleting dataset " + datasetWID + "...");
        File inputFile = new File(filename);

        log.info("Reading schema from " + filename);

        Schema schema = SchemaDocument.Factory.parse(inputFile).getSchema();

        Vector<Table> objectTables = getTablesOfType(schema, "object");
        Vector<Table> descriptorTables = getTablesOfType(schema, "object-descriptor");

        String query = "select count(*) from Entry where DataSetWID=" + datasetWID;
        ResultSet results = warehouse.executeQuery(query);
        results.next();
        log.info(results.getInt(1) + " entries for DataSetWID " + datasetWID);
        if (results.getInt(1) == 0) {
            log.warn("No entries for DataSetWID " + datasetWID);
        }
        results.close();

        // drop all objects
        deleteFromObjectTables(objectTables);

        // Drop all in object descriptor
        deleteFromDescriptorTables(descriptorTables);

        // drop Entries
        query = "delete from Entry where DataSetWID=" + datasetWID;
        warehouse.executeUpdate(query);

        // drop dataset row
        query = "delete from DataSet where WID=" + datasetWID;
        warehouse.executeUpdate(query);

        warehouse.commit();

        reportTimeElapsed();
    }

    private void deleteFromObjectTables(Vector<Table> tables) {
        String query;

        for (Table table : tables) {
            log.info("Deleting from " + table.getName());
            try {
                query = "delete from " + table.getName() + " where DataSetWID=" + datasetWID;
                int result = warehouse.executeUpdate(query);
                if (result > 0) {
                    log.info("deleted " + result + " rows from " + table.getName());
                }
            } catch (SQLException e) {
                log.error("Error populating prepared statement " + table.getName() + ": " + e.getMessage());
            }
        }
    }

    private void deleteFromDescriptorTables(Vector<Table> tables) {
        String query;
        String subSelect = "select OtherWID from Entry where DataSetWID=" + datasetWID;

        for (Table table : tables) {
            log.info("Deleting from " + table.getName());
            try {
                query = "delete from " + table.getName() + " where OtherWID in (" + subSelect + ")";
                int result = warehouse.executeUpdate(query);
                if (result > 0) {
                    log.info("deleted " + result + " rows from " + table.getName());
                }
            } catch (SQLException e) {
                log.error("Error populating prepared statement " + table.getName() + ": " + e.getMessage());
            }
        }
    }


    private Vector<Table> getTablesOfType(Schema schema, String tableType) {
        Vector<Table> list = new Vector<Table>();
        for (Table table : schema.getTableArray()) {
            if (table.getType().equals(tableType) && !table.getName().equals("Entry")) {
                list.add(table);
            }
        }
        return list;
    }

}
