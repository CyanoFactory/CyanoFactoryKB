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
import com.sri.bw.dbschema.TableDocument;
import org.apache.log4j.Logger;
import org.apache.xmlbeans.XmlException;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Hashtable;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Aug 14, 2006
 */
public class SummarizeDataSets extends LoaderMain {

    private static final Logger log = Logger.getLogger(SummarizeDataSets.class);
    private String summary = "\n\n";

    public SummarizeDataSets(String [] args) {
        super(args, "SummarizeDataSets", null);
    }

    public static void main(String[] args) throws IOException, XmlException {
        SummarizeDataSets loader = new SummarizeDataSets(args);
        loader.run();
    }

    private void run() throws IOException, XmlException {
        if (!properties.validate()) {
            return;
        }
        connectToWarehouse();
        String filename = properties.getFile();

        log.info("Using schema file: " + filename);
        SchemaDocument schemaDoc = SchemaDocument.Factory.parse(new File(filename));
        SchemaDocument.Schema schema = schemaDoc.getSchema();

        String query = "select * from DataSet";
        Warehouse warehouse = WarehouseManager.getWarehouse();
        Hashtable<String, StringBuffer> tableCounts = new Hashtable<String, StringBuffer>();

        for (TableDocument.Table table : getTablesOfType(schema, "object")) {
            tableCounts.put(table.getName(), new StringBuffer(String.format("%-28s,", table.getName())));
        }

        try {
            ResultSet results = warehouse.executeQuery(query);
            printDataSetTable(results);


            Vector<DataSetStruct> datasetList = new Vector<DataSetStruct>();
            results.beforeFirst();
            while (results.next()) {
                DataSetStruct ds = new DataSetStruct(
                        results.getString("name"),
                        results.getString("wid"),
                        results.getString("version"),
                        results.getString("changeDate"),
                        results.getString("loadDate"),
                        results.getString("loadedBy")
                );
                datasetList.add(ds);
            }
            results.close();

            for (int i = 0; i < datasetList.size(); i++) {

                String datasetwid = datasetList.get(i).wid;

                for (TableDocument.Table table : getTablesOfType(schema, "object")) {
                    query = "select count(*) from " + table.getName() + " where DataSetWID=" + datasetwid;
                    try {
                        ResultSet countResults = warehouse.executeQuery(query);
                        if (countResults.next()) {
                            StringBuffer buf = tableCounts.get(table.getName());
                            String text = String.format("%-18s,", countResults.getString(1));
                            buf.append(text);
                        }
                        countResults.close();
                    } catch (SQLException e) {
                        //ignore
                    }
                }
            }
            String text = String.format("%-28s,", "OBJECT TABLES / DATASETS");
            summary += text;
             for (int j = 0; j < datasetList.size(); j++) {
                text = String.format("%18s,", datasetList.get(j).name);
                summary += text;
            }
            println();

            text = String.format("%-28s,", "   DATASET VERSION");
            summary += text;

            for (int j = 0; j < datasetList.size(); j++) {
                text = String.format("%18s,", datasetList.get(j).version);
                summary += text;
            }
            println();

            text = String.format("%-28s,", "   DATASET WID");
            summary += text;
            for (int j = 0; j < datasetList.size(); j++) {
                text = String.format("%18s,", datasetList.get(j).wid);
                summary += text;
            }
            println();

            for (String key : tableCounts.keySet()) {
                StringBuffer buf = tableCounts.get(key);
                summary += buf + "\n";
            }
            println();

            // Consistency checks:
            query = "select distinct(WID) from DataSet order by WID";
            results = warehouse.executeQuery(query);
            summary += "\nDataSet.WIDs: ";
            while (results.next()) {
                summary += results.getString(1) + " ";
            }

            query = "select distinct(DataSetWID) from Entry order by DataSetWID";
            results = warehouse.executeQuery(query);
            summary += "\nEntry.DataSetWIDs: ";
            while (results.next()) {
                summary += results.getString(1) + " ";
            }

        } catch (SQLException e) {
            log.fatal(e);
            fail("Error executing query");
        }

        BufferedWriter out = new BufferedWriter(new FileWriter("datasets.txt"));
        out.write(summary);
        out.close();

        disconnectFromWarehouse();
    }

    private void printDataSetTable(ResultSet results) throws SQLException {
        print("Name", 35);
        print("WID", 4);
        print("Version", 12);
        print("LoadDate", 21);
        print("ChangeDate", 21);
        print("LoadedBy", 10);
        println();
        while (results.next()) {
            print(results.getString("Name"), 35);
            print(results.getString("WID"), 4);
            print(results.getString("Version"), 12);
            print(results.getString("LoadDate"), 21);
            print(results.getString("ChangeDate"), 21);
            print(results.getString("LoadedBy"), 10);
            println();
        }
    }

    private void println() {
        summary += "\n";
    }

    private void print(String text, int size) {
        String formatted = String.format("%-" + size + "s", text);
//        System.out.print(formatted);
        summary += formatted;
    }

    private Vector<TableDocument.Table> getTablesOfType(SchemaDocument.Schema schema, String tableType) {
        Vector<TableDocument.Table> list = new Vector<TableDocument.Table>();
        for (TableDocument.Table table : schema.getTableArray()) {
            if (table.getType().equals(tableType)) {
                list.add(table);
            }
        }
        return list;
    }

    public class DataSetStruct {
        public String wid;
        public String name;
        public String loadDate;
        public String changeDate;
        public String loadedBy;
        public String version;

        public DataSetStruct(String name, String wid, String version, String changeDate, String loadDate, String loadedBy) {
            this.changeDate = changeDate;
            this.loadDate = loadDate;
            this.loadedBy = loadedBy;
            this.name = name;
            this.version = version;
            this.wid = wid;
        }
    }
}
