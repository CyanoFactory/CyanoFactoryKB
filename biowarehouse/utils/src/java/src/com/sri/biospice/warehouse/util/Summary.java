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
import com.sri.biospice.warehouse.util.LoaderMain;
import com.sri.bw.dbschema.SchemaDocument;
import com.sri.bw.dbschema.TableDocument;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.util.Arrays;
import java.util.Hashtable;
import org.apache.log4j.Logger;

/**
 * @author Valerie Wagner
 *         Date: May 12, 2006
 */
public class Summary extends LoaderMain {

    private static final Logger log = Logger.getLogger(Summary.class);

    private boolean showEmpty = false;

    public static void main(String[] args) throws Exception {
        new Summary(args);

    }


    public Summary(String [] args) throws Exception {
        super(args, null, null, false);
        if (!properties.validate()) {
            return;
        }

        File outFile = new File("summary.html");
        BufferedWriter buf = new BufferedWriter(new FileWriter(outFile));


        parseCommandLine(args, null);

        connectToWarehouse();
        Warehouse warehouse = WarehouseManager.getWarehouse();
        ResultSet results;
        int totalRows = 0;
        Hashtable<String, TableDocument.Table> tableList = new Hashtable<String, TableDocument.Table>();

//        buf.write("table.column, #rows populated, max size seen, max allowed, type\n");
        try {
            buf.write("<html>");
            buf.write("<head>");
            buf.write("<title>Summary</title>");
            buf.write("</head>");
            buf.write("<body>");
            buf.write("<h1>Summary</h1>\n");
            buf.write("<h3>Table Index</h3>\n");
            String filename = properties.getFile();

            log.info("Using schema file: " + filename);
            SchemaDocument schemaDoc = SchemaDocument.Factory.parse(new File(filename));
            SchemaDocument.Schema schema = schemaDoc.getSchema();


            TableDocument.Table [] tables = schema.getTableArray();
            String [] tableNames = new String[tables.length];
            for (int i = 0; i < tableNames.length; i++) {
                tableNames[i] = tables[i].getName();
                tableList.put(tableNames[i], tables[i]);
            }
            Arrays.sort(tableNames);
                for (String name : tableNames) {
                    buf.write("<a href=\"#" + name + "\">" + name + "</a><br>\n");
                }
            // for each table
//            for (TableDocument.Table table : tables) {
                  for(int j = 0; j<tableNames.length; j++) {
                      String tableName = tableNames[j];
                      if (tableName.equals("Entry")) {
                          continue;
                      }
                      String query = "select * from " + tableName;
                      results = warehouse.executeQuery(query);
                      if (!results.next()) {
                          if (!showEmpty) {
                              continue;
                          }
                      }
                      results.previous();
                      buf.write("<h2><a name=\"" + tableName + "\">" + tableName + "<a></h2>\n");
                      buf.write("<a href=\"doc/mage-schema.html#" + tableName + "\">Table Documentation</a>");
                      buf.write("<table border=\"1\">");

                      // table header
                      buf.write("<tr>");
                      ResultSetMetaData md = results.getMetaData();
                      int numCols = md.getColumnCount();
                      for (int i = 0; i < numCols; i++) {
                          buf.write("<th>" + md.getColumnName(i + 1) + "</th>");
                      }
                      buf.write("</tr>");

                      boolean link;
                      int rowCount = 0;
                      String linkRef;
                      // for each row in the table
                      while (results.next()) {
                          rowCount++;
                          totalRows++;
                          buf.write("<tr>");
                          for (int i = 0; i < numCols; i++) {
                              Object cellObj = results.getObject(i + 1);
                              String text = cellObj == null ? "<br/>" : cellObj.toString();
                              link = false;
                              linkRef = "";

                              // If a WID column create an achor
                              String columnName = md.getColumnName(i + 1);
                              if (columnName.equals("WID")) {
                                  buf.write("<a name=\"WID" + text + "\"/>");
                              } else if (tableList.get(tableName).getColumnArray(i).getType().equals("@wid")) {
                                  link = true;
                                  linkRef = "#WID" + text;
                              } else if (columnName.equals("MAGEClass")) {
                                  link = true;
                                  linkRef = "doc/hibernate-map.hbm.html#" + text;
                              }
                              buf.write("<td>" +
                                      (link ? "<a href=\"" + linkRef + "\">" : "")
                                      + text +
                                      (link ? "</a>" : "") +
                                      "</td>\t");
                          }
                          buf.write("\n");
                          buf.write("</tr>");
                      }
                      results.close();
                      buf.write("</table>");
                      buf.write("Rows: " + rowCount + "<br/>\n");
                  }
        } catch (Exception e) {
            e.printStackTrace();
        }
        buf.write("<br/>Total Rows: " + totalRows + "<br/>");
        buf.write("</body>");
        buf.write("</html>");

        disconnectFromWarehouse();
        buf.close();
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
