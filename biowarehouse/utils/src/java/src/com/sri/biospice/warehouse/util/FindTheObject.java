/* $Id: FindTheObject.java,v 1.3 2008/10/10 22:02:20 valerie Exp $ */
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
import java.io.File;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;


/**
 * @author Valerie.Wagner@sri.com
 *         Date: Sep 13, 2004
 */
public class FindTheObject extends LoaderMain {
    private static final Logger log = Logger.getLogger(FindTheObject.class);


    public FindTheObject(String[] args) {
        super(args, "FindTheObject", "FindTheObject");
        connectToWarehouse();
    }


    public static void main(String[] args) throws SQLException {
        FindTheObject find = new FindTheObject(args);
        find.findIt();
    }


    public void findIt() throws SQLException {
        String wid = properties.getProperty("wid");
        String tableName, query;
        boolean objectFound = false;
        ResultSet results;

        String schemaFile = properties.getFile();
        File file = new File(schemaFile);
        SchemaDocument.Schema schema = null;
        try {
            schema = SchemaDocument.Factory.parse(file).getSchema();
        } catch (XmlException e) {
            log.fatal("Error parsing schema", e);
            fail("Cannot read warehouse schema");
        } catch (IOException e) {
            log.fatal("Error parsing schema", e);
            fail("Cannot read warehouse schema");
        }


        Warehouse warehouse = WarehouseManager.getWarehouse();

        for (int i = 0; i < schema.getTableArray().length && !objectFound; i++) {
            TableDocument.Table table = schema.getTableArray(i);
            if (table.getType().equals("object")) {
                tableName = table.getName();
                query = "Select * from " + tableName + " where WID=" + wid;
                log.info(query);
                results = warehouse.executeQuery(query);
                if (results.next()) {
                    objectFound = true;
                    log.info("Match in the " + tableName + " table");
                }
            }
        }

        if (!objectFound) {
            log.info("No match found");
        }

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
