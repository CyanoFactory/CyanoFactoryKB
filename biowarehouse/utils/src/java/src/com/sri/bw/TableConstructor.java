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
package com.sri.bw;

import com.sri.biospice.warehouse.schema.DataSet;
import com.sri.bw.dbschema.SchemaDocument;
import com.sri.bw.dbschema.TableDocument;
import org.apache.log4j.Logger;
import org.apache.xmlbeans.XmlException;
import java.io.*;
import java.sql.SQLException;
import java.util.Hashtable;

/**
 * @author Valerie Wagner
 *         Date: Jul 27, 2006
 */
public class TableConstructor {

    private static final Logger log = Logger.getLogger(TableConstructor.class);

    private Hashtable<String, GenericTable> name2table = new Hashtable<String, GenericTable>();
    private DataSet dataset;

    public TableConstructor(DataSet dataset, String schemaFilePath) throws FileNotFoundException {
        this(dataset, new File(schemaFilePath));
    }

    public TableConstructor(DataSet dataset, File schemaFile) throws FileNotFoundException {
        this(dataset, new FileInputStream(schemaFile));
    }

    public TableConstructor(DataSet dataset, InputStream schemaFile) {
        // todo: write unit tests for this class against both oracle and mysql
        this.dataset = dataset;
        try {
            SchemaDocument.Schema schema = SchemaDocument.Factory.parse(schemaFile).getSchema();
            ObjectDescriptorTable entryTable = new ObjectDescriptorTable(dataset, "Entry");
            for (TableDocument.Table table : schema.getTableArray()) {
                if (table.getType().equals("object")) {
                    name2table.put(table.getName(), new ObjectTable(dataset, table.getName(), entryTable));
                } else if (table.getType().equals("object-descriptor")) {
                    name2table.put(table.getName(), new ObjectDescriptorTable(dataset, table.getName()));
                } else if (table.getType().equals("associative")) {
                    name2table.put(table.getName(), new GenericTable(dataset, table.getName()));
                } else  {
                    name2table.put(table.getName(), new GenericTable(dataset, table.getName()));
                }
            }
            // todo replace with table construction exception
        } catch (Exception e){
            log.fatal("Error creating table", e);
            LoaderStatistics.errorOccurred();
        }
    }

    private GenericTable getTable(String tableName) {
        return name2table.get(tableName);
    }

    public TableInsert getInsertForTable(String tableName) {
        if (name2table.containsKey(tableName)) {
            return getTable(tableName).newInsert();
        } else {
            log.error("No table for " + tableName);
            log.error("Table list: " + name2table.keySet());
            LoaderStatistics.errorOccurred();
        }
        return null;
    }

    public DataSet getDataset() {
        return dataset;
    }
}

