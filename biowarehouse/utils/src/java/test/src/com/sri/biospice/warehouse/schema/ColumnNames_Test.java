/* $Id: ColumnNames_Test.java,v 1.3 2006/11/06 21:51:36 valerie Exp $ */
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
package com.sri.biospice.warehouse.schema;

import com.sri.biospice.warehouse.WarehouseFunctionalTest;
import com.sri.biospice.warehouse.schema.object.Reaction;
import com.sri.biospice.warehouse.schema.object.ObjectTableBase;
import com.sri.biospice.warehouse.database.table.Column;
import com.sri.biospice.warehouse.database.table.TableMetaData;
import com.sri.bw.dbschema.SchemaDocument;
import com.sri.bw.dbschema.TableDocument;
import com.sri.bw.ObjectTable;
import org.apache.log4j.Logger;
import org.apache.xmlbeans.XmlException;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.sql.SQLException;
import java.util.Vector;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Aug 12, 2004
 */
public class ColumnNames_Test extends WarehouseFunctionalTest {
    private static final Logger log = Logger.getLogger(ColumnNames_Test.class);

    public ColumnNames_Test(String testName) throws SQLException {
        super(testName);
    }

    public void atestAllTables() {
        int[] tableTypes = TableFactory.getTabletypes();
        for (int i = 0; i < tableTypes.length; i++) {
            int tableType = tableTypes[i];
            Table table = TableFactory.createTestTable(tableType);
            assertNotNull("TableFactory return null table for type " + tableType, table);
            checkTableColumnNames(table);
        }
    }

    private void checkTableColumnNames(Table table) {
        log.info("Checking table " + table.getTableName());
        TableMetaData tableMetaData = table.getTableMetaData();
        Vector columns = table.getColumns();

        for (int i = 0; i < columns.size(); i++) {
            Column column = (Column)columns.elementAt(i);
            assertEquals("Invalid column name: " + table.getTableName() + "." + column.getColumnName(),
                    true, tableMetaData.hasColumn(column.getColumnName()));
        }

        assertEquals("Wrong number of columns for table " + table.getTableName(), tableMetaData.getNumberOfColumns(), columns.size());
    }

    public void testSchemaMatch() throws IOException, XmlException, IllegalAccessException, InvocationTargetException {
        SchemaDocument.Schema schema = loadSchema();
        for (TableDocument.Table table : schema.getTableArray()) {
            if (table.getType().equals("object")) {
                String tableName = table.getName();
                String className = "com.sri.biospice.warehouse.schema.object." + tableName;
                Class clazz;
                try {
                    if ((clazz = Class.forName(className)) != null) {
                        Constructor[] constructors = clazz.getConstructors();
                        if (constructors.length == 0) {
                            fail("No constructors for this class: " + tableName);
                        } else {
                            for (Constructor constr : constructors) {
                                if (constr.getParameterTypes().length == 1) {
                                    long datasetWID = 0;
                                    log.info("Constructing " + tableName );
                                    ObjectTableBase objTable = (ObjectTableBase)constr.newInstance(datasetWID);
                                    checkTableColumnNames(objTable);
                                    break;
                                }
                            }
                        }
                    }
                } catch (NoClassDefFoundError e) {
                    log.info("No Java class defined for " + tableName);
                }catch (ClassNotFoundException e) {
                    log.info("No Java class defined for " + tableName);
                } catch (InstantiationException e) {
                    fail("Didn't properly instantiate this class: " + tableName);
                }
            }
        }

    }

    private SchemaDocument.Schema loadSchema() throws IOException, XmlException {
        File file = new File("../../../../schema/all-schema.xml");
        return SchemaDocument.Factory.parse(file).getSchema();
    }
}
