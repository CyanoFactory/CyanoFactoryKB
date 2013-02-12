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

import com.sri.biospice.warehouse.WarehouseFunctionalTest;
import com.sri.biospice.warehouse.database.WarehouseManager;
import com.sri.biospice.warehouse.schema.DataSet;
import org.apache.log4j.Logger;
import java.sql.SQLException;
import java.util.Hashtable;
import java.util.Date;

/**
 * @author Valerie Wagner
 *         Date: Jul 26, 2006
 */
public class Insert_Test extends WarehouseFunctionalTest {

    private static final Logger log = Logger.getLogger(Insert_Test.class);

    private Hashtable<String, TableInsert> tableName2Insert;
    private DataSet dataset;

    public Insert_Test(String testName) throws SQLException {
        super(testName);
    }

    protected void setUp() {
        try {
            super.setUp();
        } catch (Exception e) {
            log.error(e);
            fail(e.getMessage());
        }
        log.info("Insert_Test.info()");
        dataset = new DataSet("Insert_Test");
        dataset.setReleaseDate(new Date());
        try {
            dataset.store();
        } catch (SQLException e) {
            log.error(e);
            fail(e.getMessage());
        }
    }

    public void testConstruction() throws Exception {
        log.info("Insert_Test.testConstruction()");

        WarehouseManager.getWarehouse().commit();

        // todo: implement
//        TableConstructor tc = new TableConstructor(dataset, " ../../../../schema/all-schema.xml");
//        TableInsert insertRow = tc.getInsertForTable("Contact");
//
//        insertRow.put("firstName", "George");
//        insertRow.store();
//        ObjectDescriptorTable dbidTable = table.addDescriptor("DBID");
//        dbidTable.put("XID", "myid");
//
//        table.store();
//        table.put("random", "value");
//        table.store();
    }

}