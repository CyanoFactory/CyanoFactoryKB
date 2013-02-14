/* $Id: WarehouseFunctionalTest.java,v 1.4 2006/08/10 19:43:31 valerie Exp $ */
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
package com.sri.biospice.warehouse;

import com.sri.biospice.warehouse.database.Warehouse;
import com.sri.biospice.warehouse.database.WarehouseManager;
import com.sri.biospice.warehouse.util.WarehouseTestUtils;
import junit.framework.TestCase;
import org.apache.log4j.BasicConfigurator;
import java.sql.SQLException;
import java.sql.DatabaseMetaData;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Sep 27, 2004
 */
public class WarehouseFunctionalTest extends TestCase {
    protected Warehouse warehouse;

    public WarehouseFunctionalTest(String testName) throws SQLException {
        super(testName);
        WarehouseTestUtils.initWarehouseTest();
        warehouse = WarehouseManager.getWarehouse();
    }

    protected void setUp() throws Exception {
        BasicConfigurator.resetConfiguration();
        BasicConfigurator.configure();

        System.out.println("*******************************************************************");
        System.out.println("Test properties:");
        System.out.println("databaseDBMSType = " + WarehouseTestUtils.databaseDBMSType);
        System.out.println("databaseHost = " + WarehouseTestUtils.databaseHost);
        System.out.println("databasePort = " + WarehouseTestUtils.databasePort);
        System.out.println("databaseName = " + WarehouseTestUtils.databaseName);
        System.out.println("username = " + WarehouseTestUtils.username);
        DatabaseMetaData meta = warehouse.getMetaData();
        System.out.println("URL = " + meta.getURL());
        System.out.println("RDBMS = " + meta.getDatabaseProductName() + ", version: " + meta.getDatabaseProductVersion() );
        System.out.println("JDBC driver = " + meta.getDriverName() + ", version: " + meta.getDriverVersion() );
        System.out.println("*******************************************************************\n");
    }


    protected void commit() {
        WarehouseManager.getWarehouse().commit();
    }

}
