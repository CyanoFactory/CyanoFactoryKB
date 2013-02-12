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
package com.sri.biospice.warehouse.database.table;

import com.sri.biospice.warehouse.WarehouseFunctionalTest;
import com.sri.biospice.warehouse.database.WarehouseManager;
import com.sri.biospice.warehouse.schema.DataSet;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Calendar;
import java.util.Date;

/**
 * @author Valerie Wagner
 *         Date: Jun 8, 2005
 */
public class TimestampColumn_Test extends WarehouseFunctionalTest {

    public TimestampColumn_Test(String testName) throws SQLException {
        super(testName);
    }

    public void testUpdatingTimestamps() throws SQLException {

        // Store the dataset with a LoadTime of now
        DataSet dataset = new DataSet("test");
        Date loadDate = dataset.getLoadDate();
        dataset.store();
        System.out.println(dataset);
        System.out.println("Load date: " + loadDate + " [" + loadDate.getTime() + "]");

        // Query for LoadTime and make sure it is the same as we expected
        String query = "Select LoadDate from DataSet where WID=" + dataset.getWID();
        ResultSet results = WarehouseManager.getWarehouse().executeQuery(query);
        assertNotNull("Results are null", results);
        assertTrue("Results contained no data", results.next());
        Date storedDate = new Date(results.getTimestamp("LoadDate").getTime());
        System.out.println("Stored load date: " + storedDate + " [" + storedDate.getTime() + "]");
        checkEqualDateTime("LoadDate after first storage:", loadDate, storedDate);

        // Set ChangeDate to 20 seconds after LoadDate
        Date changeDate = new Date(loadDate.getTime() + 20000);
        System.out.println("Change date: " + changeDate + " [" + changeDate.getTime() + "]");
        dataset.setChangeDate(changeDate);
        dataset.update();

        // Query for LoadDate and ChangeDate
        query = "Select LoadDate, ChangeDate from DataSet where WID=" + dataset.getWID();
        results = WarehouseManager.getWarehouse().executeQuery(query);
        assertNotNull("Results are null", results);
        assertTrue("Results contained no data", results.next());

        // Check that LoadDate stayed the same
        storedDate = new Date(results.getTimestamp("LoadDate").getTime());
        System.out.println("Stored load date after update: " + storedDate + " [" + storedDate.getTime() + "]");
        checkEqualDateTime("LoadDate after update: ", loadDate, storedDate);

        // Check that ChangeDate got set correctly
        storedDate = new Date(results.getTimestamp("ChangeDate").getTime());
        System.out.println("Stored change date after update: " + storedDate + " [" + storedDate.getTime() + "]");
        checkEqualDateTime("ChangeDate after update: ", changeDate, storedDate);
    }

    public void testReloadingObject() throws SQLException {

        // Store the dataset with a LoadTime of now
        DataSet dataset1 = new DataSet("test");
        Date loadDate1 = dataset1.getLoadDate();
        dataset1.store();
        System.out.println("\n" + dataset1);
        System.out.println("Load date: " + loadDate1 + " [" + loadDate1.getTime() + "]");

        // Create a new DataSet object and load it
        DataSet dataset2 = new DataSet(dataset1.getWID());
        dataset2.load();
        Date loadDate2 = dataset2.getLoadDate();
        System.out.println("Load date2: " + loadDate2 + " [" + loadDate2.getTime() + "]");
        checkEqualDateTime("Reloaded object: ", loadDate1, loadDate2);
    }

    private void checkEqualDateTime(String message, Date expectedDate, Date actualDate) {
        Calendar expectedCal = Calendar.getInstance();
        expectedCal.setTime(expectedDate);

        Calendar actualCal = Calendar.getInstance();
        actualCal.setTime(actualDate);

        assertEquals(message + "Year not equal", expectedCal.get(Calendar.YEAR), actualCal.get(Calendar.YEAR));
        assertEquals(message + "Month not equal", expectedCal.get(Calendar.MONTH), actualCal.get(Calendar.MONTH));
        assertEquals(message + "Day not equal", expectedCal.get(Calendar.DAY_OF_MONTH), actualCal.get(Calendar.DAY_OF_MONTH));
        assertEquals(message + "Hour not equal", expectedCal.get(Calendar.HOUR_OF_DAY), actualCal.get(Calendar.HOUR_OF_DAY));
        assertEquals(message + "Minute not equal", expectedCal.get(Calendar.MINUTE), actualCal.get(Calendar.MINUTE));
        assertEquals(message + "Second not equal", expectedCal.get(Calendar.SECOND), actualCal.get(Calendar.SECOND));
        System.out.println(message + "Expected millis: " + expectedCal.get(Calendar.MILLISECOND));
        System.out.println(message + "Actual millis: " + actualCal.get(Calendar.MILLISECOND));
    }
}