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

import com.sri.biospice.warehouse.util.LoaderMain.DataSupport;
import junit.framework.TestCase;
import org.apache.log4j.BasicConfigurator;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;

/**
 * @author Valerie Wagner
 *         Date: May 9, 2006
 */
public class LoaderMain_Test extends TestCase {


    public LoaderMain_Test(String testName) {
        super(testName);
    }

    protected void setUp() throws Exception {
        BasicConfigurator.resetConfiguration();
        BasicConfigurator.configure();
    }

    public void testCheckDataReleaseSupported() {

        // GregorianCalendar months are zero-based
        Calendar cal1 = new GregorianCalendar(2004, 0, 1);
        Date date1 = cal1.getTime();
        Calendar cal2 = new GregorianCalendar(2005, 6, 7);
        Date date2 = cal2.getTime();
        Calendar cal3 = new GregorianCalendar(2008, 1, 28);
        Date date3 = cal3.getTime();
        Calendar cal4 = new GregorianCalendar(2009, 10, 4);
        Date date4 = cal4.getTime();

        checkDataRelease(date1, date1, "2004/01/01", DataSupport.SUPPORTED);
        checkDataRelease(date1, date1, "2004/1/1", DataSupport.SUPPORTED);
        checkDataRelease(date1, date1, "2004/1/01", DataSupport.SUPPORTED);
        checkDataRelease(date1, date1, "2005/12/01", DataSupport.UNKNOWN);
        checkDataRelease(date1, date1, "2003/1/01", DataSupport.NOT_SUPPORTED);

        checkDataRelease(date1, date2, "2004/01/01", DataSupport.SUPPORTED);
        checkDataRelease(date1, date2, "2004/1/1", DataSupport.SUPPORTED);
        checkDataRelease(date1, date2, "2004/1/01", DataSupport.SUPPORTED);
        checkDataRelease(date1, date2, "2005/12/01", DataSupport.UNKNOWN);
        checkDataRelease(date1, date2, "2003/1/01", DataSupport.NOT_SUPPORTED);

        checkDataRelease(date3, date4, "2008/2/28", DataSupport.SUPPORTED);
        checkDataRelease(date3, date4, "2009/11/4", DataSupport.SUPPORTED);
        checkDataRelease(date3, date4, "2009/3/03", DataSupport.SUPPORTED);
        checkDataRelease(date3, date4, "2009/11/5", DataSupport.UNKNOWN);
        checkDataRelease(date3, date4, "2003/1/01", DataSupport.NOT_SUPPORTED);

        checkDataRelease(date3, date4, "2008/3", DataSupport.SUPPORTED);
        checkDataRelease(date3, date4, "2009/11", DataSupport.SUPPORTED);
        checkDataRelease(date3, date4, "2009/05", DataSupport.SUPPORTED);
        checkDataRelease(date3, date4, "2009/12", DataSupport.UNKNOWN);
        checkDataRelease(date3, date4, "2003/1", DataSupport.NOT_SUPPORTED);

        checkDataRelease(date1, date1, "2004/01", DataSupport.SUPPORTED);
        checkDataRelease(date1, date1, "2004/1", DataSupport.SUPPORTED);
        checkDataRelease(date1, date1, "2004/2", DataSupport.UNKNOWN);
        checkDataRelease(date1, date1, "2003/12", DataSupport.NOT_SUPPORTED);

    }

    private void checkDataRelease(Date earliest, Date latest, String date, DataSupport supported) {
        TestLoaderMain main = new TestLoaderMain(new String []{"-r", date});
        main.setEarliestRelease(earliest);
        main.setLatestRelease(latest);
        assertEquals("Data release support wrong [" +
                earliest + "-" + latest + "] " + date,
                supported, main.checkDataReleaseSupported());
    }

    public void testCheckDataVersionSupported() {
        checkDataVersion(5, 5, "5", DataSupport.SUPPORTED);
        checkDataVersion(5.0F, 5.0F, "5", DataSupport.SUPPORTED);
        checkDataVersion(1.1F, 10.2F, "5", DataSupport.SUPPORTED);
        checkDataVersion(5, 5, "1", DataSupport.NOT_SUPPORTED);
        checkDataVersion(5.0F, 5.0F, "1", DataSupport.NOT_SUPPORTED);
        checkDataVersion(1.1F, 10.2F, "1", DataSupport.NOT_SUPPORTED);
        checkDataVersion(5, 5, "20.2", DataSupport.UNKNOWN);
        checkDataVersion(5.0F, 5.0F, "20.2", DataSupport.UNKNOWN);
        checkDataVersion(1.1F, 10.2F, "20.2", DataSupport.UNKNOWN);
//        checkDataVersion(5, 5, "", DataSupport.NOT_SUPPORTED);
//        checkDataVersion(5.0F, 5.0F, "", DataSupport.NOT_SUPPORTED);
//        checkDataVersion(1.1F, 10.2F, "", DataSupport.NOT_SUPPORTED);
//        checkDataVersion();
//        checkDataVersion();
    }

    private void checkDataVersion(float earliest, float latest, String version, DataSupport supported) {
        TestLoaderMain main = new TestLoaderMain(new String []{"-v", version});
        main.setEarliestSupportedDataVersion(earliest);
        main.setLatestSupportedDataVersion(latest);
        assertEquals("Data version not supported", supported, main.checkDataVersionSupported());
    }

}