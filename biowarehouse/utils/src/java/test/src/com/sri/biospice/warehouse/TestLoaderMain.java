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

import com.sri.biospice.warehouse.util.LoaderMain;
import java.util.Date;

/**
 * @author Valerie Wagner
 *         Date: May 9, 2006
 */
public class TestLoaderMain extends LoaderMain {

    private float earliestSupportedDataVersion;
    private float latestSupportedDataVersion;


    public TestLoaderMain(String [] args) {
        super(args, "test", "test");
    }

    public TestLoaderMain(String[] args, String loaderName,  String scriptName, boolean multipleFiles) {
        super(args, loaderName, scriptName, multipleFiles);
    }

    public TestLoaderMain(String loaderName,  boolean multipleFiles) {
        super(loaderName, multipleFiles);
    }


    public String getLoaderVersionNumber() {
        return null;
    }

    public String getLoaderBuildNumber() {
        return null;
    }



    public float getEarliestSupportedDataVersion() {
        return earliestSupportedDataVersion;
    }

    public void setEarliestSupportedDataVersion(float earliestSupportedDataVersion) {
        this.earliestSupportedDataVersion = earliestSupportedDataVersion;
    }

    public float getLatestSupportedDataVersion() {
        return latestSupportedDataVersion;
    }

    public void setLatestSupportedDataVersion(float latestSupportedDataVersion) {
        this.latestSupportedDataVersion = latestSupportedDataVersion;
    }


    public void setEarliestRelease(Date earliestSupportedDataRelease) {
        this.setEarliestSupportedDataRelease(earliestSupportedDataRelease);
    }


    public void setLatestRelease(Date latestSupportedDataRelease) {
        this.setLatestSupportedDataRelease(latestSupportedDataRelease);
    }

}
