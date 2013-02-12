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

import org.apache.log4j.Logger;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Set;
import com.sri.biospice.warehouse.util.LoaderMain;

/**
 * @author Valerie Wagner
 *         Date: Aug 1, 2006
 */
public class LoaderStatistics {

    private static final Logger log = Logger.getLogger(LoaderStatistics.class);
    private static Hashtable<GenericTable, Long> insertCounts = new Hashtable<GenericTable, Long>();

    private static int warningCount = 0;
    private static int errorCount = 0;
    private static int maxErrors = 0;
    private static boolean exitOnMaxErrors = false;

    public static void insertPerformed(GenericTable table) {
        if (insertCounts.containsKey(table)) {
            Long count = insertCounts.get(table) + 1;
            insertCounts.put(table, count);
        } else {
            insertCounts.put(table, new Long(1));
        }
    }

    public static void logStatistics() {
        int totalCount = 0;
        log.info("Insert counts: ");
        Set<GenericTable> keys = insertCounts.keySet();
        GenericTable [] tables = new GenericTable[keys.size()];
        tables = keys.toArray(tables);
        Arrays.sort(tables);
        for (GenericTable table : tables) {
            String countFormat = String.format("   %-35s: %d", table.getTableName(), insertCounts.get(table));
            log.info(countFormat);
            totalCount += insertCounts.get(table);
        }
        log.info("Total inserts: " + totalCount);
        log.info(warningCount + " warnings");
        log.info(errorCount + " errors");
    }

    public static void warningOccurred() {
        warningCount++;
    }

    public static void errorOccurred() {
        errorCount++;
        if (exitOnMaxErrors) {
            if (errorCount > maxErrors) {
                log.fatal("Maximum Number of errors allowed exceeded ("+maxErrors+").");
                throw new MaxErrorsRuntimeException("Maximum Number of errors allowed has been exceeded.  Program is exiting.");
            }
        }
    }

    public static int getErrorCount() {
        return errorCount;
    }

    public static int getWarningCount() {
        return warningCount;
    }

    /**
     * Add a check to the errorOccurred() method to see whether
     * we have exceeded the maximum number of errors allowed.
     * If so, errorOccurred() will throw a MaxErrorsRuntimeException
     * @param numErrors
     */
    public static void setMaxErrorsAllowed(int numErrors) {
        exitOnMaxErrors = true;
        maxErrors = numErrors;
    }
}
