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

import org.apache.log4j.Logger;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.ResultSetMetaData;

/**
 * @author Valerie Wagner
 *         Date: Sep 19, 2006
 */
public class DatabaseUtils {
    /**
     * Pretty-prints current row of a JDBC result set to standard out
     * Does not change position of ResultSet pointer
     * @param results ResultSet to be printed
     */
    public static void printResultSetRow(ResultSet results) {

        try {
            ResultSetMetaData md = results.getMetaData();
            System.out.println(md.getTableName(1));
            for (int i = 1; i <= md.getColumnCount(); i++) {
                System.out.println("  " +md.getColumnLabel(i) + "\t= " + results.getString(i));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }
    /**
     * Pretty-prints current row of a JDBC result set to the log
     * Does not change position of ResultSet pointer
     * @param results ResultSet to be printed
     */
    public static void logResultSetRow(Logger log, ResultSet results) {

        try {
            ResultSetMetaData md = results.getMetaData();
            log.info(md.getTableName(1));
            for (int i = 1; i <= md.getColumnCount(); i++) {
                log.info("  " +md.getColumnLabel(i) + "\t= " + results.getString(i));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }
}
