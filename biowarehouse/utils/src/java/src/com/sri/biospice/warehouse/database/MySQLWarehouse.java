/* $Id: MySQLWarehouse.java,v 1.6 2006/11/07 17:51:17 valerie Exp $ */
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
package com.sri.biospice.warehouse.database;

import com.sri.bw.MySQLWIDGenerator;
import com.sri.bw.ReservedWIDGenerator;
import com.sri.bw.WIDGenerator;
import org.apache.log4j.Logger;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Implements Warehouse functionality specific to MySQL (driver name and JDBC connection string)
 *
 * @author Valerie Wagner
 *         Date: Jun 14, 2004
 */
public class MySQLWarehouse extends MySQLDatabase implements Warehouse {
    private static final Logger log = Logger.getLogger(MySQLWarehouse.class);

    private final MySQLWIDGenerator WID_GENERATOR = new MySQLWIDGenerator("WIDTable", this);
    private final MySQLWIDGenerator SPECIAL_WID_GENERATOR = new MySQLWIDGenerator("SpecialWIDTable", this);
    private ReservedWIDGenerator reservedWIDGenerator;
    private WIDGenerator regularWIDGenerator = WID_GENERATOR;
    private WIDGenerator specialWIDGenerator = SPECIAL_WID_GENERATOR;


    MySQLWarehouse(String databaseHost, String databaseName, String databasePort)
            throws ClassNotFoundException {
        super(databaseHost, databasePort, databaseName);
    }


    /**
     * Returns the next WID in the sequence.
     *
     * @return next WID on the sequence or zero if an error occurs.
     */
    public long getNextWID() {
        return regularWIDGenerator.generate();
    }


    /**
     * Returns the next SpecialWID in the sequence.
     * If specialWID exceeds a maximum, return a normal WID instead.
     *
     * @return next SpecialWID on the sequence or zero if an error occurs.
     */
    public long getNextSpecialWID() {
        return specialWIDGenerator.generate();
    }


    /**
     * Do not call this method until the Warehouse connection is established,
     * as the reserved WID range is determined by querying the Warehouse
     * Generates
     *
     * @return The next reserved WID in the sequence
     */
    public long getNextReservedWID() {
        return getReservedWIDGenerator().generate();
    }

    public WIDGenerator getWIDGenerator() {
        return regularWIDGenerator;
    }


    /**
     * Do not call this method until the Warehouse connection is established,
     * as the reserved WID range is determined by querying the Warehouse
     * Generates
     *
     * @return The reserved WID Generator
     */
    public WIDGenerator getReservedWIDGenerator() {
        if (reservedWIDGenerator == null) {
            reservedWIDGenerator = new ReservedWIDGenerator(this);
        }
        return reservedWIDGenerator;
    }

    public void setWIDGenerator(WIDGenerator generator) {
        this.regularWIDGenerator = generator;
    }

    public void setSpecialWIDGenerator(WIDGenerator generator) {
        this.specialWIDGenerator = generator;
    }

    public WIDGenerator getSpecialWIDGenerator() {
        return specialWIDGenerator;
    }

    private void resetWIDTables() throws TableResetException {
        ResultSet results = null;
        long wid = -1;
        try {
            // Retrieve last AUTO_INCREMENT value
            String sqlLine = "select max(WID) FROM DataSet";
            results = executeQuery(sqlLine);

            if (results.next()) {
                wid = results.getLong(1);
            }
            results.close();

            if (wid <= 1) {
                log.fatal("Reset of WID Tables failed; WIDTable and SpecialWIDTable must be repaired manually");
                throw new TableResetException("Reset of WID Tables failed (max specialWID <= 1); WIDTable and SpecialWIDTable must be repaired manually");
            }

            // Set SpecialWID counter to to the largest WID + 2
            sqlLine = "ALTER TABLE SpecialWIDTable AUTO_INCREMENT=" + (wid + 2);
            executeUpdate(sqlLine);
            log.info("SpecialWIDTable counter set to " + (wid + 2));

            // Determine the likely legal value for next (regular) WID, and
            // insert into WIDTable
            wid = 0;
            sqlLine = "SELECT max(OtherWID) FROM Entry";
            results = executeQuery(sqlLine);
            if (results.next()) {
                wid = results.getLong(1);
            }

            results.close();

            if (wid <= 1) {
                log.fatal("Reset of WID Tables failed; WIDTable and SpecialWIDTable must be repaired manually");
                throw new TableResetException("Reset of WID Tables failed (max WID <= 1); WIDTable and SpecialWIDTable must be repaired manually");
            }

            // Set WID counter to to the largest WID * 2 (conservative value)
            sqlLine = "ALTER TABLE WIDTable AUTO_INCREMENT=" + (wid * 2);
            executeUpdate(sqlLine);
            log.info("WIDTable counter set to " + (wid * 2));
        } catch (SQLException ex) {
            ex.printStackTrace();
            System.exit(1);
        } finally {
            if (results != null) {
                try {
                    results.close();
                } catch (SQLException ex) {
                    ex.printStackTrace();
                }
            }
        }

    }
}
