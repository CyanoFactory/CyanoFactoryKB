/* $Id: OracleWarehouse.java,v 1.5 2006/11/07 17:51:17 valerie Exp $ */
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

import com.sri.bw.OracleWIDGenerator;
import com.sri.bw.ReservedWIDGenerator;
import com.sri.bw.WIDGenerator;
import org.apache.log4j.Logger;

/**
 * Implements Warehouse functionality specific to Oracle (driver name and JDBC connection string)
 *
 * @author Valerie Wagner
 *         Date: Jun 14, 2004
 */
public class OracleWarehouse extends OracleDatabase implements Warehouse {
    private static final Logger log = Logger.getLogger(OracleWarehouse.class);

    private final OracleWIDGenerator WID_GENERATOR = new OracleWIDGenerator("WID_sequence", this);
    private final OracleWIDGenerator SPECIAL_WID_GENERATOR = new OracleWIDGenerator("SpecialWID_sequence", this);
    private ReservedWIDGenerator reservedWIDGenerator;
    private WIDGenerator regularWIDGenerator = WID_GENERATOR;
    private WIDGenerator specialWIDGenerator = SPECIAL_WID_GENERATOR;


    OracleWarehouse(String databaseHost, String databaseName, String databasePort)
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

    public WIDGenerator getSpecialWIDGenerator() {
        return specialWIDGenerator;
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

}
