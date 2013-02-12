/* $Id: Warehouse.java,v 1.5 2006/11/07 17:51:17 valerie Exp $ */
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

import com.sri.bw.WIDGenerator;

// todo: Make the ability to check that this version of the Warehouse code base
// is the same as the actual database being run against
// e.g. If we run 3.5 code against database containing 3.0 schema, our method can check that

/**
 * Provides database interaction specific to the BioSPICE Warehouse
 *
 * @author Valerie Wagner
 *         Date: Jun 14, 2004
 */
public interface Warehouse extends Database {
    // todo: Warehouse instances should probably check whether they are
    // compatible with the version of the actual database instance

    /**
     * The maximum allowed "special WID" size
     */
    public static final int MAX_SPECIAL_WID = 999;


    /**
     * @return The next WID available for Object entries
     */
    public long getNextWID();


    /**
     * @return The next "Special" WID available for identifiying DataSets
     */
    public long getNextSpecialWID();

    /**
     * Do not call this method until the Warehouse connection is established,
     * as the reserved WID range is determined by querying the Warehouse
     * Generates
     *
     * @return The reserved WID in the sequence
     */
    public long getNextReservedWID();

    WIDGenerator getWIDGenerator();

    WIDGenerator getSpecialWIDGenerator();

    /**
     * Do not call this method until the Warehouse connection is established,
     * as the reserved WID range is determined by querying the Warehouse
     * Generates
     *
     * @return The reserved WID Generator
     */
    WIDGenerator getReservedWIDGenerator();

    public void setWIDGenerator(WIDGenerator generator);
    
    public void setSpecialWIDGenerator(WIDGenerator generator);
}
