/* $Id: ReservedWIDGenerator.java,v 1.1 2006/11/07 00:11:12 valerie Exp $ */
package com.sri.bw;

import com.sri.biospice.warehouse.database.Warehouse;
import com.sri.biospice.warehouse.util.LoaderMain;
import org.apache.log4j.Logger;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Generates reserved WIDs in the range Warehouse.MaxSpecialWID+1 to Warehouse.MaxReservedWID
 *
 * @author Valerie Wagner
 *         Date: Nov 6, 2006
 */
public class ReservedWIDGenerator implements WIDGenerator {

    private static final Logger log = Logger.getLogger(ReservedWIDGenerator.class);

    private boolean initialized = false;
    private Long reservedWIDCounter;
    private long minReservedWID;
    private long maxReservedWID;


    public ReservedWIDGenerator(Warehouse warehouse) {

        String query = "select MaxSpecialWID from Warehouse";
        try {
            ResultSet results = warehouse.executeQuery(query);
            if (results.next()) {
                minReservedWID = results.getLong(1) + 1;
            } else {
                LoaderMain.fail("Error determining reserved WID sequence. Warehouse row not populated.  Cannot continue.");
            }
            results.close();

            query = "select MaxReservedWID from Warehouse";
            results = warehouse.executeQuery(query);
            if (results.next()) {
                maxReservedWID = results.getLong(1);
            } else {
                LoaderMain.fail("Error determining reserved WID sequence. Warehouse row not populated.  Cannot continue.");
            }
            results.close();

        } catch (SQLException e) {
            log.fatal("Error querying for WID state: " + query, e);
            LoaderMain.fail("Error determining reserved WID sequence.  Cannot continue.");
        }


        reservedWIDCounter = minReservedWID;
    }

    public Long generate() {
        if (reservedWIDCounter <= maxReservedWID) {
            reservedWIDCounter++;
            return reservedWIDCounter;
        } else {
            LoaderMain.fail("Reserved WID counter has exceeded maximum allowed size of " + maxReservedWID);
            return new Long(0);
        }
    }

}
