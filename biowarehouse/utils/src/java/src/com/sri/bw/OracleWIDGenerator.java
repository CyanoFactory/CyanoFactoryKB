/* $Id: OracleWIDGenerator.java,v 1.2 2006/11/07 00:11:12 valerie Exp $ */
package com.sri.bw;

import com.sri.biospice.warehouse.database.Warehouse;
import com.sri.biospice.warehouse.util.LoaderMain;
import org.apache.log4j.Logger;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author Valerie Wagner
 *         Date: Nov 2, 2006
 */
public class OracleWIDGenerator implements WIDGenerator {

    private static final Logger log = Logger.getLogger(OracleWIDGenerator.class);
    private final Warehouse warehouse;
    private final String selectSequence;


    public OracleWIDGenerator(String sequenceName, Warehouse warehouse) {
        this.warehouse = warehouse;
        selectSequence = "select " + sequenceName + ".NextVal from dual";
    }

    public Long generate() {
        long wid = 0;
        try {
            wid = getNextWID();
        } catch (SQLException e) {
            log.fatal("Error generating WID", e);
            LoaderMain.fail("An error occurred during WID generation.  This is a serious error and the application cannot continue.");
        }
        return new Long(wid);
    }


    /**
     * Returns the next WID in the sequence.
     *
     * @return next WID on the sequence or zero if an error occurs.
     */
    private long getNextWID() throws SQLException {
        ResultSet results = null;
        long nextWID;

        try {
            // If DB is Oracle, use SEQUENCE
            results = warehouse.executeQuery(selectSequence);
            if (results.next()) {
                // Retrieve column values for this row
                nextWID = results.getLong(1);
                results.close();
                return nextWID;
            }

        }
        catch (SQLException error) {
            log.fatal("getNextWID", error);
            log.fatal(selectSequence);
            try {                                        
                results.close();
            }
            catch (SQLException error2) {
                log.fatal("getNextWID2", error2);
            }
            throw error;
        }
        log.error("returning zero as WID");
        return 0;
    }
}
