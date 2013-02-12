/* $Id: MySQLWIDGenerator.java,v 1.1 2006/11/02 20:44:40 valerie Exp $ */
package com.sri.bw;

import org.apache.log4j.Logger;
import com.sri.biospice.warehouse.database.Warehouse;
import com.sri.biospice.warehouse.database.WarehouseManager;
import com.sri.biospice.warehouse.util.LoaderMain;
import java.sql.SQLException;
import java.sql.ResultSet;

/**
 * @author Valerie Wagner
*         Date: Nov 2, 2006
*/
public class MySQLWIDGenerator implements WIDGenerator {
    private static final Logger log = Logger.getLogger(MySQLWIDGenerator.class);
    private final String deleteSQL;
    private final String insertSQL;
    private final Warehouse warehouse;

    public MySQLWIDGenerator(String widTableName, Warehouse warehouse) {
        deleteSQL = "delete from " + widTableName;
        insertSQL = "insert into " + widTableName + " () values ()";
        this.warehouse = warehouse;
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
    public long getNextWID() throws SQLException {

        // todo consider whether there will be any transaction or concurrency problems here
        ResultSet results = null;
        long nextWID;

        try {
            // todo: can we do anything to make this more efficient, does mysql use prepared statements?

            // If DB is MySQL, use AUTO_INCREMENT in WIDTable

            // Delete rows from WIDTable to prevent it from growing too large
            warehouse.executeUpdate(deleteSQL);

            // Add row to WIDTable to advance AUTO_INCREMENT
            int result = warehouse.executeUpdate(insertSQL);

            // Retrieve last AUTO_INCREMENT value
            String sqlLine = "select LAST_INSERT_ID()";
            results = warehouse.executeQuery(sqlLine);
            if (results.next()) {
                // Retrieve column values for this row
                nextWID = results.getLong(1);
                results.close();
                return nextWID;
            }

        }
        catch (SQLException error) {
            log.error(error);
            try {
                if (results != null) {
                    results.close();
                }
            }
            catch (SQLException error2) {
            }
            throw error;

        }
        return 0;
    }
}
