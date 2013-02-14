/* $Id: OracleDatabase.java,v 1.2 2007/08/29 23:44:02 nguo Exp $ */
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

import org.apache.log4j.Logger;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Implements Database functionality specific to Oracle (driver name and JDBC connection string)
 * @author Valerie Wagner
 *         Date: Apr 15, 2004
 */
public class OracleDatabase extends DatabaseProxy
{
    private static final Logger log = Logger.getLogger( OracleDatabase.class );


    OracleDatabase( String databaseHost, String databasePort, String databaseName )
            throws ClassNotFoundException
    {
        super( databaseHost, databasePort, databaseName );
    }


    /**
     * Initializes the driver name for Oracle.
     * Creates the database connection string for Oracle.
     */
    void initialize()
    {

//         - All database URLs must comply with the following generic syntax:
//           jdbc:oracle:driver_type:[username/password]@database_specifier

        //The URL syntax is as follows:
        //"jdbc:oracle:<driver>:@<db connection string>"
        //<driver>, can be 'thin' or 'oci8'
        //<db connect string>, is a Net8 name-value, denoting the TNSNAMES entry
        //Form the database connect string(TNSNAMES entry) as a name-value pair

        driverName = "oracle.jdbc.driver.OracleDriver";
        databaseURL = "jdbc:oracle:thin:@" +
                      "(DESCRIPTION=(ADDRESS=(HOST=" + databaseHost + ")" +
                      "(PROTOCOL=tcp)(PORT=" + databasePort + "))" +
                      "(CONNECT_DATA=(SID=" + databaseName + ")))";
        log.debug( "Oracle database url = " + databaseURL );

    }


    public int getDBMSType()
    {
        return Database.ORACLE;
    }


//    /**
//     * Inserts a string in a CLOB field.
//     *
//     * @param sqlLine String which is a Select query to get the CLOB field.
//     * @param value   The string with which to update the CLOB field once it is retrieved.
//     */
//
//    public void updateClob( String tableName, String columnName, long WID, String value ) throws SQLException
//    {
//        String sqlLine = "SELECT " + columnName + " FROM " + tableName +
//            " WHERE WID = " + WID + " for update nowait";
//
//        CLOB sequence = null;
//        log.debug( "Inserting a Clob" );
//
//        try
//        {
//            ResultSet rs = executeQuery( sqlLine );
//
//            rs.next();
//
//            sequence = ((OracleResultSet)rs).getCLOB( 1 );
//            rs.close();
//            rs = null;
//            int length = 0;
//            int buff_size = sequence.getBufferSize();
//
//            Writer out_clob = sequence.getCharacterOutputStream();
//            /*out_clob.write(value);
//            out_clob.flush();
//            out_clob.close();*/
//
//            int char_length = value.length();
//
//            while( (length + buff_size) < char_length )
//            {
//                out_clob.write( value, length, buff_size ); //write data[] in buff_size chunks
//                length += buff_size;
//            }
//            out_clob.write( value, length, (char_length - length) );
//            out_clob.flush();
//            out_clob.close();
//
//            commit();
//
//            /* Test to read the clob again. */
//            /*pstmt = biospiceCon.prepareStatement(sqlLine);
//            rs = pstmt.executeQuery(sqlLine);
//            rs.next();
//            sequence = ((OracleResultSet)rs).getCLOB(1);
//            if(sequence == null)
//               System.out.println(" While reading sequence is null");
//
//            rs.close();
//            String out = new String();
//                 System.out.println("The sequence length is " + sequence.length());
//            char[] buffer = new char[3000];
//                 long pos = 1;
//             sequence.getChars(pos,
//                          (new Long(sequence.length())).intValue(),
//                          buffer);
//                 String str = new String(buffer);
//                 System.out.println("The following sequence got saved: " + str); */
//
//        }
//        catch( IOException error )
//        {
//            log.error( "Problem filling clob", error );
//        }
//
//    }

     public ResultSet getTableMetaData( String tableName )
    {
        ResultSet results = null;

        try
        {
            // The "%" character matches all column names
            results = connection.getMetaData().getColumns( null, this.getSchemaName(), tableName.toUpperCase(), "%" );
            
        }
        catch( SQLException e )
        {
            log.error( "" );
        }

        return results;
    }

}
