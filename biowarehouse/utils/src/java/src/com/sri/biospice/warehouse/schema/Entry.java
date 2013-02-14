/******************************************************************************
 * Copyright 2004 SRI International.  All rights reserved.                     *
 *                                                                             *
 * The material contained in this file is confidential and proprietary to SRI  *
 * International and may not be reproduced, published, or disclosed to others  *
 * without authorization from SRI International.                               *
 *                                                                             *
 * DISCLAIMER OF WARRANTIES                                                    *
 *                                                                             *
 * SRI International MAKES NO REPRESENTATIONS OR WARRANTIES ABOUT THE          *
 * SUITABILITY OF THE SOFTWARE, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT   *
 * LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A         *
 * PARTICULAR PURPOSE, OR NON-INFRINGEMENT. SRI International SHALL NOT BE     *
 * LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A RESULT OF USING, MODIFYING *
 * OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES                            *
 *******************************************************************************/


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
package com.sri.biospice.warehouse.schema;

import com.sri.biospice.warehouse.database.table.*;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Date;
import java.util.Vector;

/**
 * An entry in the Entry table
 * @author Liviu Popescu
 * @author Valerie Wagner
 */
public class Entry extends TableBase
{
    // Columns in the Entry table
    private LongColumn otherWIDCol;
    private TimestampColumn insertDateCol;
    private TimestampColumn creationDateCol;
    private TimestampColumn modifiedDateCol;
    private StringColumn loadErrorCol;
    private IntegerColumn lineNumberCol;
    private StringColumn errorMessageCol;
    private LongColumn datasetWIDCol;

    private static TableMetaData tableMetaData;


    /**
     * Constructor for Table Entry. Takes as arguments the WID of the object it represents and the dataset WID
     * and sets the insert date to the current date and Load Error to F. If an error appears to set the flag
     * the Function @see com.sri.biospice.warehouse.schema.Entry#setLoadError()
     * @param otherWID
     * @param datasetWID
     */
    public Entry( long otherWID, long datasetWID )
    {
        super( "Entry" );
        init();
        otherWIDCol.setValue( otherWID );
        insertDateCol.setValue( new Date() );
        loadErrorCol.setValue( "F" );
        datasetWIDCol.setValue( datasetWID );
    }


    /**
     * Constructor used for entries already in the database
     * @param otherWID
     */
    public Entry( long otherWID )
    {
        super( "Entry" );
        init();
        otherWIDCol.setValue( otherWID );
    }

    public Entry( )
    {
        super( "Entry" );
        init();
    }


    private void init()
    {
        otherWIDCol = new LongColumn( "OtherWID", this );
        insertDateCol = new TimestampColumn( "InsertDate", this );
        creationDateCol = new TimestampColumn( "CreationDate", this );
        modifiedDateCol = new TimestampColumn( "ModifiedDate", this );
        loadErrorCol = new StringColumn( "LoadError", this );
        lineNumberCol = new IntegerColumn( "LineNumber", this );
        errorMessageCol = new StringColumn( "ErrorMessage", this );
        datasetWIDCol = new LongColumn( "DatasetWID", this );
    }


    protected void setTableMetaData( TableMetaData tableMetaData )
    {
        this.tableMetaData = tableMetaData;
    }


    public TableMetaData getTableMetaData()
    {
        return tableMetaData;
    }


    protected ResultSet getRow() throws SQLException
    {
        // todo: implement getRow()
        return null;
    }


    public void setCreationDate( Date creationDate )
    {
        creationDateCol.setValue( creationDate );
    }


    /**
     * Sets the creation date to be the same as the insertion date
     */
    public void setCreationDate()
    {
        creationDateCol.setValue( insertDateCol.getValue() );
    }


    public void setErrorMessage( String errorMessage )
    {
        errorMessageCol.setValue( errorMessage );
    }


    public void setLineNumber( String lineNumber )
    {
        lineNumberCol.setValue( lineNumber );
    }


    public void setLineNumber( Integer lineNumber )
    {
        lineNumberCol.setValue( lineNumber );
    }


    /**
     * Set the modified date to be the same as the insert date
     */
    public void setModifiedDate()
    {
        modifiedDateCol.setValue( insertDateCol.getValue() );
    }


    public void setModifiedDate( Date modifiedDate )
    {
        modifiedDateCol.setValue( modifiedDate );
    }


    /**
     * Set the Load error to T
     */
    public void setLoadError()
    {
        loadErrorCol.setValue( "T" );
    }

    public static Vector loadEntries(long otherWID) {
        return TableFactory.loadTables("select * from Entry where OtherWID='" + otherWID + "'", TableFactory.ENTRY);
    }
}
