/* $Id: Reaction.java,v 1.1 2006/07/07 15:03:37 valerie Exp $ */
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
package com.sri.biospice.warehouse.schema.object;

import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Jan 5, 2005
 */
public class Reaction extends ObjectTableBase
{
    // Columns in the table
    private StringColumn deltaGCol;
    private StringColumn nameCol;
    private StringColumn ecNumberCol;
    private StringColumn ecNumberProposedCol;
    private StringColumn spontaneousCol;

    // Table MetaData
    private static TableMetaData tableMetaData = null;


    /**
     * Constructor for ObjectTableBase. Constructs a new object table with minimal required
     * initialization.
     * @param datasetWID required for all objects in the database
     */
    public Reaction( long datasetWID )
    {
        super( datasetWID, "Reaction" );
        init();
    }


    /**
     * Constructor for an object table that has already been stored in the Warehouse
     * @param datasetWID
     * @param thisWID
     */
    public Reaction( long datasetWID, long thisWID )
    {
        super( datasetWID, "Reaction", thisWID );
        init();
    }


    private void init()
    {
        deltaGCol = new StringColumn( "DeltaG", this );
        nameCol = new StringColumn( "Name", this );
        ecNumberCol = new StringColumn( "ECNumber", this );
        ecNumberProposedCol = new StringColumn( "ECNumberProposed", this );
        spontaneousCol = new StringColumn( "Spontaneous", this );
    }


    protected void setTableMetaData( TableMetaData tableMetaData )
    {
        this.tableMetaData = tableMetaData;
    }


    public TableMetaData getTableMetaData()
    {
        return tableMetaData;
    }


    public String getDeltaG()
    {
        return deltaGCol.getValue();
    }


    public void setDeltaG( String deltaG )
    {
        deltaGCol.setValue( deltaG );
    }


    public String getECNumber()
    {
        return ecNumberCol.getValue();
    }


    public void setECNumber( String ecNumber )
    {
        ecNumberCol.setValue( ecNumber );
    }


    public String getECNumberProposed()
    {
        return ecNumberProposedCol.getValue();
    }


    public void setECNumberProposed( String ecNumberProposed )
    {
        ecNumberProposedCol.setValue( ecNumberProposed );
    }


    public String getSpontaneous()
    {
        return spontaneousCol.getValue();
    }


    public void setSpontaneous( String spontaneous )
    {
        spontaneousCol.setValue( spontaneous );
    }
}
