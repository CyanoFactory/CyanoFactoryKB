/* $Id: LinkingTable.java,v 1.2 2006/11/03 06:44:08 kejariwa Exp $ */
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
package com.sri.biospice.warehouse.schema.linking;

import com.sri.biospice.warehouse.database.table.LongColumn;
import com.sri.biospice.warehouse.schema.TableBase;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * @author Valerie Wagner
 *         Date: Aug 13, 2004
 */
abstract class LinkingTable extends TableBase
{
    private LongColumn wid1Col;
    private LongColumn wid2Col;


    public LinkingTable( String table1Name, String table2Name )
    {
        super( table1Name + table2Name );
        wid1Col = new LongColumn( table1Name, this );
        wid2Col = new LongColumn( table2Name, this );
    }

    //AK - this constructor was added so that we can define the table name
    //and so the table is just the two columns concatenated
    public LinkingTable(String linkingTableName, String table1Name, String table2Name )
    {
        super(linkingTableName );
        wid1Col = new LongColumn( table1Name, this );
        wid2Col = new LongColumn( table2Name, this );
    }

    protected void setWID1( long wid )
    {
        wid1Col.setValue( wid );
    }


    protected void setWID2( long wid )
    {
        wid2Col.setValue( wid );
    }


    protected ResultSet getRow() throws SQLException
    {
        // todo: implement getRow()
        return null;
    }


}
