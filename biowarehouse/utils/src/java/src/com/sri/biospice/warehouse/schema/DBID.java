/* $Id: DBID.java,v 1.1 2006/07/07 15:03:36 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.table.LongColumn;
import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;
import com.sri.biospice.warehouse.schema.object.ObjectTable;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Vector;

/**
 * An entry in the DBID table
 * @author Priyanka Gupta
 * @author Valerie Wagner
 */

public class DBID extends TableBase
{
    // Columns of the table
    private LongColumn otherWIDCol;
    private StringColumn xidCol;
    private StringColumn typeCol;
    private StringColumn versionCol;

    // Table metadata
    private static TableMetaData tableMetaData;

    // Enumerated values used by the Type column
//    public static final String ENUM_TYPE_ACCESSION = "Accession";
//    public static final String ENUM_TYPE_GUID = "GUID";

    /**
     * Enumerated values allowed for column <i>DBID.Type</i>.
     */
    public enum Type {
        Accession,
        GUID
    }

    // Hashcode used for comparing DBIDs when placed in a hash-based data structure
    private int hashCode = 0;


    /**
     * Constructor for the class. Used to store a unique
     * identifier for an object from this database.
     * @param xid      An identifier for an object in this DataSet.
     * @param type     The type of the identifier.
     *                 i.e what the identifier represents.
     * @param version  The version of the xid.
     * @param otherWid The wid of the warehouse object that this
     *                 gets linked from.
     * @deprecated Use the typed version instead
     */
    public DBID( String xid, String type, String version, long otherWid )
    {
        super( "DBID" );
        init();
        xidCol.setValue( xid );
        typeCol.setValue( type );
        versionCol.setValue( version );
        otherWIDCol.setValue( otherWid );
    }

    /**
     * Constructor for the class. Used to store a unique
     * identifier for an object from this database.
     * @param xid      An identifier for an object in this DataSet.
     * @param type     The type of the identifier.
     *                 i.e what the identifier represents.
     * @param version  The version of the xid.
     * @param otherWid The wid of the warehouse object that this
     *                 gets linked from.
     */
    public DBID( String xid, Type type, String version, long otherWid )
    {
        super( "DBID" );
        init();
        xidCol.setValue( xid );
        typeCol.setValue( type );
        versionCol.setValue( version );
        otherWIDCol.setValue( otherWid );
    }


    public DBID()
    {
        super( "DBID" );
        init();
    }


    /**
     * Initializes columns
     */
    private void init()
    {
        otherWIDCol = new LongColumn( "OtherWID", this );
        xidCol = new StringColumn( "XID", this );
        typeCol = new StringColumn( "Type", this );
        versionCol = new StringColumn( "Version", this );
    }


    public TableMetaData getTableMetaData()
    {
        return tableMetaData;
    }


    protected void setTableMetaData( TableMetaData tableMetaData )
    {
        this.tableMetaData = tableMetaData;
    }


    protected ResultSet getRow() throws SQLException
    {
        return null;
    }


    public Enum [] getAllowedValues( String columnName )
    {
        if( columnName.equalsIgnoreCase( typeCol.getColumnName() ) )
        {
            return Type.values();
        }

        return null;
    }


    public boolean equals( Object obj )
    {
        if( obj == this )
        {
            return true;
        }
        else if( !(obj instanceof DBID) )
        {
            return false;
        }
        DBID dbid = (DBID)obj;
        return dbid.otherWIDCol.equals( this.otherWIDCol ) &&
               dbid.typeCol.equals( this.typeCol ) &&
               dbid.versionCol.equals( this.versionCol ) &&
               dbid.xidCol.equals( this.xidCol );
    }


    public int hashCode()
    {
        if( hashCode == 0 )
        {
            int result = 23;
            result = 37 * result + otherWIDCol.hashCode();
            result = 37 * result + typeCol.hashCode();
            result = 37 * result + versionCol.hashCode();
            result = 37 * result + xidCol.hashCode();
            hashCode = result;
        }
        return hashCode;
    }

    public long getOtherWID() {
        return otherWIDCol.getValue();
    }

    public String getXID() {
        return xidCol.getValue();
    }

    public String getType() {
        return typeCol.getValue();
    }

    public String getVersion() {
        return versionCol.getValue();
    }

    public static Vector loadDBIDs(long otherWID) {
        return TableFactory.loadTables("select * from DBID where OtherWID='" + otherWID + "'", TableFactory.DBID);
    }
}
