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
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Apr 7, 2005
 */
public class TermRelationship extends TableBase {

    public static final String ENUM_RELATIONSHIP_IS_A = "is_a";
    public static final String ENUM_RELATIONSHIP_PART_OF = "part_of";

    private static TableMetaData tableMetaData;
    private LongColumn termWIDCol;
    private LongColumn relatedTermWIDCol;
    private StringColumn relationshipCol;

    public TermRelationship() {
        super("TermRelationship");
        init();
    }

    public TermRelationship(long termWID, long relatedTermWID, String relationship) {
        super("TermRelationship");
        init();
        termWIDCol.setValue(termWID);
        relatedTermWIDCol.setValue(relatedTermWID);
        relationshipCol.setValue(relationship);
    }

    private void init() {
        termWIDCol = new LongColumn("TermWID", this);
        relatedTermWIDCol = new LongColumn("RelatedTermWID", this);
        relationshipCol = new StringColumn("Relationship", this);
    }

    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }

    /**
     * todo not implemented
     *
     * @return
     * @throws SQLException
     */
    protected ResultSet getRow() throws SQLException {
        return null;
    }

    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }

    public static Vector loadTermRelationships(long termWID) {
        return TableFactory.loadTables("select * from TermRelationship where TermWID='" + termWID + "'", TableFactory.TERM_RELATIONSHIP);
    }
}
