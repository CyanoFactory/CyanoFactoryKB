/* $Id: CommentTable.java,v 1.1 2006/07/07 15:03:36 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.table.ClobColumn;
import com.sri.biospice.warehouse.database.table.LongColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;
import com.sri.biospice.warehouse.schema.object.ObjectTable;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Vector;

/**
 * An entry in the CommentTable table
 *
 * @author Valerie Wagner
 *         Date: Aug 13, 2004
 */
public class CommentTable extends TableBase {
    private static TableMetaData tableMetaData;
    private LongColumn otherWIDCol;
    private ClobColumn commentCol;


    public CommentTable(ObjectTable refTable, String comment) {
        super("CommentTable");
        init();
        this.setComment(comment);
        this.setOtherWID(refTable.getWID());
    }

    public CommentTable() {
        super("CommentTable");
        init();
    }


    private void init() {
        otherWIDCol = new LongColumn("OtherWID", this);
        commentCol = new ClobColumn("Comm", this);
    }


    public String getComment() {
        return commentCol.getValue();
    }


    public void setComment(String comm) {
        commentCol.setValue(comm);
    }


    public long getOtherWID() {
        return otherWIDCol.getValue();
    }


    public void setOtherWID(long otherWID) {
        otherWIDCol.setValue(otherWID);
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
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

    public static Vector loadCommentTables(long otherWID) {
        return TableFactory.loadTables("select * from CommentTable where OtherWID='" + otherWID + "'", TableFactory.COMMENT_TABLE);
    }
}
