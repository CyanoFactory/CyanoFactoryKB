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
package com.sri.bw.database.ddl;

import com.sri.bw.dbschema.*;

/**
 * @author Valerie Wagner
 *         Date: Feb 21, 2006
 */
public class OracleDDLCreation extends DDLCreation {


    public OracleDDLCreation(String[] args) {
        super(args, "oracle");
    }

    protected String getEndTable() {
        return "";
    }

    protected String getFKLine(TableDocument.Table table, ColumnDocument.Column column, ForeignKeyDocument.ForeignKey fkey) {
        return "ALTER TABLE " + table.getName() +
                "\nADD CONSTRAINT " + fkey.getName() +
                "\nFOREIGN KEY (" + column.getName() + ")" +
                "\nREFERENCES " + fkey.getToTable() + " (" + fkey.getToColumn() + ") ON DELETE CASCADE;";
    }


    protected String getPKConstraint(PrimaryKeyDocument.PrimaryKey pk) {
        return "  CONSTRAINT " + pk.getName() + " PRIMARY KEY (" + pk.getColumn() + ") USING INDEX TABLESPACE INDEXES";
    }

    protected String getDropConstraint() {
        return "constraint";
    }

    public String autoIncrement(String autoIncrement) {
        return "";
    }

    protected String getDropSequence(SequenceDocument.Sequence sequence) {
        return "drop SEQUENCE " + sequence.getName() + ";";
    }

    protected String getDropFKStatement(TableDocument.Table table, ForeignKeyDocument.ForeignKey key) {
        return "alter table " + table.getName() +
                " drop " + getDropConstraint() + " " + key.getName() + ";\n";
    }

    protected String getType(DtypeDocument.Dtype type) {
        return type.getOracle();
    }

    protected String getIndexEnding() {
        return " TABLESPACE INDEXES";
    }

    protected void processSequence(SequenceDocument.Sequence sequence) {
        String[] comments = sequence.getCommentArray();
        for (int i = 0; i < comments.length; i++) {
            String comment = comments[i];
            output("-- " + comment);
        }
        output("--");
//        CREATE SEQUENCE WID_sequence INCREMENT BY 1 START WITH 1000 MINVALUE 1000;
        output("CREATE SEQUENCE " + sequence.getName() + " INCREMENT BY " + sequence.getIncrement() +
                " START WITH " + sequence.getStart() + " MINVALUE " + sequence.getMinValue() + ";\n");
    }

    protected String startDestroyScript() {
        return "";
    }

    protected String finishDestroyScript() {

        String endDestroy = "-- Destroy the capability to flash restore tables, reclaiming space.\n" +
                "-- This only works for Oracle 9+\n" +
                "PURGE RECYCLEBIN;";
        return endDestroy;
    }

}
