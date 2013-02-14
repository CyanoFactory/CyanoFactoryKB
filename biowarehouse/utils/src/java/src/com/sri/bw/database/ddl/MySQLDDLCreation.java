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
public class MySQLDDLCreation extends DDLCreation {

    public MySQLDDLCreation(String [] args) {
        super(args, "mysql");
    }

    protected String getEndTable() {
        return " TYPE=INNODB";
    }

    protected String getFKLine(TableDocument.Table table, ColumnDocument.Column column, ForeignKeyDocument.ForeignKey fkey) {
        return "ALTER TABLE " + table.getName() +
                "\nADD CONSTRAINT " + fkey.getName() +
                "\nFOREIGN KEY (" + column.getName() + ")" +
                "\nREFERENCES " + fkey.getToTable() + " (" + fkey.getToColumn() + ") ON DELETE CASCADE;";
    }


    protected String getPKConstraint(PrimaryKeyDocument.PrimaryKey pk) {
        return "  CONSTRAINT " + pk.getName() + " PRIMARY KEY (" + pk.getColumn() + ")";
    }

    protected String getDropConstraint() {
        return "foreign key";
    }

    protected void processSequence(SequenceDocument.Sequence sequence) {
    }

    protected String startDestroyScript() {
        return "SET FOREIGN_KEY_CHECKS = 0;\n\n";
    }

    protected String finishDestroyScript() {
        return "SET FOREIGN_KEY_CHECKS = 1;\n";
    }

    protected String getDropSequence(SequenceDocument.Sequence sequence) {
        return "";
    }

    protected String getDropFKStatement(TableDocument.Table table, ForeignKeyDocument.ForeignKey key) {
        return "";
    }

    protected String getType(DtypeDocument.Dtype type) {
        return type.getMysql();
    }

    protected String getIndexEnding() {
        return "";
    }

    public String autoIncrement(String autoIncrement) {
        if (autoIncrement != null && autoIncrement.equals("true")) {
            return " AUTO_INCREMENT";
        }
        return "";
    }

}
