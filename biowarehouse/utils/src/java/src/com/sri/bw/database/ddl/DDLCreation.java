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
import com.sri.bw.dbschema.ColumnDocument.Column;
import com.sri.bw.dbschema.Enum;
import com.sri.bw.dbschema.ForeignKeyDocument.ForeignKey;
import com.sri.bw.dbschema.TableDocument.Table;
import com.sri.biospice.warehouse.database.DatabaseProxy;
import org.apache.xmlbeans.XmlAnySimpleType;
import org.apache.xmlbeans.XmlObject;
import java.io.File;
import java.io.PrintStream;
import java.util.Hashtable;
import java.util.Vector;

/**
 * Translates DB schema in XML format to dialect-specific RDBMS DDL scripts
 *
 * @author Valerie Wagner
 *         Date: Jul 28, 2005
 */
public abstract class DDLCreation {

    private PrintStream createScriptOS;
    private PrintStream indexScriptOS;
    private PrintStream destroyScriptOS;
    private Vector<String> foreignKeyDefinitionLines = new Vector<String>();
    private String dbmsName;
    private Hashtable<String, String> columnTypes = null;
    private SchemaDocument.Schema schemaRoot;

    public static void main(String[] args) {
        // Oracle
        DDLCreation ddloracle = new OracleDDLCreation(args);
        ddloracle.translate();

        // MySQL
        DDLCreation ddlmysql = new MySQLDDLCreation(args);
        ddlmysql.translate();
    }


    public DDLCreation(String[] args, String dbmsName) {

        if (args.length != 2) {
            System.err.println("Wrong number of arguments.");
            System.err.println("Usage: DDLCreation inputfile output-directory");
            System.exit(1);
        }
        String inputFilename = args[0];
        String outputDir = args[1];
        try {
            this.dbmsName = dbmsName;
            // Set up output for DDL Script
            File createFile = new File(outputDir, "warehouse-" + dbmsName + "-create.sql");
            if (!createFile.exists()) {
                createFile.createNewFile();
            }
            createScriptOS = new PrintStream(createFile);

            // Set up output for Index Script
            File indexFile = new File(outputDir, "warehouse-" + dbmsName + "-index.sql");
            if (!indexFile.exists()) {
                indexFile.createNewFile();
            }
            indexScriptOS = new PrintStream(indexFile);

            // Set up output for Destroy Script
            File destroyFile = new File(outputDir, "warehouse-" + dbmsName + "-destroy.sql");
            if (!destroyFile.exists()) {
                destroyFile.createNewFile();
            }
            destroyScriptOS = new PrintStream(destroyFile);

            // Read the input file
            File inputFile = new File(inputFilename);
            System.out.println("Translating file " + inputFile + " to:");
            System.out.println(createFile);
            System.out.println(indexFile);
            System.out.println(destroyFile);

            schemaRoot = SchemaDocument.Factory.parse(inputFile).getSchema();

            // Set up the type map
            columnTypes = new Hashtable<String, String>();
            DatatypesDocument.Datatypes[] types = schemaRoot.getDatatypesArray();
            for (DtypeDocument.Dtype type : types[0].getDtypeArray()) {
                columnTypes.put(type.getNotation(), getType(type));
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Translates XML schema into DDL scripts
     */
    private void translate() {
        // Traverse through nodes in the schema document
        XmlObject[] nodes = schemaRoot.selectPath("*");
        for (int i = 0; i < nodes.length; i++) {
            XmlObject node = nodes[i];
            if (node instanceof TableDocument.Table) {
                processTable((TableDocument.Table)node);
            } else if (node instanceof SequenceDocument.Sequence) {
                processSequence((SequenceDocument.Sequence)node);
            } else if (node instanceof CommandDocument.Command) {
                CommandDocument.Command cmd = (CommandDocument.Command)node;
                if (cmd.isSetDialect()) {
                    if (cmd.getDialect().equals(dbmsName)) {
                        output(cmd.getStringValue() + "\n");
                    }
                } else {
                    output(cmd.getStringValue() + "\n");
                }
            } else if (node instanceof XmlAnySimpleType) {
                output(((XmlAnySimpleType)node).getStringValue() + "\n");
            } else {
//                System.out.println("unknown node type: " + node.getClass());
            }
        }

        writeForeignKeys();

        writeDestroyScript();

        writeCommit();

        // Close output files
        createScriptOS.flush();
        createScriptOS.close();
        indexScriptOS.flush();
        indexScriptOS.close();
        destroyScriptOS.flush();
        destroyScriptOS.close();
    }

    private void writeCommit() {
        createScriptOS.println("\ncommit;");
        indexScriptOS.println("\ncommit;");
        destroyScriptOS.println("\ncommit;");
    }

    private void writeDestroyScript() {
        destroyScriptOS.print(startDestroyScript());
        for (Table table : schemaRoot.getTableArray()) {
            for (Column column : table.getColumnArray()) {
                for (ForeignKey key : column.getForeignKeyArray()) {
                    destroyScriptOS.print(getDropFKStatement(table, key));
                }
            }
        }
        for (Table table : schemaRoot.getTableArray()) {
            destroyScriptOS.println("drop table " + table.getName() + ";");
        }

        for (SequenceDocument.Sequence sequence : schemaRoot.getSequenceArray()) {
            destroyScriptOS.println(getDropSequence(sequence));
        }

        destroyScriptOS.println(finishDestroyScript());
    }


    private void processTable(TableDocument.Table table) {
        Vector<String> columnDefinitionLines = new Vector<String>();
        Vector<String> otherConstraintLines = new Vector<String>();
        Vector<String> primaryKeyDefinitionLine = new Vector<String>();

        String line;

        output("CREATE TABLE " + table.getName());
        output("(");

        // Iterate through columns
        ColumnDocument.Column[] columns = table.getColumnArray();
        for (int i = 0; i < columns.length; i++) {

            // Column definitions
            // WaterSolubility         CHAR(1),                  -- 'T' if soluble in water, else 'F'
            ColumnDocument.Column column = columns[i];
            line = "  " + column.getName() + "  ";
            line += determineType(column);
            line += columnRequired(column.getRequired());
            line += autoIncrement(column.getAutoIncrement());
            columnDefinitionLines.add(line);

            /* Foreign key definitions
              ALTER TABLE TT_PhysicalArrayDesign
              ADD CONSTRAINT FK_1
              FOREIGN KEY (surfaceType_ID)
              REFERENCES TT_OntologyEntry (ID)
             */
            for (ForeignKeyDocument.ForeignKey fkey : column.getForeignKeyArray()) {
                line = getFKLine(table, column, fkey);
                if (line != null) {
                    foreignKeyDefinitionLines.add(line);
                }
            }
        }

        // Other constraints
        String[] constraints = table.getConstraintArray();
        for (int i = 0; i < constraints.length; i++) {
            String constraint = constraints[i];
            line = "  " + constraint;
            otherConstraintLines.add(line);
        }

        // Primary Key definitions
        // CONSTRAINT PK_Chemical1 PRIMARY KEY (WID) USING INDEX TABLESPACE INDEXES
        if (table.isSetPrimaryKey()) {
            PrimaryKeyDocument.PrimaryKey pk = table.getPrimaryKey();
            primaryKeyDefinitionLine.add(getPKConstraint(pk));
        }

        int totalLines = columnDefinitionLines.size() + otherConstraintLines.size() + primaryKeyDefinitionLine.size();
        int lineIndex = 1;
        for (String oneLine : columnDefinitionLines) {
            if (lineIndex++ < totalLines) {
                oneLine += ",";
            }
            output(oneLine);
        }
        output("  --");

        for (String oneLine : otherConstraintLines) {
            if (lineIndex++ < totalLines) {
                oneLine += ",";
            }
            output(oneLine);
        }

        for (String oneLine : primaryKeyDefinitionLine) {
            if (lineIndex++ < totalLines) {
                oneLine += ",";
            }
            output(oneLine);
        }

        // SD comments
//        String[] sds = table.getSdArray();
//        for (int i = 0; i < sds.length; i++) {
//            String sd = sds[i];
//            output("  -- " + sd);
//        }

        output(")" + getEndTable() + ";");
        output("");

        // Index
        IndexDocument.Index[] indexes = table.getIndexArray();
        for (int i = 0; i < indexes.length; i++) {
            IndexDocument.Index index = indexes[i];
            //CREATE INDEX CHEMICAL_DWID_NAME ON Chemical(DataSetWID, Name) TABLESPACE INDEXES;
            String indexStmt = "CREATE INDEX " + index.getName() + " ON " + table.getName() + "(" + getIndexColumns(index) + ")" + getIndexEnding() + ";\n";
            if (index.isSetInitial() && index.getInitial().equalsIgnoreCase("true")) {
                createScriptOS.print(indexStmt);
            } else {
                indexScriptOS.print(indexStmt);
            }
        }
        if (indexes.length > 0) {
            indexScriptOS.println();
        }
        output("");

        // Add enumerated values, if any
        for (Column col : table.getColumnArray()) {
            for (Enumeration enumset : col.getEnumerationArray()) {
                for (Enum rest : enumset.getRestrictionArray()) {
                    output("INSERT INTO Enumeration (TableName, ColumnName, Value, Meaning) VALUES ('" +
                            table.getName() + "', '" +
                            col.getName() + "', '" +
                            rest.getValue() + "', '" +
                            DatabaseProxy.escapeSingleQuotes(rest.getDescription()) + "');");
                }
            }
        }
        output("");
    }

    private String getIndexColumns(IndexDocument.Index index) {
        if (index.getColumns() != null) {
            return index.getColumns();
        }
        for (IndexDocument.Index.Variant variant : index.getVariantArray()) {
            if (variant.getDialect().equals(dbmsName)) {
                return variant.getColumns();
            }
        }

        return null;
    }


    private String determineType(ColumnDocument.Column column) {
        String type = column.getType();

        if (columnTypes.containsKey(type)) {
            type = columnTypes.get(type);
            if (column.isSetLength()) {
                type += "(" + column.getLength() + ")";
            } else if (column.isSetPrecision() && column.isSetScale()) {
                type += "(" + column.getPrecision() + "," + column.getScale() + ")";
            }
        }

        return type;
    }

    private String columnRequired(String required) {
        if ("true".equals(required)) {
            return "  NOT NULL";
        } else {
            return "";
        }
    }


    private void writeForeignKeys() {
        for (String line : foreignKeyDefinitionLines) {
            output(line);
            output("");
        }
    }

    protected void output(String text) {
//        System.out.println(text);
        createScriptOS.print(text + "\n");
        createScriptOS.flush();
    }

    protected abstract String getType(DtypeDocument.Dtype type);

    protected abstract String getIndexEnding();

    protected abstract String autoIncrement(String autoIncrement);

    protected abstract String getFKLine(TableDocument.Table table, ColumnDocument.Column column, ForeignKeyDocument.ForeignKey fkey);

    protected abstract String getEndTable();

    protected abstract String getPKConstraint(PrimaryKeyDocument.PrimaryKey pk);

    protected abstract String getDropConstraint();

    protected abstract void processSequence(SequenceDocument.Sequence sequence);

    protected abstract String startDestroyScript();

    protected abstract String finishDestroyScript();

    protected abstract String getDropSequence(SequenceDocument.Sequence sequence);

    protected abstract String getDropFKStatement(Table table, ForeignKey key);


}
