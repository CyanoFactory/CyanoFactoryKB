/* $Id: Validation.java,v 1.1 2006/10/27 17:23:03 valerie Exp $ */
package com.sri.biospice.warehouse.util;

import com.sri.bw.dbschema.ColumnDocument;
import com.sri.bw.dbschema.ForeignKeyDocument;
import com.sri.bw.dbschema.SchemaDocument;
import com.sri.bw.dbschema.TableDocument;
import com.sri.bw.element.*;
import org.apache.log4j.Logger;
import org.apache.log4j.BasicConfigurator;
import org.apache.xmlbeans.XmlCursor;
import org.apache.xmlbeans.XmlException;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.File;
import java.util.Enumeration;
import java.util.Properties;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Oct 27, 2006
 */
public class Validation {

    private static final Logger log = Logger.getLogger(Validation.class);

    public static final int MAX_IDENTIFIER_LENGTH = 30;

    public static void main(String[] args) {
        BasicConfigurator.configure();
        if (args.length != 1) {
            System.out.println("Usage: Validation path-to-schema-file");
            System.exit(1);
        }
        String filename = args[0];
        if (LoaderMain.fileIsReadable(filename)) {

            try {
                SchemaDocument schemaDoc = SchemaDocument.Factory.parse(new File(filename));
                Validation validation = new Validation();
                validation.validateDatabaseSchema(schemaDoc.getSchema());
            } catch (XmlException e) {
                fail("Error reading schema file", e);
            } catch (IOException e) {
                fail("Error reading schema file", e);                
            }
        } else {
            System.out.println("Cannot read file: " + filename);
        }
    }

    private void validateElementMapping(ModelDocument.Model model, SchemaDocument.Schema schema) {
        for (ElementDocument.Element element : model.getElementArray()) {
            if (element.isSetSkip2() && element.getSkip2().equals("true")) {
                continue;
            }
            String elementName = element.getName();
            String tableName = element.getTable();
            TableDocument.Table elementTable = null;
            Vector<String> associationNames = new Vector<String>();
            if (!SchemaUtils.hasTable(schema, tableName)) {
                fail("MAGE element maps to a DB table that does not exist: " + element.getName() + "-->" + tableName);
                continue;
            } else {
                elementTable = SchemaUtils.findTable(schema, tableName);
                for (Attribute att : element.getAttributeArray()) {
                    associationNames.add(att.getName());
                    if (!(att.isSetSkip() && att.getSkip().equals("true"))) {
                        if (!SchemaUtils.hasColumn(elementTable, att.getColumn())) {
                            fail("MAGE attribute maps to a DB column that does not exist: " +
                                    element.getName() + "." + att.getName() + "-->" + tableName + "." + att.getColumn());
                        }
                    }
                }
            }

            for (OneToOne one2one : element.getOneToOneArray()) {
                tableName = one2one.getTarget();
                associationNames.add(one2one.getName());
                if (!(one2one.isSetSkip() && one2one.getSkip().equals("true"))) {
                    if (!SchemaUtils.hasTable(schema, tableName)) {
                        fail("MAGE " + elementName + " one2one maps to a DB table that does not exist: " + one2one.getName() + "-->" + tableName);
                    } else if (!SchemaUtils.hasColumn(elementTable, one2one.getColumn())) {
                        fail("MAGE " + elementName + " one2one maps to a DB column that does not exist: " +
                                one2one.getName() + "-->" + elementTable.getName() + "." + one2one.getColumn());
                    }
                }
            }

            for (OneToMany one2many : element.getOneToManyArray()) {
                tableName = one2many.getTarget();
                associationNames.add(one2many.getName());
                if (!(one2many.isSetSkip() && one2many.getSkip().equals("true"))) {
                    TableDocument.Table table = SchemaUtils.findTable(schema, tableName);
                    if (!SchemaUtils.hasTable(schema, tableName)) {
                        fail("MAGE " + elementName + " one2many maps to a DB table that does not exist: " + one2many.getName() + "-->" + tableName);
                    } else if (!SchemaUtils.hasColumn(table, one2many.getColumn())) {
                        fail("MAGE " + elementName + " one2many maps to a DB column that does not exist: " +
                                one2many.getName() + "-->" + table.getName() + "." + one2many.getColumn());
                    }
                }
            }

            for (ManyToOne many2one : element.getManyToOneArray()) {
                tableName = many2one.getTarget();
                associationNames.add(many2one.getName());
                if (!(many2one.isSetSkip() && many2one.getSkip().equals("true"))) {
                    if (!SchemaUtils.hasTable(schema, tableName)) {
                        fail("MAGE " + elementName + " many2one maps to a DB table that does not exist: " + many2one.getName() + "-->" + tableName);
                    }
                    if (!SchemaUtils.hasColumn(elementTable, many2one.getColumn())) {
                        fail("MAGE " + elementName + " many2one maps to a DB column that does not exist: " +
                                many2one.getName() + "-->" + elementTable.getName() + "." + many2one.getColumn());
                    }
                }
            }

            for (ManyToMany many2many : element.getManyToManyArray()) {
                // Ensure we don't accidentally use this
                elementTable = null;
                if (!SchemaUtils.hasTable(schema, many2many.getTable())) {
                    fail("MAGE " + elementName + " many2many maps to a DB table that does not exist: " + many2many.getName() + "-->" + tableName);
                }
                TableDocument.Table linkingTable = SchemaUtils.findTable(schema, many2many.getTable());
                if (!SchemaUtils.hasColumn(linkingTable, many2many.getColumn())) {
                    fail("MAGE " + elementName + " many2many maps to a DB column that does not exist: " +
                            many2many.getName() + "-->" + linkingTable.getName() + "." + many2many.getColumn());
                }
                if (!SchemaUtils.hasColumn(linkingTable, many2many.getCrosscolumn())) {
                    fail("MAGE " + elementName + " many2many maps to a DB column that does not exist: " +
                            many2many.getName() + "-->" + linkingTable.getName() + "." + many2many.getCrosscolumn());
                }
            }

            checkForDuplicateNames(associationNames);
        }
    }

    /**
     * Check that column and table names don't exceed allowed length (30)
     * Check for duplicate column and table  names
     */
    private void validateDatabaseSchema(SchemaDocument.Schema schema) {

        log.info("Checking that XML schema meets validation requirements...");
        TableDocument.Table[] tables = schema.getTableArray();
        Vector<String> tableNames = new Vector<String>();

        for (TableDocument.Table table : tables) {
            tableNames.add(table.getName());
            if (table.getName().length() > MAX_IDENTIFIER_LENGTH) {
                fail("Table name exceeds max length: " + table.getName());
            }
            if (databaseReservedWords.contains(table.getName())) {
                fail("Table name identifier is a reserved word: " + table.getName());
            }

            if (table.getType().equals("object")) {
                if (!table.isSetPrimaryKey()) {
                    fail("Object table does not have primary key defined: " + table.getName());
                }
            }

            ColumnDocument.Column[] columns = table.getColumnArray();
            Vector<String> columnNames = new Vector<String>();
            for (ColumnDocument.Column column : columns) {
                columnNames.add(table.getName() + "." + column.getName());
                if (column.getName().length() > MAX_IDENTIFIER_LENGTH) {
                    fail("Column name exceeds max length: " + table.getName() + "." + column.getName());
                }
                if (databaseReservedWords.contains(column.getName())) {
                    fail("Column name identifier is a reserved word: " + table.getName() + "." + column.getName());
                }

                if (column.getCommentArray().length == 0) {
                    warn(table.getName() + "." + column.getName() + " does not have a comment");
                }

                // Make sure foreign keys point to objects that exist
                for (ForeignKeyDocument.ForeignKey key : column.getForeignKeyArray()) {
                    TableDocument.Table keyTable = SchemaUtils.findTable(schema, key.getToTable());
                    if (keyTable == null) {
                        fail("Foreign Key points to table that does not exist: " + key.getName() + ", " + key.getToTable());
                    } else {
                        if (!SchemaUtils.hasColumn(keyTable, key.getToColumn())) {
                            fail("Foreign Key points to a column that does not exist: " + key.getName() +
                                    ", " + key.getToTable() + "." + key.getToColumn());
                        }
                    }
                }

                // Make sure WID columns that should have a foreign key have one
                // todo: add back in.  determine other WID-type columns and their behavior
//                if (column.getType().equals(WID_TYPE)) {
//                    if (table.isSetPrimaryKey() && table.getPrimaryKey().getColumn().equals(column.getName())) {
//                        // ok -- primary keys don't have to map to a foreign key
//                    } else if (column.getForeignKeyArray().length == 0) {
//                        fail("Column is WID type yet does not define foreign key: " + table.getName() + "." + column.getName());
//                    }
//                }
            }

            checkForDuplicateNames(columnNames);
        }
        checkForDuplicateNames(tableNames);

        // check for globally unique foreign key names
        Vector<String> fkNames = new Vector<String>();
        XmlCursor cursor = schema.newCursor();
        cursor.selectPath(".//foreignKey");
        while (cursor.toNextSelection()) {
            fkNames.add(((ForeignKeyDocument.ForeignKey)cursor.getObject()).getName());
        }
        checkForDuplicateNames(fkNames);
        log.info("Validation successful.");
    }

    private void warn(String message) {
        log.warn(message);
    }

    /**
     * Checks for duplicate names in a list of names
     *
     * @param names The list of names
     */
    private void checkForDuplicateNames(Vector<String> names) {
        for (int i = 0; i < names.size(); i++) {
            for (int j = i + 1; j < names.size(); j++) {
                if (names.get(i).equals(names.get(j))) {
                    fail("Duplicate name: " + names.get(i));
                }
            }
        }
    }


    private static void fail(String message) {
        log.fatal(message);
        System.exit(1);
    }

    private static void fail(String message, Exception e) {
        log.fatal(message, e);
        System.exit(1);
    }


    public static Vector<String> loadReservedWords() {
        loadReservedWordsFile("oracle-reserved-words.txt");
        loadReservedWordsFile("mysql-reserved-words.txt");
        return databaseReservedWords;
    }

    private static void loadReservedWordsFile(String filename) {
        InputStream owis = Validation.class.getResourceAsStream("/com/sri/bw/mage/translation/" + filename);

        if (owis == null) {
            fail("Cannot find reserved words file " + filename);
            return;
        }

        try {
            Properties props = new Properties();
            props.load(owis);
            Enumeration keys = props.keys();
            while (keys.hasMoreElements()) {
                String word = (String)keys.nextElement();
                databaseReservedWords.add(word);
            }
        } catch (FileNotFoundException e) {
            fail("Unable to load reserved words file " + filename, e);
        } catch (IOException e) {
            fail("Unable to load reserved words file " + filename, e);
        }
    }

    private static Vector<String> databaseReservedWords = new Vector<String>();


}
