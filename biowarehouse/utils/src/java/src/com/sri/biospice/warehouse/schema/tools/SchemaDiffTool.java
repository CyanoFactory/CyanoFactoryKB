/*
 * SchemaDiffTool.java
 *
 * Created on November 15, 2006, 3:52 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package com.sri.biospice.warehouse.schema.tools;

import com.sri.biospice.warehouse.schema.tools.types.Column;
import com.sri.biospice.warehouse.schema.tools.types.Command;
import com.sri.biospice.warehouse.schema.tools.types.Datatypes;
import com.sri.biospice.warehouse.schema.tools.types.Dtype;
import com.sri.biospice.warehouse.schema.tools.types.ForeignKey;
import com.sri.biospice.warehouse.schema.tools.types.Index;
import com.sri.biospice.warehouse.schema.tools.types.Index.Variant;
import com.sri.biospice.warehouse.schema.tools.types.PrimaryKey;
import com.sri.biospice.warehouse.schema.tools.types.Schema;
import com.sri.biospice.warehouse.schema.tools.types.SdKey;
import com.sri.biospice.warehouse.schema.tools.types.Sequence;
import com.sri.biospice.warehouse.schema.tools.types.Table;
import com.sri.biospice.warehouse.schema.tools.types.Tabletypes;
import com.sri.biospice.warehouse.schema.tools.types.Ttype;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;

/**
 *
 * @author gpetry
 */
public class SchemaDiffTool {
    private enum Indent {
        DTYPE ("       "),
        TTYPE ("       "),
        COMMAND ("         "),
        SEQUENCE ("          "),
        TABLE ("       "),
        PRIMARY_KEY ("            "),
        COLUMN ("               "),
        CONSTRAINT ("            "),
        COMMENT ("         "),
        FOREIGN_KEY ("            "),
        SD_KEY ("       "),
        INDEX ("              "),
        VARIANT ("         "),
        ENUMERATION ("        "),
        RESTRICTION ("             ");

        Indent(String value) {
            this.value = value;
        }

        private final String value;

        public String getValue() {
            return value;
        }
    }

    Schema baseSchema = null;
    Schema compSchema = null;

    private static Logger logger = Logger.getLogger(SchemaDiffTool.class.getName());

    StringBuffer changeBuffer = new StringBuffer("");
    
    /**
     * Creates a new instance of SchemaDiffTool
     */
    public SchemaDiffTool(String file1, String file2) throws InstantiationException {
        this(file1, file2, null);
    }

    /**
     * Creates a new instance of SchemaDiffTool
     */
    public SchemaDiffTool(String file1, String file2, String logFile) throws InstantiationException {
        if (logFile != null) {
            try {
                FileHandler fh = new FileHandler(logFile);
                fh.setFormatter(new SimpleFormatter());
                fh.setLevel(Level.INFO);
                logger.setUseParentHandlers(false);
                logger.addHandler(fh);
                ConsoleHandler ch = new ConsoleHandler();
                ch.setLevel(Level.WARNING);
                logger.addHandler(ch);
                System.out.println("Comparison Details will be written to: "+logFile);
            } catch (SecurityException ex) {
                ex.printStackTrace();
                throw new InstantiationException("Security Exception: "+ex.getMessage());
            } catch (IOException ex) {
                ex.printStackTrace();
                throw new InstantiationException("IO Exception: "+ex.getMessage());
            }
        }

        JAXBContext jc;
        try {
            jc = JAXBContext.newInstance("com.sri.biospice.warehouse.schema.tools.types");

            Unmarshaller u = jc.createUnmarshaller();
            baseSchema = (Schema)u.unmarshal(new FileInputStream(file1));
            compSchema = (Schema)u.unmarshal(new FileInputStream(file2));
        } catch (JAXBException ex) {
            ex.printStackTrace();
            throw new InstantiationException("JAXB Error: "+ex.getMessage());
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
            throw new InstantiationException("FileNotFound Error: "+ex.getMessage());
        }

        
        changeBuffer.append("\nBase File:     "+file1+"\n");
        changeBuffer.append("Compared File: "+file2+"\n\n");
        compareSchemas();
        logger.log(Level.INFO, changeBuffer.toString());
    }

    private void compareSchemas() {
        changeBuffer.append("Base File: ");
        if (!baseSchema.getName().equalsIgnoreCase(compSchema.getName())) {
            changeBuffer.append("Schema Name changed from "+baseSchema.getName()+" to "+compSchema.getName());
        }
        if (!baseSchema.getVersion().equalsIgnoreCase(compSchema.getVersion())) {
            changeBuffer.append("Schema Version changed from "+baseSchema.getVersion()+" to "+compSchema.getVersion());
        }

        // Compare Warehouse Types
        // Compare DataTypes
        compareDatatypes();

        // Compare TableTypes
        compareTabletypes();

        //Compare Commands
        compareCommandElements();

        // Compare Sequences
        compareSequenceElements();

        //Compare Tables
        compareTables();
    }

    private Datatypes getDatatypes(Schema schema) {
        List elements = schema.getDatatypesOrTabletypesOrComment();
        for (Object obj: elements) {
            if (obj instanceof Datatypes) {
                return (Datatypes)obj;
            }
        }
        return null;
    }

    private void compareDatatypes() {
        Datatypes base = getDatatypes(baseSchema);
        Datatypes comp = getDatatypes(compSchema);
        List<Dtype> l = base.getDtype();
        List<Dtype> l2 = comp.getDtype();
        StringBuffer buff = new StringBuffer("********* Datatypes *************");
        for (Dtype dtype: l) {
            boolean matches = false;
            String notation = dtype.getNotation();
            String idStr = "\nDtype: notation: \""+notation+"\"";
            String details = "";
            for (Dtype other: l2) {
                if (other.getNotation().equals(notation)) {
                    matches = true;
                    boolean mods = false;
                    if (!dtype.getMysql().equals(other.getMysql())) {
                        mods = true;
                        details = details.concat("\n"+Indent.DTYPE.getValue()+"mysql=\""+dtype.getMysql()+"\" CHANGED TO \""+other.getMysql()+"\"");
                    }
                    if (!dtype.getOracle().equals(other.getOracle())) {
                        mods = true;
                        details = details.concat("\n"+Indent.DTYPE.getValue()+"oracle=\""+dtype.getOracle()+"\" CHANGED TO \""+other.getOracle()+"\"");
                    }
                    if (!dtype.getDescription().equalsIgnoreCase(other.getDescription())) {
                        mods = true;
                        details = details.concat("\n"+Indent.DTYPE.getValue()+"description=\""+dtype.getDescription()+"\" CHANGED TO \""+other.getDescription()+"\"");
                    }
                    l2.remove(other);
                    if (mods) {
                        buff.append(idStr+" MODIFIED "+details+"\n");
                    }
                    break;
                }
            }
            if (!matches) {
                buff.append(idStr+" REMOVED\n");
//                buff.append("\n"+Indent.DTYPE.getValue()+"mysql=\""+dtype.getMysql()+"\"");
//                buff.append("\n"+Indent.DTYPE.getValue()+"oracle=\""+dtype.getOracle()+"\"");
//                buff.append("\n"+Indent.DTYPE.getValue()+"hibernate=\""+dtype.getHibernate()+"\"");
//                buff.append("\n"+Indent.DTYPE.getValue()+"description=\""+dtype.getDescription()+"\"");
            }
        }
        for (Dtype other: l2) {
            buff.append("\nDtype: notation: \""+other.getNotation()+"\" ADDED");
            buff.append("\n"+Indent.DTYPE.getValue()+"mysql=\""+other.getMysql()+"\"");
            buff.append("\n"+Indent.DTYPE.getValue()+"oracle=\""+other.getOracle()+"\"");
            buff.append("\n"+Indent.DTYPE.getValue()+"description=\""+other.getDescription()+"\"\n");
        }

        changeBuffer.append(buff.toString());
    }

    private Tabletypes getTabletypes(Schema schema) {
        List elements = schema.getDatatypesOrTabletypesOrComment();
        for (Object obj: elements) {
            if (obj instanceof Tabletypes) {
                return (Tabletypes)obj;
            }
        }
        return null;
    }

    private void compareTabletypes() {
        int indentLevel = 3;
        Tabletypes base = getTabletypes(baseSchema);
        Tabletypes comp = getTabletypes(compSchema);
        List<Ttype> l = base.getTtype();
        List<Ttype> l2 = comp.getTtype();
        StringBuffer buff = new StringBuffer("\n********* Tabletypes *************");
        for (Ttype ttype: l) {
            boolean matches = false;
            String name = ttype.getName();
            String details = "\nTtype: name: \""+name+"\"";
            for (Ttype other: l2) {
                if (other.getName().equals(name)) {
                    matches = true;
                    boolean mods = false;
                    String s = compareComments(Indent.TTYPE, ttype.getComment(), other.getComment());
                    if (s.length() != 0) {
                        buff.append(details+" MODIFIED\n"+s);
                    }
                    l2.remove(other);
                    break;
                }
            }
            if (!matches) {
                buff.append(details+" REMOVED\n");
//                for (String comment: ttype.getComment()) {
//                    buff.append(Indent.TTYPE.getValue()+"comment=\""+comment+"\"\n");
//                }
            }
        }
        for (Ttype other: l2) {
            buff.append("\nTtype: name: \""+other.getName()+"\" ADDED");
            for (String comment: other.getComment()) {
                buff.append("\n"+Indent.TTYPE.getValue()+"comment=\""+comment+"\"\n");
            }
        }

        changeBuffer.append(buff.toString());
    }

    private List<Command> getCommands(Schema schema) {
        ArrayList<Command> elements = new ArrayList<Command>();
        List schemaElements = schema.getDatatypesOrTabletypesOrComment();
        for (Object obj: schemaElements) {
            if (obj instanceof Command) {
                elements.add((Command)obj);
            }
        }
        return elements;
    }

    private void compareCommandElements() {
        List<Command> base = getCommands(baseSchema);
        List<Command> comp = getCommands(compSchema);
        ArrayList<Command> temp = new ArrayList<Command>();
        StringBuffer buff = new StringBuffer("\n********* Commands *************");

        for (Command command: base) {
            boolean matches = false;
            String value = command.getValue();
            String dialect = command.getDialect();
            String details = "\nCommand: \""+value+"\"";
            for (Command other: comp) {
                if (value.equalsIgnoreCase(other.getValue())) {
                    matches = true;
                    temp.add(other);
                    String otherdialect = other.getDialect();
                    if ((dialect == null) && (otherdialect == null)) {
                        // No Change
                    } else if (dialect == null) {
                        buff.append(details+" MODIFIED\n");
                        buff.append(Indent.COMMAND.getValue()+"dialect=\""+otherdialect+"\" ADDED\n");
                    } else if (otherdialect == null) {
                        buff.append(details+" MODIFIED\n");
                        buff.append(Indent.COMMAND.getValue()+"dialect=\""+dialect+"\" REMOVED\n");
                    } else {
                        if (!dialect.equals(otherdialect)) {
                            buff.append(details+" MODIFIED\n");
                            buff.append(Indent.COMMAND.getValue()+"dialect=\""+dialect+"\" CHANGED TO \"\""+otherdialect+"\"\n");
                        }
                    }
                }
            }

            if (!matches) {
                buff.append(details);
                if (dialect != null) {
                    buff.append(" dialect=\""+dialect+"\"");
                }
                buff.append(" REMOVED\n");
            }
        }

        for (Command other: comp) {
            boolean matches = false;
            for (Command command: temp) {
                if (other.getValue().equalsIgnoreCase(command.getValue())) {
                    matches = true;
                }
            }

            if (!matches) {
                buff.append("\nCommand: \""+other.getValue()+"\""+((other.getDialect() != null) ? " dialect=\""+other.getDialect()+"\"" : "")+" ADDED\n");
            }
        }

        changeBuffer.append(buff.toString());
    }

    private List<Sequence> getSequences(Schema schema) {
        ArrayList<Sequence> elements = new ArrayList<Sequence>();
        List schemaElements = schema.getDatatypesOrTabletypesOrComment();
        for (Object obj: schemaElements) {
            if (obj instanceof Sequence) {
                elements.add((Sequence)obj);
            }
        }
        return elements;
    }

    private void compareSequenceElements() {
        int indentLevel = 5;
        List<Sequence> base = getSequences(baseSchema);
        List<Sequence> comp = getSequences(compSchema);
        ArrayList<Sequence> temp = new ArrayList<Sequence>();
        StringBuffer buff = new StringBuffer("\n********* Sequences *************");

        for (Sequence sequence: base) {
            boolean matches = false;
            String name = sequence.getName();
            String start = sequence.getStart();
            String minValue = sequence.getMinValue();
            String increment = sequence.getIncrement();
            String idStr = "\nSequence: name: \""+name+"\"";
            String details = "";
            for (Sequence other: comp) {
                if (name.equalsIgnoreCase(other.getName())) {
                    matches = true;
                    temp.add(other);
                    boolean mod = false;
                    if (!start.equals(other.getStart())) {
                        mod = true;
                        details = details.concat("\n"+Indent.SEQUENCE.getValue()+"start=\""+start+"\" CHANGED TO \""+other.getStart()+"\"");
                    }
                    if (!minValue.equals(other.getMinValue())) {
                        mod = true;
                        details = details.concat("\n"+Indent.SEQUENCE.getValue()+"minValue=\""+minValue+"\" CHANGED TO \""+other.getMinValue()+"\"");
                    }
                    if (!increment.equals(other.getIncrement())) {
                        mod = true;
                        details = details.concat("\n"+Indent.SEQUENCE.getValue()+"increment=\""+increment+"\" CHANGED TO \""+other.getIncrement()+"\"");
                    }
                    String s = compareComments(Indent.SEQUENCE, sequence.getComment(), other.getComment());
                    if (s.length() > 0) {
                        mod = true;
                        details = details.concat("\n"+s);
                    }
                    if (mod) {
                        buff.append(idStr+" MODIFIED "+details);
                    }
                }
            }

            if (!matches) {
                buff.append(idStr+" REMOVED\n");
            }
        }

        for (Sequence other: comp) {
            boolean matches = false;
            for (Sequence sequence: temp) {
                if (other.getName().equalsIgnoreCase(sequence.getName())) {
                    matches = true;
                }
            }

            if (!matches) {
                buff.append("\nSequence: name: \""+other.getName()+"\" ADDED");
                buff.append("\n"+Indent.SEQUENCE.getValue()+"start=\""+other.getStart()+"\"");
                buff.append("\n"+Indent.SEQUENCE.getValue()+"minValue=\""+other.getMinValue()+"\"");
                buff.append("\n"+Indent.SEQUENCE.getValue()+"increment=\""+other.getIncrement()+"\"\n");
            }
        }

        changeBuffer.append(buff.toString());
    }

    private List<Table> getTables(Schema schema) {
        ArrayList<Table> elements = new ArrayList<Table>();
        List schemaElements = schema.getDatatypesOrTabletypesOrComment();
        for (Object obj: schemaElements) {
            if (obj instanceof Table) {
                elements.add((Table)obj);
            }
        }
        return elements;
    }

    private void compareTables() {
        int indentLevel = 3;
        List<Table> base = getTables(baseSchema);
        List<Table> comp = getTables(compSchema);
        ArrayList<Table> temp = new ArrayList<Table>();
        StringBuffer buff = new StringBuffer("\n********* Tables *************");

        for (Table table: base) {
            boolean matches = false;
            String name = table.getName();
            String type = table.getType();
            PrimaryKey pKey = table.getPrimaryKey();
            String idStr = "\nTable: name: \""+name+"\"";
            String details = "";
            for (Table other: comp) {
                if (name.equalsIgnoreCase(other.getName())) {
                    matches = true;
                    temp.add(other);
                    boolean mod = false;
                    if (!type.equals(other.getType())) {
                        mod = true;
                        details = details.concat("\n"+Indent.TABLE.getValue()+"type=\""+type+"\" CHANGED TO \""+other.getType()+"\"\n");
                    }

                    String pk = comparePrimaryKeys(Indent.TABLE, table.getPrimaryKey(), other.getPrimaryKey());
                    if (pk.length() > 0) {
                        //Good
                        mod = true;
                        details = details.concat("\n"+pk);
                    }

                    String col = compareColumns(table.getColumn(), other.getColumn());
                    if (col.length() > 0) {
                        mod = true;
                        details = details.concat("\n"+col);
                    }

                    String c = compareConstraints(Indent.TABLE, table.getConstraint(), other.getConstraint());
                    if (c.length() > 0) {
                        //Good
                        mod = true;
                        details = details.concat(c);
                    }

                    String fk = compareForeignKeys(Indent.TABLE, table.getForeignKey(), other.getForeignKey());
                    if (fk.length() > 0) {
                        //Good
                        mod = true;
                        details = details.concat("\n"+fk);
                    }

                    String idx = compareIndexes(table.getIndex(), other.getIndex());
                    if (idx.length() > 0) {
                        //Good
                        mod = true;
                        details = details.concat("\n"+idx);
                    }

                    String s = compareComments(Indent.TABLE, table.getComment(), other.getComment());
                    if (s.length() > 0) {
                        //Good
                        mod = true;
                        details = details.concat(s);
                    }
                    if (mod) {
                        buff.append(idStr+" MODIFIED "+details);
                    }
                }
            }

            if (!matches) {
                buff.append(idStr+" REMOVED\n");
            }
        }

        for (Table other: comp) {
            boolean matches = false;
            for (Table table: temp) {
                if (other.getName().equalsIgnoreCase(table.getName())) {
                    matches = true;
                }
            }

            if (!matches) {
                // New table
                buff.append("\nTable: name: \""+other.getName()+"\" ADDED");
                if (other.getPrimaryKey() != null) {
                    buff.append("\n"+Indent.TABLE.getValue()+"PrimaryKey: name=\""+other.getPrimaryKey().getName()+"\"");
                    buff.append("\n"+Indent.TABLE.getValue()+Indent.PRIMARY_KEY.getValue()+"column=\""+other.getPrimaryKey().getColumn()+"\"\n");
                }
                for (Column column: other.getColumn()) {
                    buff.append("\n"+Indent.TABLE.getValue()+"Column name=\""+column.getName()+"\"n");
                    buff.append("\n"+Indent.COLUMN.getValue()+"type=\""+column.getType()+"\"");
                    buff.append("\n"+Indent.COLUMN.getValue()+"required=\""+column.getRequired()+"\"\n");
                    if (column.getLength() != null) {
                        buff.append(Indent.COLUMN.getValue()+"length=\""+column.getLength()+"\"\n");
                    }
                    if (column.getPrecision() != null) {
                        buff.append(Indent.COLUMN.getValue()+"precision=\""+column.getPrecision()+"\"\n");
                    }
                    if (column.getScale() != null) {
                        buff.append(Indent.COLUMN.getValue()+"scale=\""+column.getScale()+"\"\n");
                    }
                    if (column.getAutoIncrement() != null) {
                        buff.append(Indent.COLUMN.getValue()+"auto-increment=\""+column.getAutoIncrement()+"\"\n");
                    }
                    if (column.getSize() != null) {
                        buff.append(Indent.COLUMN.getValue()+"size=\""+column.getSize()+"\"\n");
                    }

                    for (Object obj: column.getContent()) {
                        if (obj instanceof ForeignKey) {
                            ForeignKey fkey = (ForeignKey) obj;
                            buff.append(Indent.COLUMN.getValue()+"ForeignKey: name=\""+fkey.getName()+"\"\n");
                            if (fkey.getColumn() != null) {
                                buff.append(Indent.COLUMN.getValue()+Indent.FOREIGN_KEY.getValue()+"column=\""+fkey.getColumn()+"\"\n");
                            }
                            buff.append(Indent.COLUMN.getValue()+Indent.FOREIGN_KEY.getValue()+"totable=\""+fkey.getToTable()+"\"\n");
                            buff.append(Indent.COLUMN.getValue()+Indent.FOREIGN_KEY.getValue()+"tocolumn=\""+fkey.getToColumn()+"\"\n");
                        } else if (obj instanceof SdKey) {
                            SdKey skey = (SdKey) obj;
                            buff.append(Indent.COLUMN.getValue()+"SdKey: name=\""+skey.getName()+"\"\n");
                            buff.append(Indent.COLUMN.getValue()+Indent.SD_KEY.getValue()+"totable=\""+skey.getToTable()+"\"\n");
                            buff.append(Indent.COLUMN.getValue()+Indent.SD_KEY.getValue()+"tocolumn=\""+skey.getToColumn()+"\"\n");
                        } else if (obj instanceof com.sri.biospice.warehouse.schema.tools.types.Enumeration) {
                            buff.append(Indent.COLUMN.getValue()+"Enumerations:\n");
                            com.sri.biospice.warehouse.schema.tools.types.Enumeration ens = (com.sri.biospice.warehouse.schema.tools.types.Enumeration) obj;
                            for (com.sri.biospice.warehouse.schema.tools.types.Enum e: ens.getRestriction()) {
                                buff.append(Indent.COLUMN.getValue()+Indent.ENUMERATION.getValue()+"Restriction: value=\""+e.getValue()+"\"\n");
                                buff.append(Indent.COLUMN.getValue()+Indent.ENUMERATION.getValue()+Indent.RESTRICTION.getValue()+"description=\""+e.getDescription()+"\"\n");
                            }
                        } else if (obj instanceof JAXBElement) {
                            JAXBElement jobj = (JAXBElement) obj;
                            if (jobj.getDeclaredType().equals(String.class)) {
                                if (jobj.getName().getLocalPart().equalsIgnoreCase("comment")) {
                                    String str = (String)jobj.getValue();
                                    buff.append(Indent.COLUMN.getValue()+"Comment: \""+str.trim()+"\"\n");
                                }
                            } else if (jobj.getDeclaredType().equals(com.sri.biospice.warehouse.schema.tools.types.Enumeration.class)) {
                                buff.append(Indent.COLUMN.getValue()+"Enumerations:\n");
                                com.sri.biospice.warehouse.schema.tools.types.Enumeration ens = (com.sri.biospice.warehouse.schema.tools.types.Enumeration) jobj.getValue();
                                for (com.sri.biospice.warehouse.schema.tools.types.Enum e: ens.getRestriction()) {
                                    buff.append(Indent.COLUMN.getValue()+Indent.ENUMERATION.getValue()+"Restriction: value=\""+e.getValue()+"\"\n");
                                    buff.append(Indent.COLUMN.getValue()+Indent.ENUMERATION.getValue()+Indent.RESTRICTION.getValue()+"description=\""+e.getDescription()+"\"\n");
                                }
                            } else {
                                logger.log(Level.SEVERE, "UNKNOWN: TYPE"+jobj.getDeclaredType().getName());
                                logger.log(Level.SEVERE, "UNKNOWN: NME"+jobj.getName());
                                logger.log(Level.SEVERE, "UNKNOWN: VAL"+jobj.getValue());
                            }
                        } else {
                            if (obj instanceof String) {
                                String str = (String)obj;
                                if (str.trim().length() > 0) {
                                    buff.append(Indent.COLUMN.getValue()+"Comment: \""+str.trim()+"\"\n");
                                }
                            } else {
                                logger.log(Level.SEVERE, "  TA UNKNOWN ["+obj.getClass().getName()+"] -> "+obj.toString());
                            }
                        }
                    }
                }
                //Constraints
                for (String constraint: other.getConstraint()) {
                    buff.append(Indent.TABLE.getValue()+"Constraint: \""+constraint+"\"\n");
                }

                //Foreign Key
                for (ForeignKey fkey: other.getForeignKey()) {
                    buff.append(Indent.TABLE.getValue()+"ForeignKey: name=\""+fkey.getName()+"\"\n");
                    if (fkey.getColumn() != null) {
                        buff.append(Indent.TABLE.getValue()+Indent.FOREIGN_KEY.getValue()+"column=\""+fkey.getColumn()+"\"\n");
                    }
                    buff.append(Indent.TABLE.getValue()+Indent.FOREIGN_KEY.getValue()+"totable=\""+fkey.getToTable()+"\"\n");
                    buff.append(Indent.TABLE.getValue()+Indent.FOREIGN_KEY.getValue()+"tocolumn=\""+fkey.getToColumn()+"\"\n");
                }

                //Index
                for (Index idx: other.getIndex()) {
                    buff.append(Indent.TABLE.getValue()+"Index: name=\""+idx.getName()+"\"\n");
                    if (idx.getColumns() != null) {
                        buff.append(Indent.INDEX.getValue()+"columns=\""+idx.getColumns()+"\"\n");
                    }
                    if (idx.getColumns() != null) {
                        buff.append(Indent.INDEX.getValue()+"initial=\""+idx.getInitial()+"\"\n");
                    }
                    for (Variant v: idx.getVariant()) {
                        buff.append(Indent.INDEX.getValue()+"Variant: dialect=\""+v.getDialect()+"\" columns=\""+v.getColumns()+"\"\n");
                    }
                }

                //Comment
                for (String comment: other.getComment()) {
                    buff.append(Indent.TABLE.getValue()+"Comment: \""+comment+"\"\n");
                }
            }
        }

        changeBuffer.append(buff.toString());
    }

    private String comparePrimaryKeys(Indent indent, PrimaryKey base, PrimaryKey other) {
        StringBuffer buff = new StringBuffer("");
        if (base == null) {
            if (other != null) {
                buff.append(indent.getValue()+"PrimaryKey: name=\""+other.getName()+"\" column=\""+other.getColumn()+"\" ADDED");
            }
        } else {
            if (other != null) {
                String name = base.getName();
                String column = base.getColumn();
                // Ignore Content - Primary Keys normally do not have content
                if (name.equalsIgnoreCase(other.getName())) {
                    if (!column.equalsIgnoreCase(other.getColumn())) {
                        buff.append(indent.getValue()+"PrimaryKey: name=\""+other.getName()+"\" MODIFIED\n");
                        buff.append(indent.getValue()+Indent.PRIMARY_KEY.getValue()+"column=\""+column+"\" CHANGED TO \""+other.getColumn()+"\"");
                    }
                } else {
                    buff.append(indent.getValue()+"PrimaryKey: REPLACED\n");
                    buff.append(indent.getValue()+"    OLD name=\""+name+"\"\n");
                    buff.append(indent.getValue()+Indent.PRIMARY_KEY.getValue()+"column=\""+column+"\"\n");
                    buff.append(indent.getValue()+"    NEW name=\""+other.getName()+"\"\n");
                    buff.append(indent.getValue()+Indent.PRIMARY_KEY.getValue()+"column=\""+other.getColumn()+"\"");
                }
            } else {
                buff.append(indent.getValue()+"PrimaryKey: name=\""+base.getName()+"\" column=\""+base.getColumn()+"\" REMOVED");
            }
        }

        return buff.toString();
    }

    private String compareColumns(List<Column> base, List<Column> comp) {
        StringBuffer buff = new StringBuffer("");
        ArrayList<Column> temp = new ArrayList<Column>();
        for (Column column: base) {
            boolean matches = false;
            String name = column.getName();
            String type = column.getType();
            String req = column.getRequired();
            String length = column.getLength();
            String auto = column.getAutoIncrement();
            String scale = column.getScale();
            String precision = column.getPrecision();
            String size = column.getSize();

            String idStr = Indent.TABLE.getValue()+"Column: name=\""+name+"\"";
            String details = "";
            for (Column other: comp) {
                if (name.equalsIgnoreCase(other.getName())) {
                    matches = true;
                    temp.add(other);
                    boolean mod = false;
                    if (!type.equals(other.getType())) {
                        mod = true;
                        details = details.concat("\n"+Indent.COLUMN.getValue()+"type=\""+type+"\" CHANGED TO \""+other.getType()+"\"\n");
                    }

                    if (!req.equals(other.getRequired())) {
                        mod = true;
                        details = details.concat("\n"+Indent.COLUMN.getValue()+"required=\""+req+"\" CHANGED TO \""+other.getRequired()+"\"\n");
                    }

                    if (length == null) {
                        if (other.getLength() != null) {
                            mod = true;
                            details = details.concat("\n"+Indent.COLUMN.getValue()+"length=\""+other.getLength()+"\" ADDED\n");
                        }
                    } else {
                        if (other.getLength() == null) {
                            mod = true;
                            details = details.concat("\n"+Indent.COLUMN.getValue()+"length=\""+length+"\" REMOVED\n");
                        } else {
                            if (!length.equals(other.getLength())) {
                                mod = true;
                                details = details.concat("\n"+Indent.COLUMN.getValue()+"length=\""+length+"\" CHANGED TO \""+other.getLength()+"\"\n");
                            }
                        }
                    }

                    if (precision == null) {
                        if (other.getPrecision() != null) {
                            mod = true;
                            details = details.concat("\n"+Indent.COLUMN.getValue()+"precision=\""+other.getPrecision()+"\" ADDED\n");
                        }
                    } else {
                        if (other.getPrecision() == null) {
                            mod = true;
                            details = details.concat("\n"+Indent.COLUMN.getValue()+"precision=\""+precision+"\" REMOVED\n");
                        } else {
                            if (!precision.equals(other.getPrecision())) {
                                mod = true;
                                details = details.concat("\n"+Indent.COLUMN.getValue()+"precision=\""+precision+"\" CHANGED TO \""+other.getPrecision()+"\"\n");
                            }
                        }
                    }

                    if (scale == null) {
                        if (other.getScale() != null) {
                            mod = true;
                            details = details.concat("\n"+Indent.COLUMN.getValue()+"scale=\""+other.getScale()+"\" ADDED\n");
                        }
                    } else {
                        if (other.getScale() == null) {
                            mod = true;
                            details = details.concat("\n"+Indent.COLUMN.getValue()+"scale=\""+scale+"\" REMOVED\n");
                        } else {
                            if (!scale.equals(other.getScale())) {
                                mod = true;
                                details = details.concat("\n"+Indent.COLUMN.getValue()+"scale=\""+scale+"\" CHANGED TO \""+other.getScale()+"\"\n");
                            }
                        }
                    }

                    if (auto == null) {
                        if (other.getAutoIncrement() != null) {
                            mod = true;
                            details = details.concat("\n"+Indent.COLUMN.getValue()+"auto-increment=\""+other.getAutoIncrement()+"\" ADDED\n");
                        }
                    } else {
                        if (other.getAutoIncrement() == null) {
                            mod = true;
                            details = details.concat("\n"+Indent.COLUMN.getValue()+"auto-increment=\""+auto+"\" REMOVED\n");
                        } else {
                            if (!auto.equals(other.getAutoIncrement())) {
                                mod = true;
                                details = details.concat("\n"+Indent.COLUMN.getValue()+"auto-increment=\""+auto+"\" CHANGED TO \""+other.getAutoIncrement()+"\"\n");
                            }
                        }
                    }

                    if (size == null) {
                        if (other.getSize() != null) {
                            mod = true;
                            details = details.concat("\n"+Indent.COLUMN.getValue()+"size=\""+other.getSize()+"\" ADDED\n");
                        }
                    } else {
                        if (other.getSize() == null) {
                            mod = true;
                            details = details.concat("\n"+Indent.COLUMN.getValue()+"size=\""+size+"\" REMOVED\n");
                        } else {
                            if (!size.equals(other.getSize())) {
                                mod = true;
                                details = details.concat("\n"+Indent.COLUMN.getValue()+"size=\""+size+"\" CHANGED TO \""+other.getSize()+"\"\n");
                            }
                        }
                    }

                    String fk = compareForeignKeys(Indent.COLUMN, getForeignKeys(column.getContent()), getForeignKeys(other.getContent()));
                    if (fk.length() > 0) {
                        //Good
                        mod = true;
                        details = details.concat(fk);
                    }

                    String sk = compareSdKeys(Indent.COLUMN, getSdKeys(column.getContent()), getSdKeys(other.getContent()));
                    if (sk.length() > 0) {
                        //Good
                        mod = true;
                        details = details.concat(sk);
                    }

                    String e = compareEnumerations(Indent.COLUMN, getEnumeration(column.getContent()), getEnumeration(other.getContent()));
                    if (e.length() > 0) {
                        mod = true;
                        details = details.concat("\n\t        Enumerations:  MODIFIED\n"+e);
                    }

                    String c = compareComments(Indent.COLUMN, getComments(column.getContent()), getComments(other.getContent()));
                    if (c.length() > 0) {
                        mod = true;
                        details = details.concat(c);
                    }
                    if (mod) {
                        buff.append(idStr+" MODIFIED "+details);
                    }
                }
            }

            if (!matches) {
                buff.append(idStr+" REMOVED\n");
            }
        }

        for (Column other: comp) {
            boolean matches = false;
            for (Column key: temp) {
                if (other.getName().equalsIgnoreCase(key.getName())) {
                    matches = true;
                }
            }

            if (!matches) {
                buff.append(Indent.TABLE.getValue()+"Column: name=\""+other.getName()+"\" ADDED\n");
                buff.append(Indent.COLUMN.getValue()+"type=\""+other.getType()+"\"\n");
                buff.append(Indent.COLUMN.getValue()+"required=\""+other.getRequired()+"\"\n");
                if (other.getLength() != null) {
                    buff.append(Indent.COLUMN.getValue()+"length=\""+other.getLength()+"\"\n");
                }
                if (other.getScale() != null) {
                    buff.append(Indent.COLUMN.getValue()+"scale=\""+other.getScale()+"\"\n");
                }
                if (other.getPrecision() != null) {
                    buff.append(Indent.COLUMN.getValue()+"precision=\""+other.getPrecision()+"\"\n");
                }
                if (other.getAutoIncrement() != null) {
                    buff.append(Indent.COLUMN.getValue()+"auto-increment=\""+other.getAutoIncrement()+"\"\n");
                }
                if (other.getSize() != null) {
                    buff.append(Indent.COLUMN.getValue()+"size=\""+other.getSize()+"\"\n");
                }

                for (Object obj :other.getContent()) {
                    if (obj instanceof ForeignKey) {
                        ForeignKey fkey = (ForeignKey) obj;
                        buff.append(Indent.COLUMN.getValue()+"ForeignKey name=\""+fkey.getName()+"\"\n");
                        if (fkey.getColumn() != null) {
                            buff.append(Indent.COLUMN.getValue()+Indent.FOREIGN_KEY.getValue()+"column=\""+fkey.getColumn()+"\"\n");
                        }
                        buff.append(Indent.COLUMN.getValue()+Indent.FOREIGN_KEY.getValue()+"totable=\""+fkey.getToTable()+"\"\n");
                        buff.append(Indent.COLUMN.getValue()+Indent.FOREIGN_KEY.getValue()+"tocolumn=\""+fkey.getToColumn()+"\"\n");
                    } else if (obj instanceof SdKey) {
                        SdKey skey = (SdKey) obj;
                        buff.append(Indent.COLUMN.getValue()+"SdKey name=\""+skey.getName()+"\"\n");
                        buff.append(Indent.COLUMN.getValue()+Indent.SD_KEY.getValue()+"totable=\""+skey.getToTable()+"\"\n");
                        buff.append(Indent.COLUMN.getValue()+Indent.SD_KEY.getValue()+"tocolumn=\""+skey.getToColumn()+"\"\n");
                    } else if (obj instanceof com.sri.biospice.warehouse.schema.tools.types.Enumeration) {
                        buff.append(Indent.COLUMN.getValue()+"Enumerations:\n");
                        com.sri.biospice.warehouse.schema.tools.types.Enumeration ens = (com.sri.biospice.warehouse.schema.tools.types.Enumeration) obj;
                        for (com.sri.biospice.warehouse.schema.tools.types.Enum e: ens.getRestriction()) {
                            buff.append(Indent.COLUMN.getValue()+Indent.ENUMERATION.getValue()+"Restriction: value=\""+e.getValue()+"\"\n");
                            buff.append(Indent.COLUMN.getValue()+Indent.ENUMERATION.getValue()+Indent.RESTRICTION.getValue()+"description=\""+e.getDescription()+"\"\n");
                        }
                    } else if (obj instanceof JAXBElement) {
                        JAXBElement jobj = (JAXBElement) obj;
                        if (jobj.getDeclaredType().equals(String.class)) {
                            if (jobj.getName().getLocalPart().equalsIgnoreCase("comment")) {
                                String str = (String)jobj.getValue();
                                buff.append(Indent.COLUMN.getValue()+"Comment: \""+str.trim()+"\"\n");
                            }
                        } else if (jobj.getDeclaredType().equals(com.sri.biospice.warehouse.schema.tools.types.Enumeration.class)) {
                        buff.append(Indent.COLUMN.getValue()+"Enumerations:\n");
                            com.sri.biospice.warehouse.schema.tools.types.Enumeration ens = (com.sri.biospice.warehouse.schema.tools.types.Enumeration) jobj.getValue();
                            for (com.sri.biospice.warehouse.schema.tools.types.Enum e: ens.getRestriction()) {
                                buff.append(Indent.COLUMN.getValue()+Indent.ENUMERATION.getValue()+"Restriction: value=\""+e.getValue()+"\"\n");
                                buff.append(Indent.COLUMN.getValue()+Indent.ENUMERATION.getValue()+Indent.RESTRICTION.getValue()+"description=\""+e.getDescription()+"\"\n");
                            }
                        }
                    } else {
                        //Skip
                        if (obj instanceof String) {
                            String str = (String)obj;
                            if (str.trim().length() > 0) {
                                buff.append(Indent.COLUMN.getValue()+"Comment: \""+str.trim()+"\"\n");
                            }
                        } else {
                            logger.log(Level.SEVERE, "\t   Cl Add   UNKNOWN Element ["+obj.getClass().getName()+"] "+obj.toString());
                        }
                    }
                }
            }
        }
        return buff.toString();
    }

    private String compareConstraints(Indent indent, List<String> base, List<String> comp) {
        return compareStringLists(indent, "Constraint", base, comp);
    }

    private List<String> getComments(List list) {
        ArrayList<String> tmp = new ArrayList<String>();
        for (Object obj: list) {

            if (obj instanceof JAXBElement) {
                JAXBElement jobj = (JAXBElement) obj;
                if (jobj.getDeclaredType().equals(String.class)) {
                    String jos = jobj.getName().getLocalPart();
                    if (jos.equalsIgnoreCase("comment")) {
                        String str = (String)jobj.getValue();
                        tmp.add(str.trim());
                    }
                }
            } else if (obj instanceof String) {
                String str = (String)obj;
                if (str.trim().length() > 0) {
                    tmp.add(str.trim());
                }
            }
        }

        return tmp;
    }

    private String compareComments(Indent indent, List<String> base, List<String> comp) {
        return compareStringLists(indent, "Comment", base, comp);
    }

    private List<ForeignKey> getForeignKeys(List list) {
        ArrayList<ForeignKey> tmp = new ArrayList<ForeignKey>();
        for (Object obj: list) {
            if (obj instanceof ForeignKey) {
                tmp.add((ForeignKey)obj);
            }
        }

        return tmp;
    }

    private String compareForeignKeys(Indent indent, List<ForeignKey> base, List<ForeignKey> comp) {
        ArrayList<ForeignKey> temp = new ArrayList<ForeignKey>();
        StringBuffer buff = new StringBuffer("");
        if (base != null) {
            if (comp != null) {
                for (ForeignKey fKey: base) {
                    boolean matches = false;
                    String name = fKey.getName();
                    String column = fKey.getColumn();
                    String tocolumn = fKey.getToColumn();
                    String totable = fKey.getToTable();
                    String idStr = "\n"+indent.getValue()+"ForeignKey: name: \""+name+"\"";
                    String details = "";
                    for (ForeignKey other: comp) {
                        if (name.equalsIgnoreCase(other.getName())) {
                            matches = true;
                            boolean mod = false;
                            temp.add(other);
                            if (column != null) {
                                if (other.getColumn() != null) {
                                    if (!column.equals(other.getColumn())) {
                                        mod = true;
                                        details = details.concat("\n"+indent.getValue()+Indent.FOREIGN_KEY.getValue()+"column=\""+column+"\" CHANGED TO \""+other.getColumn()+"\"\n");
                                    }
                                } else {
                                    mod = true;
                                    details = details.concat("\n"+indent.getValue()+Indent.FOREIGN_KEY.getValue()+"column=\""+column+"\" REMOVED\n");
                                }
                            } else {
                                if (other.getColumn() != null) {
                                    mod = true;
                                    details = details.concat("\n"+indent.getValue()+Indent.FOREIGN_KEY.getValue()+"column=\""+other.getColumn()+"\" ADDED\n");
                                }
                            }
                            if (!totable.equals(other.getToTable())) {
                                mod = true;
                                details = details.concat("\n"+indent.getValue()+Indent.FOREIGN_KEY.getValue()+"totable=\""+totable+"\" CHANGED TO \""+other.getToTable()+"\"\n");
                            }
                            if (!tocolumn.equals(other.getToColumn())) {
                                mod = true;
                                details = details.concat("\n"+indent.getValue()+Indent.FOREIGN_KEY.getValue()+"tocolumn=\""+tocolumn+"\" CHANGED TO \""+other.getToColumn()+"\"\n");
                            }
                            if (mod) {
                                buff.append(idStr+" MODIFIED "+details);
                            }
                        }
                    }

                    if (!matches) {
                        buff.append(idStr+" REMOVED\n");
                    }
                }

                for (ForeignKey other: comp) {
                    boolean matches = false;
                    for (ForeignKey key: temp) {
                        if (other.getName().equalsIgnoreCase(key.getName())) {
                            matches = true;
                        }
                    }

                    if (!matches) {
                        buff.append("\n"+indent.getValue()+"ForeignKey: name: \""+other.getName()+"\" ADDED");
                        if (other.getColumn() != null) {
                            buff.append("\n"+indent.getValue()+Indent.FOREIGN_KEY.getValue()+"column=\""+other.getColumn()+"\"");
                        }
                        buff.append("\n"+indent.getValue()+Indent.FOREIGN_KEY.getValue()+"totable=\""+other.getToTable()+"\"");
                        buff.append("\n"+indent.getValue()+Indent.FOREIGN_KEY.getValue()+"tocolumn=\""+other.getToColumn()+"\"\n");
                    }
                }
            } else {
                for (ForeignKey key: base) {
                    buff.append(indent.getValue()+"ForeignKey: name: \""+key.getName()+"\" REMOVED\n");
                }
            }
        } else {
            if (comp != null) {
                for (ForeignKey other: comp) {
                    buff.append("\n"+indent.getValue()+"ForeignKey: name: \""+other.getName()+"\" ADDED");
                    if (other.getColumn() != null) {
                        buff.append("\n"+indent.getValue()+Indent.FOREIGN_KEY.getValue()+"column=\""+other.getColumn()+"\"");
                    }
                    buff.append("\n"+indent.getValue()+Indent.FOREIGN_KEY.getValue()+"totable=\""+other.getToTable()+"\"");
                    buff.append("\n"+indent.getValue()+Indent.FOREIGN_KEY.getValue()+"tocolumn=\""+other.getToColumn()+"\"\n");
                }
            }
        }

        return buff.toString();
    }

    private List<SdKey> getSdKeys(List list) {
        ArrayList<SdKey> tmp = new ArrayList<SdKey>();
        for (Object obj: list) {
            if (obj instanceof SdKey) {
                tmp.add((SdKey)obj);
            }
        }

        return tmp;
    }

    private String compareSdKeys(Indent indent, List<SdKey> base, List<SdKey> comp) {
        ArrayList<SdKey> temp = new ArrayList<SdKey>();
        StringBuffer buff = new StringBuffer("");
        if (base != null) {
            if (comp != null) {
                for (SdKey key: base) {
                    boolean matches = false;
                    String name = key.getName();
                    String tocolumn = key.getToColumn();
                    String totable = key.getToTable();
                    String idStr = "\n"+indent.getValue()+"SdKey: name: \""+name+"\"";
                    String details = "";
                    for (SdKey other: comp) {
                        if (name.equalsIgnoreCase(other.getName())) {
                            matches = true;
                            temp.add(other);
                            boolean mod = false;
                            if (!totable.equals(other.getToTable())) {
                                mod = true;
                                details = details.concat(indent.getValue()+Indent.SD_KEY.getValue()+"totable=\""+totable+"\" CHANGED TO \""+other.getToTable()+"\"\n");
                            }
                            if (!tocolumn.equals(other.getToColumn())) {
                                mod = true;
                                details = details.concat(indent.getValue()+Indent.SD_KEY.getValue()+"tocolumn=\""+tocolumn+"\" CHANGED TO \""+other.getToColumn()+"\"\n");
                            }
                            if (mod) {
                                buff.append(idStr+" MODIFIED\n"+details);
                            }
                        }
                    }

                    if (!matches) {
                        buff.append(idStr+" REMOVED\n");
                    }
                }

                for (SdKey other: comp) {
                    boolean matches = false;
                    for (SdKey key: temp) {
                        if (other.getName().equalsIgnoreCase(key.getName())) {
                            matches = true;
                        }
                    }

                    if (!matches) {
                        buff.append("\n"+indent.getValue()+"SdKey: name: \""+other.getName()+"\" ADDED\n");
                        buff.append(indent.getValue()+Indent.SD_KEY.getValue()+"totable=\""+other.getToTable()+"\"\n");
                        buff.append(indent.getValue()+Indent.SD_KEY.getValue()+"tocolumn=\""+other.getToColumn()+"\"\n");
                    }
                }
            } else {
                for (SdKey key: base) {
                    buff.append("\n"+indent.getValue()+"SdKey: name: \""+key.getName()+"\" REMOVED\n");
                }
            }
        } else {
            if (comp != null) {
                for (SdKey other: comp) {
                    buff.append("\n"+indent.getValue()+"SdKey: name: \""+other.getName()+"\" ADDED\n");
                    buff.append(indent.getValue()+Indent.SD_KEY.getValue()+"totable=\""+other.getToTable()+"\"\n");
                    buff.append(indent.getValue()+Indent.SD_KEY.getValue()+"tocolumn=\""+other.getToColumn()+"\"\n");
                }
            }
        }
        return buff.toString();
    }

    private String compareIndexes(List<Index> base, List<Index> comp) {
        StringBuffer buff = new StringBuffer("");
        ArrayList<Index> temp = new ArrayList<Index>();
        for (Index index: base) {
            boolean matches = false;
            String name = index.getName();
            String columns = index.getColumns();
            String initial = index.getInitial();
            String idStr = Indent.TABLE.getValue()+"Index: name=\""+name+"\"";
            String details = "";
            for (Index other: comp) {
                if (name.equalsIgnoreCase(other.getName())) {
                    matches = true;
                    temp.add(other);
                    boolean mod = false;
                    if (columns != null) {
                        if (other.getColumns() != null) {
                            if (!columns.equals(other.getColumns())) {
                                mod = true;
                                details = details.concat("\n"+Indent.INDEX.getValue()+"columns=\""+columns+"\" CHANGED TO \""+other.getColumns()+"\"\n");
                            }
                        } else {
                            mod = true;
                            details = details.concat("\n"+Indent.INDEX.getValue()+"columns=\""+columns+"\" REMOVED\n");
                        }
                    } else {
                        if (other.getColumns() != null) {
                            mod = true;
                            details = details.concat("\n"+Indent.INDEX.getValue()+"columns=\""+other.getColumns()+"\" ADDED\n");
                        }
                    }

                    if (initial != null) {
                        if (other.getInitial() != null) {
                            if (!initial.equals(other.getInitial())) {
                                mod = true;
                                details = details.concat("\n"+Indent.INDEX.getValue()+"initial=\""+initial+"\" CHANGED TO \""+other.getInitial()+"\"\n");
                            }
                        } else {
                            mod = true;
                            details = details.concat("\n"+Indent.INDEX.getValue()+"initial=\""+initial+"\" REMOVED\n");
                        }
                    } else {
                        if (other.getInitial() != null) {
                            mod = true;
                            details = details.concat("\n"+Indent.INDEX.getValue()+"initial=\""+other.getInitial()+"\" ADDED\n");
                        }
                    }
                    String s = compareVariants(Indent.INDEX, index.getVariant(), other.getVariant());
                    if (s.length() > 0) {
                        mod = true;
                        details = details.concat("\n"+s);
                    }
                    if (mod) {
                        buff.append(idStr+" MODIFED "+details);
                    }
                }
            }

            if (!matches) {
                buff.append(idStr+" REMOVED\n");
            }
        }

        for (Index other: comp) {
            boolean matches = false;
            for (Index key: temp) {
                if (other.getName().equalsIgnoreCase(key.getName())) {
                    matches = true;
                }
            }

            if (!matches) {
                buff.append(Indent.TABLE.getValue()+"Index: name=\""+other.getName()+"\" ADDED\n");
                if (other.getColumns() != null) {
                    buff.append(Indent.INDEX.getValue()+"columns=\""+other.getColumns()+"\"\n");
                }
                for (Variant variant: other.getVariant()) {
                    buff.append(Indent.INDEX.getValue()+"Variant: dialect=\""+variant.getDialect()+"\"\n");
                    buff.append(Indent.INDEX.getValue()+Indent.VARIANT.getValue()+"columns=\""+variant.getColumns()+"\"\n");
                }
            }
        }
        return buff.toString();
    }

    private String compareVariants(Indent indent, List<Variant> base, List<Variant> comp) {
        ArrayList<Variant> temp = new ArrayList<Variant>();
        StringBuffer buff = new StringBuffer("");
        for (Variant variant: base) {
            boolean matches = false;
            String dialect = variant.getDialect();
            String columns = variant.getColumns();
            String idStr = indent.getValue()+"Variant: dialect=\""+dialect+"\"";
            String details = "";
            for (Variant other: comp) {
                if (dialect.equalsIgnoreCase(other.getDialect())) {
                    matches = true;
                    temp.add(other);
                    if (!columns.equals(other.getColumns())) {
                        buff.append(idStr+"\n"+indent.getValue()+Indent.VARIANT.getValue()+"columns=\""+columns+"\" CHANGED TO \""+other.getColumns()+"\"\n");
                    }
                }
            }

            if (!matches) {
                buff.append(idStr+" columns=\""+columns+"\" REMOVED\n");
            }
        }

        for (Variant other: comp) {
            boolean matches = false;
            for (Variant variant: temp) {
                if (other.getDialect().equalsIgnoreCase(variant.getDialect())) {
                    matches = true;
                }
            }

            if (!matches) {
                buff.append(indent.getValue()+"Variant: dialect: \""+other.getDialect()+"\" ADDED\n");
                buff.append(indent.getValue()+Indent.VARIANT.getValue()+"columns=\""+other.getColumns()+"\"\n");
            }
        }
        return buff.toString();
    }

    private com.sri.biospice.warehouse.schema.tools.types.Enumeration getEnumeration(List list) {
        com.sri.biospice.warehouse.schema.tools.types.Enumeration tmp = null;
        for (Object obj: list) {
            if (obj instanceof com.sri.biospice.warehouse.schema.tools.types.Enumeration) {
                tmp = (com.sri.biospice.warehouse.schema.tools.types.Enumeration)obj;
            } else if (obj instanceof JAXBElement) {
                JAXBElement jobj = (JAXBElement) obj;
                if (jobj.getDeclaredType().equals(com.sri.biospice.warehouse.schema.tools.types.Enumeration.class)) {
                    tmp = (com.sri.biospice.warehouse.schema.tools.types.Enumeration)jobj.getValue();
                }
            }
        }
        if (tmp != null) {
            for (com.sri.biospice.warehouse.schema.tools.types.Enum e:tmp.getRestriction()) {
            }
        }
        return tmp;
    }

    private String compareEnumerations(Indent indent, com.sri.biospice.warehouse.schema.tools.types.Enumeration base, com.sri.biospice.warehouse.schema.tools.types.Enumeration comp) {
        StringBuffer buff = new StringBuffer("");
        ArrayList<com.sri.biospice.warehouse.schema.tools.types.Enum> tmp = new ArrayList<com.sri.biospice.warehouse.schema.tools.types.Enum>();
        if (base != null) {
            if (comp != null) {
                for (com.sri.biospice.warehouse.schema.tools.types.Enum enm: base.getRestriction()) {
                    boolean matches = false;
                    String value = enm.getValue();
                    if (value == null) {
                        continue;
                    }

                    String idStr = indent.getValue()+Indent.ENUMERATION.getValue()+"Restriction: value: \""+value+"\"";
                    String details = "";
                    for (com.sri.biospice.warehouse.schema.tools.types.Enum other: comp.getRestriction()) {
                        if (value.equalsIgnoreCase(other.getValue())) {
                            matches = true;
                            tmp.add(enm);
                            String desc = enm.getDescription();
                            boolean mod = false;
                            if (desc != null) {
                                if (other.getDescription() != null) {
                                    if (!desc.equals(other.getDescription())) {
                                        mod = true;
                                        details = details.concat(indent.getValue()+Indent.ENUMERATION.getValue()+Indent.RESTRICTION.getValue()+"description=\""+desc+"\" CHANGED TO \""+other.getDescription()+"\"\n");
                                    }
                                } else {
                                    mod = true;
                                    details = details.concat(indent.getValue()+Indent.ENUMERATION.getValue()+Indent.RESTRICTION.getValue()+"description=\""+desc+"\" REMOVED\n");
                                }
                            } else {
                                if (other.getDescription() != null) {
                                    mod = true;
                                    details = details.concat(indent.getValue()+Indent.ENUMERATION.getValue()+Indent.RESTRICTION.getValue()+"description=\""+other.getDescription()+"\" ADDED\n");
                                }
                            }
                            if (mod) {
                                buff.append(idStr+" MODIFIED\n"+details);
                            }
                        }
                    }
                    if (!matches) {
                        buff.append(idStr+" REMOVED \n");
                    }
                }

                for (com.sri.biospice.warehouse.schema.tools.types.Enum other: comp.getRestriction()) {
                    boolean matches = false;
                    String value = other.getValue();
                    for (com.sri.biospice.warehouse.schema.tools.types.Enum enm: tmp) {
                        if (value.equalsIgnoreCase(enm.getValue())) {
                            matches = true;
                        }
                    }

                    if (!matches) {
                        buff.append(indent.getValue()+Indent.ENUMERATION.getValue()+"Restriction: value=\""+value+"\" description=\""+other.getDescription()+"\" ADDED\n");
                    }
                }
            } else {
                for (com.sri.biospice.warehouse.schema.tools.types.Enum enm: comp.getRestriction()) {
                    buff.append(indent.getValue()+Indent.ENUMERATION.getValue()+"Restriction: value=\""+enm.getValue()+"\" REMOVED\n");
                }
            }
        } else {
            if (comp != null) {
                for (com.sri.biospice.warehouse.schema.tools.types.Enum other: comp.getRestriction()) {
                    buff.append(indent.getValue()+Indent.ENUMERATION.getValue()+"Restriction: value=\""+other.getValue()+"\" description=\""+other.getDescription()+"\" ADDED\n");
                }
            }
        }
        return buff.toString();
    }

    private String compareStringLists(Indent indent, String elementName, List<String> base, List<String> comp) {
        StringBuffer buff = new StringBuffer("");
        if (base != null) {
            if (comp != null) {
                for (String val: base) {
                    boolean matches = false;
                    for (String other: comp) {
                        if (val.equalsIgnoreCase(other)) {
                            matches = true;
                        }
                    }
                    if (!matches) {
                        buff.append(indent.getValue()+elementName+": \""+val+"\" REMOVED\n");
                    }
                }
                for (String val: comp) {
                    boolean matches = false;
                    for (String other: base) {
                        if (val.equalsIgnoreCase(other)) {
                            matches = true;
                        }
                    }

                    if (!matches) {
                        buff.append(indent.getValue()+elementName+": \""+val+"\" ADDED\n");
                    }
                }
            } else {
                for (String val: base) {
                    buff.append(indent.getValue()+elementName+": \""+val+"\" REMOVED\n");
                }
            }
        } else {
            if (comp != null) {
                for (String val: comp) {
                    buff.append(indent.getValue()+elementName+": \""+val+"\" ADDED\n");
                }
            }
        }
        return buff.toString();
    }

    private static void usage() {
        StringBuffer buff = new StringBuffer("\n\nUsage: XmlSchemaDiffTool file1 file2 [log file]\n");
        buff.append("       where file1 is the base file for comparison\n");
        buff.append("       where file2 is the file to compare with\n");
        buff.append("       where log file is the name of the file to write the results to - OPTIONAL\n\n");
        logger.log(Level.SEVERE, buff.toString());
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String f1 = "data/oldSchema.xml";
        String f2 = "data/newSchema.xml";
        String lf = null;
        if (args.length ==2) {
            f1 = args[0];
            f2 = args[1];
        } else if (args.length == 3) {
            f1 = args[0];
            f2 = args[1];
            lf = args[2];
        } else {
            usage();
            System.exit(1);
        }

        try {
            SchemaDiffTool xd = new SchemaDiffTool(f1, f2, lf);
        } catch (InstantiationException ex) {
            ex.printStackTrace();
        }
    }

}
