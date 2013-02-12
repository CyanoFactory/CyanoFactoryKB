/* $Id: EnzymaticReaction.java,v 1.1 2006/07/07 15:03:37 valerie Exp $ */
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

import com.sri.biospice.warehouse.database.table.LongColumn;
import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;

/**
 * @author Valerie.Wagner@sri.com
 *         Date: Sep 21, 2004
 */
public class EnzymaticReaction extends ObjectTableBase {
    // Metadata about the EnzymaticReaction table
    private static TableMetaData tableMetaData;

    // Columns in the EnzymaticReaction table
    private LongColumn reactionWIDCol;
    private LongColumn proteinWIDCol;
    private LongColumn complexWIDCol;
    private StringColumn reactionDirectionCol;

    // Enumerated values for the ReactionDirection column
//    public static final String ENUM_REACTION_DIRECTION_reversible = "reversible";
//    public static final String ENUM_REACTION_DIRECTION_physiol_left_to_right = "physiol-left-to-right";
//    public static final String ENUM_REACTION_DIRECTION_physiol_right_to_left = "physiol-right-to-left";
//    public static final String ENUM_REACTION_DIRECTION_irreversible_left_to_right = "irreversible-left-to-right";
//    public static final String ENUM_REACTION_DIRECTION_irreversible_right_to_left = "irreversible-right-to-left";

    /**
     * Enumerated values allowed for column <i>EnzymaticReaction.ReactionDirection</i>.
     */
    public enum ReactionDirection {
        reversible("reversible"),
        physiol_left_to_right("physiol-left-to-right"),
        physiol_right_to_left("physiol-right-to-left"),
        irreversible_left_to_right("irreversible-left-to-right"),
        irreversible_right_to_left("irreversible-right-to-left");

        private String direction;

        ReactionDirection(String dir) {
            this.direction = dir;
        }

        public String toString() {
            return direction;
        }
    }

    /**
     * Creates a new entry in the EnzymaticReaction table
     *
     * @param datasetWID
     */
    public EnzymaticReaction(long datasetWID) {
        super(datasetWID, "EnzymaticReaction");
        init();
    }


    /**
     * Represents an already-existing entry in the EnzymaticReaction table
     *
     * @param datasetWID
     * @param enzymaticReactionWID
     */
    public EnzymaticReaction(long datasetWID, long enzymaticReactionWID) {
        super(datasetWID, "EnzymaticReaction", enzymaticReactionWID);
        init();
    }


    private void init() {
        reactionWIDCol = new LongColumn("ReactionWID", this);
        proteinWIDCol = new LongColumn("ProteinWID", this);
        complexWIDCol = new LongColumn("ComplexWID", this);
        reactionDirectionCol = new StringColumn("ReactionDirection", this);
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    public void setReactionWID(long reactionWID) {
        reactionWIDCol.setValue(reactionWID);
    }


    public void setProteinWID(long proteinWID) {
        proteinWIDCol.setValue(proteinWID);
    }


    public void setComplexWID(long complexWID) {
        complexWIDCol.setValue(complexWID);
    }


    public void setReactionDirection(ReactionDirection direction) {
        reactionDirectionCol.setValue(direction);
    }


    /**
     * For testing purposes
     *
     * @param columnName Name of the columne to get the allowed enumerated values for
     * @return A Vector of the allowed enumerated values
     */
    public Enum[] getAllowedValues(String columnName) {
        if (columnName.equalsIgnoreCase(reactionDirectionCol.getColumnName())) {
            return ReactionDirection.values();
        }
        return null;
    }

}
