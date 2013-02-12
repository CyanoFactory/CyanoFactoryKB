/*
    #=========================================================================
    # Copyright 2006 SRI International.  All rights reserved.
    #
    # The material contained in this file is confidential and proprietary to SRI
    # International and may not be reproduced, published, or disclosed to others
    # without authorization from SRI International.
    #
    # DISCLAIMER OF WARRANTIES
    #
    # SRI International MAKES NO REPRESENTATIONS OR WARRANTIES ABOUT THE
    # SUITABILITY OF THE SOFTWARE, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT
    # LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
    # PARTICULAR PURPOSE, OR NON-INFRINGEMENT. SRI International SHALL NOT BE
    # LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A RESULT OF USING, MODIFYING
    # OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES
    #=========================================================================
*/
package com.sri.biospice.warehouse.schema;

import com.sri.biospice.warehouse.database.table.LongColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;
import com.sri.biospice.warehouse.database.table.IntegerColumn;
import com.sri.biospice.warehouse.schema.object.Interaction;
import com.sri.biospice.warehouse.schema.object.NucleicAcid;
import com.sri.biospice.warehouse.schema.object.Protein;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * InteractionParticipant
 *
 * @author David Dunkley
 * @version 0.1 Nov 10, 2006
 */
public class InteractionParticipant extends TableBase {
    // Columns of the table
    private LongColumn interactionWIDCol;
    private LongColumn otherWIDCol;
    private IntegerColumn coefficientCol;

    // Table metadata
    private static TableMetaData tableMetaData;

    // Hashcode used for comparing DBIDs when placed in a hash-based data structure
    private int hashCode = 0;


    public InteractionParticipant(Interaction interaction, Protein protein) {
        super("InteractionParticipant");
        init();
        interactionWIDCol.setValue(interaction.getWID());
        otherWIDCol.setValue(protein.getWID());
    }

    public InteractionParticipant(Interaction interaction, NucleicAcid nucleicAcid) {
        super("InteractionParticipant");
        init();
        interactionWIDCol.setValue(interaction.getWID());
        otherWIDCol.setValue(nucleicAcid.getWID());
    }

    /**
     * Initializes columns
     */
    private void init() {
        interactionWIDCol = new LongColumn("InteractionWID", this);
        otherWIDCol = new LongColumn("OtherWID", this);
        coefficientCol = new IntegerColumn("Coefficient", this);
    }


    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


    protected ResultSet getRow() throws SQLException {
        return null;
    }

    public long getInteractionWID() {
        return interactionWIDCol.getValue();
    }

    public long getOtherWID() {
        return otherWIDCol.getValue();
    }

    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        InteractionParticipant that = (InteractionParticipant) o;

        return interactionWIDCol.equals(that.interactionWIDCol) && otherWIDCol.equals(that.otherWIDCol);
    }

    public int hashCode() {
        if (hashCode == 0) {
            hashCode = interactionWIDCol.hashCode();
            hashCode = 31 * hashCode + otherWIDCol.hashCode();
        }
        return hashCode;
    }
}
