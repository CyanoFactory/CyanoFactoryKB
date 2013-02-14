/* $Id: Support.java,v 1.2 2007/11/30 00:09:29 nguo Exp $ */
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

import com.sri.biospice.warehouse.database.table.FloatColumn;
import com.sri.biospice.warehouse.database.table.LongColumn;
import com.sri.biospice.warehouse.database.table.StringColumn;
import com.sri.biospice.warehouse.database.table.TableMetaData;

/**
 * An entry in the Support table
 * @author Nan Guo
 */
public class Support extends ObjectTableBase {
  // Columns
  private LongColumn otherWIDCol;
  private StringColumn typeCol;
  private StringColumn evidenceTypeCol;
  private FloatColumn confidenceCol;
  
  private static TableMetaData tableMetaData = null;
  /**
   * Creates a new entry in the Support table
   * @param datasetWid
   */
  public Support( long datasetWid )
  {
    super( datasetWid, "Support" );
    init();
  }
  

  /**
   * Represents an already-existing entry in the Support table
   * @param datasetWid
   * @param supportWID
   */
  public Support( long datasetWid, long supportWID )
  {
    super( datasetWid, "Support", supportWID );
    init();
  }

  private void init()
  {
    otherWIDCol = new LongColumn("OtherWID", this);
    typeCol = new StringColumn( "Type", this );
    evidenceTypeCol = new StringColumn( "EvidenceType", this );
    confidenceCol = new FloatColumn("Confidence", this);
  }

  public long getOtherWID() {
      return otherWIDCol.getValue();
  }


  public void setOtherWID(long otherWID) {
      otherWIDCol.setValue(otherWID);
  }
  
  public String getType() {
      return typeCol.getValue();
  }


  public void setType(String type) {
      typeCol.setValue(type);
  }
  
  public String getEvidenceType() {
    return evidenceTypeCol.getValue();
  }
  
  public void setEvidenceType(String evidenceType) {
    evidenceTypeCol.setValue(evidenceType);
  }

  
  public float getConfidence() {
      return confidenceCol.getValue();
  }


  public void setConfidence(float confidence) {
      confidenceCol.setValue(confidence);
  }

  public TableMetaData getTableMetaData()
  {
    return tableMetaData;
  }


  protected void setTableMetaData( TableMetaData tableMetaData )
  {
    this.tableMetaData = tableMetaData;
  }
}
