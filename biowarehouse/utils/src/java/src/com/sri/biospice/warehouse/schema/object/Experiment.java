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

import com.sri.biospice.warehouse.database.table.*;
import java.util.Date;
import java.sql.Timestamp;

/**
 * Created by IntelliJ IDEA.
 * User: kejariwa
 * Date: Oct 31, 2006
 * Time: 10:16:33 AM
 * To change this template use File | Settings | File Templates.
 */
public class Experiment extends ObjectTableBase {
    private StringColumn typeCol;
    private LongColumn contactWIDCol;
    private LongColumn archiveWIDCol;
    private TimestampColumn startDateCol;
    private TimestampColumn endDateCol;
    private ClobColumn descriptionCol;
    private LongColumn groupWIDCol;
    private StringColumn groupTypeCol;
    private IntegerColumn groupSizeCol;
    private IntegerColumn groupIndexCol;
    private IntegerColumn timePointCol;
    private StringColumn timeUnitCol;
    private LongColumn bioSourceWIDCol;

    private static TableMetaData tableMetaData = null;

    public Experiment(long datasetWid) {
        super(datasetWid, "Experiment");
        init();
    }

    /**
     * For an entry loaded from the Gene table in the database
     *
     * @param datasetWID
     * @param thisWID    The WID of the entry in the Gene table
     */
    public Experiment(long datasetWID, long thisWID) {
        super(datasetWID, "Experiment", thisWID);
        init();
    }

    private void init() {
        typeCol = new StringColumn("Type", this);
        contactWIDCol = new LongColumn("ContactWID",this);
        archiveWIDCol = new LongColumn("ArchiveWID",this);
        startDateCol = new TimestampColumn("StartDate",this);
        endDateCol = new  TimestampColumn("EndDate",this);
        descriptionCol = new ClobColumn("Description",this);
        groupWIDCol = new LongColumn("GroupWID",this);
        groupTypeCol = new StringColumn("GroupType",this);
        groupSizeCol = new IntegerColumn("GroupSize",this);
        groupIndexCol = new IntegerColumn("GroupIndex",this);
        timePointCol = new IntegerColumn("TimePoint",this);
        timeUnitCol = new StringColumn("TimeUnit",this);
        bioSourceWIDCol = new LongColumn("BioSourceWID",this);
    }


    public Long getBioSourceWID() {
        return bioSourceWIDCol.getValue();
    }

    public void setBioSourceWID(Long bioSourceWID) {
        this.bioSourceWIDCol.setValue(bioSourceWID);
    }

    public Long getContactWID() {
        return contactWIDCol.getValue();
    }

    public void setContactWID(Long contactWID) {
        this.contactWIDCol.setValue(contactWID);
    }

    public Long getArchiveWID() {
        return archiveWIDCol.getValue();
    }

    public void setArchiveWID(Long archiveWID) {
        this.archiveWIDCol.setValue(archiveWID);
    }

    public Date getStartDate() {
        return startDateCol.getValue();
    }

    public void setStartDate(Date startDate) {
        this.startDateCol.setValue(startDate);
    }

    public Date getEndDate() {
        return endDateCol.getValue();
    }

    public void setEndDate(Date endDateCol) {
        this.endDateCol.setValue(endDateCol);
    }

    public String getDescription() {
        return descriptionCol.getValue();
    }

    public void setDescription(String description) {
        this.descriptionCol.setValue(description);
    }

    public Long getGroupWID() {
        return groupWIDCol.getValue();
    }

    public void setGroupWID(Long groupWID) {
        this.groupWIDCol.setValue(groupWID);
    }

    public String getGroupType() {
        return groupTypeCol.getValue();
    }

    public void setGroupType(String groupType) {
        this.groupTypeCol.setValue(groupType);
    }

    public Integer getGroupSize() {
        return groupSizeCol.getValue();
    }

    public void setGroupSize(Integer groupSize) {
        this.groupSizeCol.setValue(groupSize);
    }

    public Integer getGroupIndex() {
        return groupIndexCol.getValue();
    }

    public void setGroupIndex(Integer groupIndex) {
        this.groupIndexCol.setValue(groupIndex);
    }

    public Integer getTimePoint() {
        return timePointCol.getValue();
    }

    public void setTimePoint(Integer timePoint) {
        this.timePointCol.setValue(timePoint);
    }

    public String getTimeUnit() {
        return timeUnitCol.getValue();
    }

    public void setTimeUnit(String timeUnit) {
        this.timeUnitCol.setValue(timeUnit);
    }

    public String getType() {
        return typeCol.getValue();
    }

    public void setType(String type) {
        this.typeCol.setValue(type);
    }

    public TableMetaData getTableMetaData() {
        return tableMetaData;
    }


    protected void setTableMetaData(TableMetaData tableMetaData) {
        this.tableMetaData = tableMetaData;
    }


}
