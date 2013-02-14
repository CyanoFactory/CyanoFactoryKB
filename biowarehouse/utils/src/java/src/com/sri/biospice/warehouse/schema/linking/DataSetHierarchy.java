package com.sri.biospice.warehouse.schema.linking;

import com.sri.biospice.warehouse.database.table.TableMetaData;
import com.sri.biospice.warehouse.schema.DataSet;

/**
 * @author Valerie Wagner
 *         Date: Jul 16, 2008
 */
public class DataSetHierarchy extends LinkingTable
{
    private static TableMetaData tableMetaData;


    /**
     * Used for a non-hierarchical dataset (one that 'contains' itself).
     * @param dataset
     */
    public DataSetHierarchy( DataSet dataset ) {
        this(dataset, dataset);
    }

    public DataSetHierarchy( DataSet dataset1, DataSet dataset2 )
    {
        super( "DataSetHierarchy", "SuperWID", "SubWID" );
        this.setWID1( dataset1.getWID() );
        this.setWID2( dataset2.getWID() );
    }


    protected void setTableMetaData( TableMetaData tableMetaData )
    {
        this.tableMetaData = tableMetaData;
    }


    public TableMetaData getTableMetaData()
    {
        return tableMetaData;
    }
}
