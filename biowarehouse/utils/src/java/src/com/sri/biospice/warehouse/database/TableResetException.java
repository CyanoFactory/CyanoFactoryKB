/*
 * TableResetException.java
 *
 * Created on October 26, 2006, 10:46 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package com.sri.biospice.warehouse.database;

/**
 *
 * @author gpetry
 */
public class TableResetException extends java.lang.Exception {
    
    /**
     * Creates a new instance of <code>TableResetException</code> without detail message.
     */
    public TableResetException() {
    }
    
    
    /**
     * Constructs an instance of <code>TableResetException</code> with the specified detail message.
     * @param msg the detail message.
     */
    public TableResetException(String msg) {
        super(msg);
    }
}
