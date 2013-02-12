/* $Id: XMLUtil.java,v 1.1 2006/07/07 15:03:38 valerie Exp $ */
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
package com.sri.biospice.warehouse.util;

import org.w3c.dom.CharacterData;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * This class has some utility functions
 * which makes searches in the DOM tree easier.
 * @author Priyanka Gupta
 * @author Valerie Wagner
 */
public class XMLUtil
{
    /**
     * Goes thru the tree rooted at element, and searches for a child
     * with tag name : tagName. There might be multiple nodes with same
     * tagName, and the index decided which one you are requesting.
     * Note: In most cases you would expect to see only one and would be
     * interested in the one for which the index is 0.
     * @param element Root of the DOM Tree which is to be searched.
     * @param tagName The tag Name of the Node that is to be searched.
     * @param index   The index of the Node in the list of matching Nodes.
     * @return Returns the data associated with that node.
     */
    public static String getCharacterDataByTagName( Element element, String tagName, int index )
    {
        if( (element == null) ||
            (tagName == null) )
        {
            return null;
        }

        String data = null;

        NodeList list = element.getElementsByTagName( tagName );
        if( (list != null) &&
            (list.getLength() > index) )
        {
            Node listNode = list.item( index );
            if( listNode != null && listNode.getFirstChild() != null )
            {
                data = ((CharacterData)listNode.getFirstChild()).getData();
            }
        }

        return data;
    }


    /**
     * Goes thru the tree rooted at element, and searches for a child
     * with tag name : tagName. There might be multiple nodes with same
     * tagName, and the index decided which one you are requesting.
     * Note: In most cases you would expect to see only one and would be
     * interested in the one for which the index is 0.
     * @param element       Root of the DOM Tree which is to be searched.
     * @param tagName       The tag Name of the Node that is to be searched.
     * @param attributeName The attribute name that we want the data from.
     * @param index         The index of the Node in the list of matching Nodes.
     * @return Returns the data associated with the attributeName of the node tagName
     */
    public static String getAttributeByTagName( Element element, String tagName, String attributeName, int index )
    {
        if( (element == null) ||
            (tagName == null) )
        {
            return null;
        }

        String data = null;

        NodeList list = element.getElementsByTagName( tagName );
        if( (list != null) &&
            (list.getLength() > index) )
        {
            Element listNode = (Element)list.item( index );
            if( listNode != null )
            {
                data = listNode.getAttribute( attributeName );
            }
        }
        return data;
    }


    /**
     * Searchs for the first matching subelement
     * @param element     The root of the tree to be searched
     * @param elementName The name of the subelement to be matched
     * @return The first Element matching elementName, or null if no match
     */
    public static Element searchForElement( Element element, String elementName, int index )
    {
        Element foundElement = null;
        NodeList list = element.getElementsByTagName( elementName );
        if( list != null && list.getLength() > 0 )
        {
            foundElement = (Element)list.item( index );
        }

        return foundElement;
    }


    /**
     * Searchs for the first matching subelement
     * @param element     The root of the tree to be searched
     * @param elementName The name of the subelement to be matched
     * @return The first Element matching elementName, or null if no match
     */
    public static Element searchForElement( Element element, String elementName )
    {
        return searchForElement( element, elementName, 0 );
    }
}
