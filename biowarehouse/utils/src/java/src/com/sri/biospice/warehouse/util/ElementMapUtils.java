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

import com.sri.bw.LoaderStatistics;
import com.sri.bw.element.*;
import org.apache.log4j.Logger;
import org.apache.xmlbeans.XmlCursor;
import javax.xml.namespace.QName;
import java.util.Vector;

/**
 * @author Valerie Wagner
 *         Date: Aug 27, 2006
 */
public class ElementMapUtils {

    private static final Logger log = Logger.getLogger(ElementMapUtils.class);

    public static final int NOT_FOUND = -100;

    public static Vector<OneToOne> findOneToOnes(ModelDocument.Model model, String targetTableName) {
        Vector<OneToOne> list = new Vector<OneToOne>();
        for (ElementDocument.Element element : model.getElementArray()) {
            for (OneToOne oneToOne : element.getOneToOneArray()) {
                if (targetTableName.equals(oneToOne.getTarget())) {
                    list.add(oneToOne);
                }
            }
        }

        return list;
    }
    
    public static Vector<OneToMany> findOneToMany(ModelDocument.Model model, String name) {
        Vector<OneToMany> list = new Vector<OneToMany>();
        for (ElementDocument.Element elem : model.getElementArray()) {
            for (OneToMany one2many : elem.getOneToManyArray()) {
                if (name.equals(one2many.getName())) {
                    list.add(one2many);
                }
            }
        }
        return list;
    }

    public static ElementDocument.Element findElement(ModelDocument.Model model, String elementName) {
        for (ElementDocument.Element elem : model.getElementArray()) {
            if (elementName.equals(elem.getName())) {
                return elem;
            }
        }
        return null;
    }

    public static OneToMany findOneToMany(ElementDocument.Element element, String name) {
        for (OneToMany one2many : element.getOneToManyArray()) {
            if (name.equals(one2many.getName())) {
                return one2many;
            }
        }
        return null;
    }

    public static void removeAttribute(ElementDocument.Element element, String name) {
        int index = findIndexOfAttribute(element, name);
        if (index > -1) {
            element.removeAttribute(index);
        } else {
            log.warn("Attribute " + name + " not found");
            LoaderStatistics.warningOccurred();
        }
    }

    private static int findIndexOfAttribute(ElementDocument.Element element, String name) {
        for (int i = 0; i < element.getAttributeArray().length; i++) {
            if (name.equals(element.getAttributeArray(i).getName())) {
                return i;
            }
        }
        return -1;
    }

    public static Attribute findAttribute(ElementDocument.Element element, String name) {
        for (Attribute att : element.getAttributeArray()) {
            if (name.equals(att.getName())) {
                return att;
            }
        }
        return null;
    }

    public static void removeOneToMany(ElementDocument.Element element, String name) {
        int index = findIndexOfOneToMany(element, name);
        if (index > -1) {
            element.removeOneToMany(index);
        } else {
            log.warn("OneToMany " + name + " not found");
            LoaderStatistics.warningOccurred();
        }
    }

    public static int findIndexOfOneToMany(ElementDocument.Element element, String name) {
        for (int i = 0; i < element.getOneToManyArray().length; i++) {
            if (name.equals(element.getOneToManyArray(i).getName())) {
                return i;
            }
        }
        return -1;
    }


    public static int findAttributeIndex(ElementDocument.Element element, String attributeName) {
        for (int i = 0; i < element.getAttributeArray().length; i++) {
            if (attributeName.equals(element.getAttributeArray(i).getName())) {
                return i;
            }
        }
        return NOT_FOUND;
    }

    public static boolean hasAttribute(ElementDocument.Element element, String attributeName) {
        return findAttributeIndex(element, attributeName) != NOT_FOUND;
    }

    public static int findOneToOneIndex(ElementDocument.Element element, String assocName) {
        for (int i = 0; i < element.getOneToOneArray().length; i++) {
            if (assocName.equals(element.getOneToOneArray(i).getName())) {
                return i;
            }
        }
        return NOT_FOUND;
    }

    public static OneToOne findOneToOne(ElementDocument.Element element, String assocName) {
        int index = findOneToOneIndex(element, assocName);
        if (index == NOT_FOUND) {
            return null;
        }
        return element.getOneToOneArray(index);
    }

    public static void removeElement(ModelDocument.Model model, String elementName) {
        int index = findElementIndex(model, elementName);
        if (index != NOT_FOUND) {
            model.removeElement(index);
        } else {
            log.warn("Element " + elementName + " not found");
            LoaderStatistics.warningOccurred();
        }
    }

    private static int findElementIndex(ModelDocument.Model model, String elementName) {
        for (int i = 0; i < model.getElementArray().length; i++) {
            if (elementName.equals(model.getElementArray(i).getName())) {
                return i;
            }
        }
        return NOT_FOUND;
    }

    public static void removeManyToOne(ElementDocument.Element element, String many2OneName) {
        int index = findManyToOneIndex(element, many2OneName);
        if (index != NOT_FOUND) {
            element.removeManyToOne(index);
        } else {
            log.warn("ManyToOne " + many2OneName + " not found");
            LoaderStatistics.warningOccurred();
        }
    }

    private static int findManyToOneIndex(ElementDocument.Element element, String many2OneName) {
        for (int i = 0; i < element.getManyToOneArray().length; i++) {
            if (many2OneName.equals(element.getManyToOneArray(i).getName())) {
                return i;
            }
        }
        return NOT_FOUND;
    }

    public static void removeAllWithAttribute(ModelDocument.Model model, String attributeName, String attributeValue) {

        QName targetAttribute = new QName(attributeName);
        XmlCursor cursor = model.newCursor();
        cursor.push();
        while (cursor.hasNextToken()) {
            cursor.toNextToken();
            String targetText = cursor.getAttributeText(targetAttribute);
            if (targetText != null) {
                if (targetText.equals(attributeValue)) {
                    cursor.removeXml();
                }
            }
        }
        cursor.pop();
    }

    public static void capitalizeAttributeColumns(ElementDocument.Element element) {
        for (Attribute att : element.getAttributeArray()) {
            att.setColumn(capFirst(att.getColumn()));
        }
    }


    public static String capFirst(String name) {
        return name.substring(0, 1).toUpperCase() + name.substring(1, name.length());
    }


    public static void removeOneToMany(ElementDocument.Element parent, OneToMany one2many) {
        removeOneToMany(parent, one2many.getName());
    }

    public static void removeOneToOne(ElementDocument.Element element, OneToOne one2one) {
        int index = findOneToOneIndex(element, one2one.getName());
        if (index != NOT_FOUND) {
            element.removeOneToOne(index);
        } else {
            log.warn("OneToOne " + one2one.getName() + " not found");
            LoaderStatistics.warningOccurred();
        }
    }

    public static void changeAttributeToSkip(ElementDocument.Element element, String attName) {
        Attribute att = findAttribute(element, attName);
        att.setSkip("true");
        att.setColumn("");
    }

    public static void changeElementToSkip(ModelDocument.Model model, String elementName) {
        ElementDocument.Element element = findElement(model, elementName);
        element.setSkip2("true");
        element.setTable("");
        element.newCursor().removeXmlContents();
    }

    public static void changeAllAttributesToSkip(ModelDocument.Model model, String attributeName, String attributeValue) {
        QName targetAttribute = new QName(attributeName);
        XmlCursor cursor = model.newCursor();
        cursor.push();
        while (cursor.hasNextToken()) {
            cursor.toNextToken();
            String targetText = cursor.getAttributeText(targetAttribute);
            if (targetText != null) {
                if (targetText.equals(attributeValue)) {
                    cursor.removeXmlContents();
                    // todo: figure out how to remove other attributes or set them to blank
                    cursor.setAttributeText(new QName("skip"), "true");
                }
            }
        }
        cursor.pop();
    }

    public static void changeOneToManyToSkip(ElementDocument.Element element, String name) {
        OneToMany one2many = findOneToMany(element, name);
        one2many.setSkip("true");
        one2many.setTarget("");
        one2many.setColumn("");
    }

    public static void changeManyToOneToSkip(ElementDocument.Element element, String many2OneName) {
        ManyToOne index = findManyToOne(element, many2OneName);
        index.setSkip("true");
        index.setTarget("");
        index.setColumn("");
    }

    private static ManyToOne findManyToOne(ElementDocument.Element element, String many2OneName) {
        int index = findManyToOneIndex(element, many2OneName);
        if (index != NOT_FOUND) {
            return element.getManyToOneArray(index);
        }
        return null;
    }

    public static OneToMany changeOneToOneIntoOneToMany(ElementDocument.Element element, String assocName) {
        OneToOne one2one = findOneToOne(element, assocName);
        OneToMany one2many = element.addNewOneToMany();
        if (one2one.isSetSkip()) {
            one2many.setSkip(one2one.getSkip());
        }
        one2many.setColumn(one2one.getColumn());
        one2many.setTarget(one2one.getTarget());
        one2many.setName(one2one.getName());
        one2many.setDefaultArray(one2one.getDefaultArray());
        removeOneToOne(element, one2one);
        return one2many;
    }

    public static OneToMany changeOneToOneIntoOneToMany(ElementDocument.Element element, OneToOne one2one) {
        return changeOneToOneIntoOneToMany(element, one2one.getName());
    }

    public static void changeTargetColumn(ModelDocument.Model model, String targetName, String newColumn) {

        QName targetAttribute = new QName("target");
        QName columnAttribute = new QName("column");
        XmlCursor cursor = model.newCursor();
        cursor.push();
        while (cursor.hasNextToken()) {
            cursor.toNextToken();
            String targetText = cursor.getAttributeText(targetAttribute);
            if (targetText != null) {
                if (targetText.equals(targetName)) {
                    cursor.setAttributeText(columnAttribute, newColumn);
                }
            }
        }
        cursor.pop();
    }
    public static void changeTargets(ModelDocument.Model model, String targetName, String newTarget) {

        QName targetAttribute = new QName("target");
        XmlCursor cursor = model.newCursor();
        cursor.push();
        while (cursor.hasNextToken()) {
            cursor.toNextToken();
            String targetText = cursor.getAttributeText(targetAttribute);
            if (targetText != null) {
                if (targetText.equals(targetName)) {
                    cursor.setAttributeText(targetAttribute, newTarget);
                }
            }
        }
        cursor.pop();
    }

    public static void addElementToModel(ModelDocument.Model model, ElementDocument.Element addElement) {
        ElementDocument.Element[] elements = model.getElementArray();
        Vector<ElementDocument.Element> elementsVector = new Vector<ElementDocument.Element>(elements.length + 1);
        for (ElementDocument.Element elem : elements) {
            elementsVector.add(elem);
        }
        elementsVector.add(addElement);
        elements = elementsVector.toArray(elements);
        model.setElementArray(elements);

    }
}
