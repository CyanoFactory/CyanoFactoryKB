<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    <!--BGCOLOR="#CCCCFF"-->
    <xsl:template match="/schema">
        <xsl:result-document href="index.html" method="html">
            <html>
                <head>
                    <link href="style.css" rel="stylesheet" type="text/css"/>
                    <title>SchemaDoc: Table Index</title>
                </head>
                <body>
                    <h1>BioWarehouse Schema</h1>
                    <a href="#Index"><h4>Table Index</h4></a>
                       The following diagram represents the relationships between the core tables of the BioWarehouse schema
                    (excluding the MAGE schema extension).
                    Click on a table name to view its documentation, or see the <a href="#Index">Table Index</a> below for the complete
                    list of BioWarehouse tables.  Object tables are colored yellow, and linking (associative) tables are colored blue.
                    <img src="core-schema-diagram.png" alt="schema diagram" usemap="#coreschema"/>
                    <h1><a name="Index"/>Table Index</h1>
                    <br/>
                    <table border="1">
                        <tr>
                            <th>Table</th>
                            <th>Comment</th>
                        </tr>
                        <xsl:for-each select="table">
                            <xsl:sort select="@name"/>
                            <tr>
                                <td>
                                    <a href="{@name}.html">
                                        <xsl:value-of select="@name"/>
                                    </a>
                                </td>
                                <td>
                                    <xsl:for-each select="comment">
                                        <xsl:value-of select="."/>
                                        <br/>
                                    </xsl:for-each>
                                </td>
                            </tr>
                        </xsl:for-each>
                    </table>

                    <!-- Summary information -->
                    <h2>Summary:</h2>
                    <xsl:value-of select="count(table)"/>
                    <xsl:text> Total tables</xsl:text>
                    <br/>
                    <xsl:value-of select="count(table[@type='object'])"/>
                    <xsl:text> Object tables</xsl:text>
                    <br/>
                    <xsl:value-of select="count(table[@type='associative'])"/>
                    <xsl:text> Linking (associative) tables</xsl:text>
                    <br/>
                    <xsl:value-of
                            select="count(table) - count(table[@type='object']) - count(table[@type='associative'])"/>
                    <xsl:text> Other tables</xsl:text>
                    <br/>

                </body>
            </html>

            <!--Create page for each table-->
            <xsl:for-each select="table">
                <xsl:sort select="@type" order="descending"/>
                <xsl:sort select="@name"/>
                <xsl:call-template name="one.table"/>
            </xsl:for-each>

        </xsl:result-document>
    </xsl:template>

    <!--Create page for one table-->
    <xsl:template name="one.table" match="table">
        <xsl:result-document href="{@name}.html" method="html">
            <html>
                <head>
                    <link href="style.css" rel="stylesheet" type="text/css"/>
                    <title>SchemaDoc:
                        <xsl:value-of select="@name"/>
                    </title>
                </head>
                <body>
                    <a href="index.html">Table Index</a>

                    <h1>
                        <a name="{@name}"/>
                        <xsl:value-of select="@name"/>
                    </h1>
                    <hr/>
                    <xsl:for-each select="comment">
                        <xsl:value-of select="."/>
                        <br/>
                    </xsl:for-each>
                    <h2>Columns</h2>
                    <table border="1">
                        <tr>
                            <th>Column</th>
                            <th>MySQL Type</th>
                            <th>Oracle Type</th>
                            <th>Nullable</th>
                            <th>Description</th>
                        </tr>
                        <xsl:for-each select="column">
                            <tr>
                                <!--Column-->
                                <td>
                                    <xsl:choose>
                                        <xsl:when test="foreignKey">
                                            <a href="{foreignKey/@toTable}.html">
                                                <xsl:value-of select="@name"/>
                                            </a>
                                        </xsl:when>
                                        <xsl:otherwise>
                                            <xsl:value-of select="@name"/>
                                        </xsl:otherwise>
                                    </xsl:choose>
                                </td>
                                <!--MySQL Type-->
                                <td>
                                    <xsl:value-of select="/schema/datatypes/dtype[@notation=current()/@type]/@mysql"/>
                                </td>
                                <!--Oracle Type-->
                                <td>
                                    <xsl:value-of select="/schema/datatypes/dtype[@notation=current()/@type]/@oracle"/>
                                    <xsl:call-template name="type.size"/>
                                </td>
                                <!--Nullable-->
                                <td>
                                    <xsl:choose>
                                        <xsl:when test="@required='true'">No</xsl:when>
                                        <xsl:otherwise>
                                            <xsl:text>Yes</xsl:text>
                                        </xsl:otherwise>
                                    </xsl:choose>
                                </td>
                                <!--Description-->
                                <td>
                                    <xsl:for-each select="comment">
                                        <xsl:value-of select="."/>
                                        <br/>
                                    </xsl:for-each>
                                    <xsl:for-each select="enumeration">
                                        <br/>
                                        <div style="margin-left: 40px;">
                                        <b><i>Enumerated Values:</i></b><br/>
                                        <xsl:for-each select="restriction">
                                               <b><xsl:value-of select="@value"/></b> -
                                               <xsl:value-of select="@description"/> <br/>
                                        </xsl:for-each>
                                        </div>
                                    </xsl:for-each>
                                </td>
                            </tr>
                        </xsl:for-each>
                    </table>
                    <xsl:if test="column/sdKey">
                        <h2>Denormalized References</h2>
                        <table border="1">
                            <tr>
                                <th>Column</th>
                                <th>References</th>
                            </tr>
                            <xsl:for-each select="column/sdKey">
                                <tr>
                                    <td>
                                        <xsl:value-of select="../@name"/>
                                    </td>
                                    <td>
                                        <a href="{@toTable}.html">
                                            <xsl:value-of select="@toTable"/>
                                        </a>
                                        :
                                        <xsl:value-of select="@toColumn"/>
                                    </td>
                                </tr>
                            </xsl:for-each>
                        </table>
                    </xsl:if>
                    <h2>Referenced By</h2>
                    <xsl:variable name="fkeys" select="/schema/table/column/foreignKey[@toTable=current()/@name]"/>
                    <xsl:variable name="sdkeys" select="/schema/table/column/sdKey[@toTable=current()/@name]"/>
                    <xsl:choose>
                        <xsl:when test="count($fkeys) = 0 and count($sdkeys) = 0">
                            <xsl:text>None.</xsl:text>
                        </xsl:when>
                        <xsl:otherwise>

                            <table border="1">
                                <tr>
                                    <th>Table</th>
                                    <th>Column</th>
                                </tr>
                                <!--for each foreign key that maps to this table-->
                                <xsl:for-each select="$fkeys">
                                    <tr>
                                        <td>
                                            <!--table name-->
                                            <a href="{../../@name}.html">
                                                <xsl:value-of select="../../@name"/>
                                            </a>
                                        </td>
                                        <td>
                                            <!--column name-->
                                            <xsl:value-of select="../@name"/>
                                        </td>
                                    </tr>
                                </xsl:for-each>
                                <xsl:for-each select="$sdkeys">
                                    <tr>
                                        <td>
                                            <!--table name-->
                                            <a href="{../../@name}.html">
                                                <xsl:value-of select="../../@name"/>
                                            </a>
                                        </td>
                                        <td>
                                            <!--column name-->
                                            <xsl:value-of select="../@name"/>
                                        </td>
                                    </tr>
                                </xsl:for-each>
                            </table>

                        </xsl:otherwise>
                    </xsl:choose>

                    <h2>Other Constraints</h2>
                    <xsl:choose>
                        <xsl:when test="constraint">
                            <xsl:for-each select="constraint">
                                <xsl:value-of select="."/>
                                <br/>
                            </xsl:for-each>

                        </xsl:when>
                        <xsl:otherwise>None.</xsl:otherwise>
                    </xsl:choose>

                    <h2>Indexes</h2>
                    <xsl:choose>
                        <xsl:when test="index">
                            <table border="1">
                                <tr>
                                    <th>Name</th>
                                    <th>Columns</th>
                                </tr>
                                <xsl:for-each select="index">
                                    <tr>
                                        <td>
                                            <xsl:value-of select="@name"/>
                                        </td>
                                        <td>
                                            <xsl:value-of select="@columns"/>
                                        </td>
                                    </tr>
                                </xsl:for-each>
                            </table>
                        </xsl:when>
                        <xsl:otherwise>None.</xsl:otherwise>
                    </xsl:choose>

                    <br/>
                    <br/>
                </body>
            </html>
        </xsl:result-document>
    </xsl:template>

    <xsl:template name="type.size" match="column">
        <xsl:choose>
            <xsl:when test="@length">
                <xsl:text>(</xsl:text>
                <xsl:value-of select="@length"/>
                <xsl:text>)</xsl:text>
            </xsl:when>
            <xsl:when test="@precision">
                <xsl:text>(</xsl:text>
                <xsl:value-of select="@precision"/>
                <xsl:text>,</xsl:text>
                <xsl:value-of select="@scale"/>
                <xsl:text>)</xsl:text>
            </xsl:when>
        </xsl:choose>
    </xsl:template>

</xsl:stylesheet>