<!-- schema-to-dot.xsl - translates schema into Graphviz graph .dot file -->
<!-- Translated into diagram by Graphviz 'neato' program, http://www.graphviz.org/ -->
<!-- Valerie Wagner, November 2006 -->
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    <xsl:output method="text" />
    <xsl:template match="/schema">

<!-- See list of colors at: http://www.graphviz.org/doc/info/colors.html -->

digraph coreschema {
    graph [overlap=false, start=239487, splines=true];
    edge [ arrowsize=0.5];
    node [shape=box, fontsize=10];
    {node [color=lightgoldenrod2, style=filled]
        <xsl:for-each select="table">
            <xsl:choose>
                <xsl:when test="@type='object'">
                    <xsl:value-of select="@name"/>
                    [URL="<xsl:value-of select="@name"/>.html" ]
                </xsl:when>
            </xsl:choose>
        </xsl:for-each>
    }

    {node [color=lightblue3, style=filled]
        <xsl:for-each select="table">
            <xsl:choose>
                <xsl:when test="@type='associative'">
                    <xsl:value-of select="@name"/>
                    [URL="<xsl:value-of select="@name"/>.html" ]
                </xsl:when>
            </xsl:choose>
        </xsl:for-each>
    }

        {node [shape=box]
            <xsl:for-each select="table">
                <xsl:choose>
                    <xsl:when test="@type!='associative' and @type!='object'">
                        <xsl:value-of select="@name"/>
                        [URL="<xsl:value-of select="@name"/>.html" ]
                    </xsl:when>
                </xsl:choose>
            </xsl:for-each>
        }



        <!-- List connections -->
        <xsl:for-each select="table">
            <xsl:variable name="theTable" select="."/>
            <xsl:for-each select="column/foreignKey">
                <xsl:choose>
                    <xsl:when test="@toTable='DataSet'"/>
                    <xsl:otherwise>
                          "<xsl:value-of select="$theTable/@name"/>"
                          <xsl:text> -> </xsl:text>
                          "<xsl:value-of select="@toTable"/>"
                        ;
                    </xsl:otherwise>
                </xsl:choose>
            </xsl:for-each>
        </xsl:for-each>
}
    </xsl:template>
</xsl:stylesheet>
