<!-- schema-merge.xsl -->
<!-- This stylesheet is used to merge the 'core' and 'mage' components
of the biowarehouse schema into one document
(see the Ant 'merge-schemas' target)
author: Valerie Wagner, November 2006
-->
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    <xsl:output method="xml" indent="yes"/>
    <xsl:template match="/">
        <schema>
            <xsl:apply-templates select="schema"/>
            <xsl:apply-templates select="document('mage-schema.xml')/schema"/>
        </schema>
    </xsl:template>
    <xsl:template match="schema">
        <xsl:copy-of select="*"/>
    </xsl:template>
</xsl:stylesheet>