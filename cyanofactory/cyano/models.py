from django.db import models
import datetime
import math
from django.utils import timezone
from django.utils.safestring import mark_safe
from django.core.urlresolvers import reverse

from biosql.helpers import get_database_item
from db_xref.views import get_database_url_from_organism

class Feature(object):
    def __init__(self, organism, feature):
        self.feature = feature
        self.organism = organism
        self.seq = feature.extract(self.organism.seq)

    @property
    def name(self):
        return self.feature.qualifiers["gene"][0]

    @property
    def start(self):
        return min(self.feature.location.start, self.feature.location.end)

    @property
    def end(self):
        return max(self.feature.location.start, self.feature.location.end)
    
    @property
    def strand(self):
        return self.feature.strand
    
    @property
    def qualifiers(self, value):
        return self.feature.qualifiers[value]
    
    def get_absolute_url(self):
        return reverse("cyano.views.gene", args=[self.organism.annotations["accessions"][0], self.name])
    
    def has_qualifier(self, typ):
        return typ in self.feature.qualifiers
    
    def __getitem__(self, key):
        return self.feature.qualifiers[key]
    
    def __len__(self):
        if self.seq == None:
            raise Exception()
        return len(self.seq)
    
    def __str__(self):
        return self.seq

class Gene:
    def __init__(self, organism, gene):
        """
        Creates a new Gene object.
        organism: Organism with the gene (Bioentry)
        gene: Name of the gene
        """
        self.organism = organism
        self.gene_name = gene
        self.record = get_database_item(organism)
        self.gene = None
        self.cds = None
        
        for item in self.record.features:
            if item.type == "gene":
                if "gene" in item.qualifiers:
                    if item.qualifiers["gene"][0] == self.gene_name:
                        self.gene = Feature(self.record, item)
            elif item.type == "CDS":
                if "gene" in item.qualifiers:
                    if item.qualifiers["gene"][0] == self.gene_name:
                        self.cds = Feature(self.record, item)        
        
        
        #if self.gene == None:
            # error
        
        self.fieldsets = [#[
            #("Database", {'fields': [self.get_databases_html()]}),
            #("Moep", {"fields": ["aaa", "bbb"]})
            #]
            
            #('Type', ['model_type']),
            ('Name', [
                {'name': 'Name', 'data': gene},
                {'name': 'Synonyms', 'data': lambda : ", ".join(self.gene.qualifiers.get("gene_synonym", [""])[0].split("; "))},
                {'name': 'Cross References', 'data' : "", 'name': 'protein_complexes'},
                ]),
            ('Structure', [
                {'name': "Structure", 'data': self.get_as_html_structure()},
                {'name': 'Sequence', 'data': self.get_as_html_sequence()},
                {'name': 'Synonyms', 'data': lambda : ", ".join(self.gene.qualifiers.get("gene_synonym", [""])[0].split("; "))},
                {'name': 'Cross References', 'data' : "", 'name': 'protein_complexes'},
                ]),
            ('Comments', [
                {"name": "References",
                 "data": mark_safe(self.get_databases_html())}]),
            #('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
            ]
    
    def get_databases_html(self):
        result = ""
        for x in self.gene["db_xref"]:
            result += '<a href="%s">%s</a>' % (get_database_url_from_organism(x), x)
            result += '<br />'
        return result
    
    def get_as_html_sequence(self):
        from cyano.helpers import format_sequence_as_html
        
        #direction = CHOICES_DIRECTION[[x[0] for x in CHOICES_DIRECTION].index(self.direction)][1]        
        
        #return 'Chromosome: <a href="%s">%s</a>, Coordinate: %s (nt), Length: %s (nt), Direction: %s, G/C content: %.1f%%, Sequence: %s' % (
        #    self.chromosome.get_absolute_url(), self.chromosome.wid, 
        #    self.coordinate, self.length, direction, 
        #    self.get_gc_content() * 100,
        #    format_sequence_as_html(self.get_sequence()))
        #return 'Sequence: %s' %        
        return 'Sequence: %s' % (format_sequence_as_html(self.gene))

    def get_as_html_structure(self):
        return self._get_as_html_structure(zoom = 1, 
            start_coordinate = self.gene.start - 2500, 
            end_coordinate = self.gene.start + len(self.gene) + 2500, 
            #highlight_wid = [self.wid]
            )

    def _get_as_html_structure(self, start_coordinate = None, end_coordinate = None, highlight_wid = None, zoom = 0):
        if zoom == 0:
            return self.get_as_html_structure_global()
        else:
            return self.get_as_html_structure_local(start_coordinate = start_coordinate, end_coordinate = end_coordinate, highlight_wid = highlight_wid)
    
    def get_as_html_structure_global(self):
        ntPerSegment = 1e4
        segmentHeight = 27
        geneHeight = 10
        featureHeight = 5
        nSegments = int(math.ceil(self.length / ntPerSegment))
        W = 636
        H = segmentHeight * nSegments
        segmentLeft = 35
        segmentW = W - 4 - segmentLeft
        chrTop = -12
        
        #style
        colors = ['3d80b3', '3db34a', 'cd0a0a', 'e78f08']
        style = ''
        for i in range(len(colors)):
            style += '.color-%s{fill:#%s; stroke: #%s;}' % (i, colors[i], colors[i], )
    
        #chromosome
        chrStyle = '\
            .chr text{fill:#222; text-anchor:end; alignment-baseline:middle; font-size:10px}\
            .chr line{stroke:#666; stroke-width:0.5px;}\
        '
        
        chr = ''
        for i in range(nSegments):
            x1 = segmentLeft
            x2 = segmentLeft + ((min(self.length, (i+1) * ntPerSegment) - 1) % ntPerSegment) / ntPerSegment * segmentW
            y = chrTop + (i + 1) * segmentHeight
            chr += '<text x="%s" y="%s">%d</text>' % (segmentLeft - 2, y, i * ntPerSegment + 1)
            chr += '<line x1="%s" x2="%s" y1="%s" y2="%s"/>' % (x1, x2, y, y)
        
        #genes
        geneStyle = '\
            .genes g polygon{stroke-width:1px; fill-opacity:0.5;}\
            .genes g text{text-anchor:middle; alignment-baseline:middle; font-size:8px; fill: #222}\
        '
            
        genes = ''        
        genesList = self.genes.all()
        nTus = 0
        iTUs = {}
        tus = []
        for i in range(len(genesList)):
            gene = genesList[i]
            tu = gene.transcription_units.all()[0]
            if iTUs.has_key(tu.wid):
                iTu = iTUs[tu.wid]
            else:
                tus.append(tu)
                iTu = nTus
                iTUs[tu.wid] = iTu
                nTus += 1
                
            iSegment = math.floor((gene.coordinate - 1) / ntPerSegment)
            
            if gene.direction == 'f':
                x1 = segmentLeft + ((gene.coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
                x3 = min(segmentLeft + segmentW, x1 + gene.length / ntPerSegment * segmentW)
                x2 = max(x1, x3 - 5)
            else:
                x3 = segmentLeft + ((gene.coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
                x1 = min(segmentLeft + segmentW, x3 + gene.length / ntPerSegment * segmentW)
                x2 = min(x1, x3 + 5)
                
            y2 = chrTop + (iSegment + 1) * segmentHeight - 2
            y1 = y2 - geneHeight
            
            if math.fabs(x3 - x1) > len(gene.wid) * 5:
                label = gene.wid
            else:
                label = ''
                
            if gene.name:
                tip_title = gene.name
            else:
                tip_title = gene.wid
            tip_content = 'Transcription unit: %s' % tu.name
            tip_title = tip_title.replace("'", "\'")
            tip_content = tip_content.replace("'", "\'")
                
            genes += '<g>\
                <a xlink:href="%s">\
                    <polygon class="color-%s" points="%s,%s %s,%s %s,%s %s,%s %s,%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);"/>\
                </a>\
                <a xlink:href="%s">\
                    <text x="%s" y="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);">%s</text>\
                </a>\
                </g>' % (                
                gene.get_absolute_url(),
                iTu % len(colors),
                x1, y1,
                x2, y1,
                x3, (y1 + y2) / 2, 
                x2, y2,
                x1, y2,
                tip_title, tip_content,
                gene.get_absolute_url(),
                (x1 + x3) / 2, (y1 + y2) / 2 + 1, 
                tip_title, tip_content,
                label,                
                )
        
        #promoters
        promoterStyle = '.promoters rect{fill:#%s; opacity:0.5}' % (colors[0], )
        tfSiteStyle = '.tfSites rect{fill:#%s; opacity:0.5}' % (colors[1], )
        promoters = ''        
        tfSites = ''
        for tu in tus:
            if tu.promoter_35_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_35_coordinate
                tu_length = tu.promoter_35_length
                
                iSegment = math.floor((tu_coordinate - 1) / ntPerSegment)
            
                x = segmentLeft + ((tu_coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
                w = max(1, x1 - min(segmentLeft + segmentW, x1 + tu_length / ntPerSegment * segmentW))
                
                y = chrTop + (iSegment + 1) * segmentHeight + 2
                    
                if tu.name:
                    tip_title = tu.name
                else:
                    tip_title = tu.wid
                tip_title = tip_title.replace("'", "\'")
                    
                promoters += '<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);"/></a>' % (
                    tu.get_absolute_url(), x, y, w, featureHeight, tip_title, 'Promoter -35 box')
                    
            if tu.promoter_10_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_10_coordinate
                tu_length = tu.promoter_10_length
                
                iSegment = math.floor((tu_coordinate - 1) / ntPerSegment)
            
                x = segmentLeft + ((tu_coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
                w = max(1, x1 - min(segmentLeft + segmentW, x1 + tu_length / ntPerSegment * segmentW))
                
                y = chrTop + (iSegment + 1) * segmentHeight + 2
                
                if tu.name:
                    tip_title = tu.name
                else:
                    tip_title = tu.wid
                tip_title = tip_title.replace("'", "\'")
                    
                promoters += '<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);"/></a>' % (
                    tu.get_absolute_url(), x, y, w, featureHeight, tip_title, 'Promoter -10 box')
    
            for tr in tu.transcriptional_regulations.all():
                if tr.binding_site is not None:    
                    iSegment = math.floor((tr.binding_site.coordinate - 1) / ntPerSegment)
            
                    x = segmentLeft + ((tr.binding_site.coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
                    w = max(1, x1 - min(segmentLeft + segmentW, x1 + tr.binding_site.length / ntPerSegment * segmentW))
                    
                    y = chrTop + (iSegment + 1) * segmentHeight + 2    
                    
                    if tr.name:
                        tip_title = tr.name
                    else:
                        tip_title = tr.wid
                    tip_title = tip_title.replace("'", "\'")
                        
                    tfSites += '<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);"/></a>' % (
                        tr.get_absolute_url(), x, y, w, featureHeight, tip_title, 'Transcription factor binding site')
        
        #features
        featureStyle = '.features rect{fill:#%s;}' % (colors[2], )
        features = ''
        types = {}
        nTypes = 0
        for feature in self.features.all():        
            iSegment = math.floor((feature.coordinate - 1) / ntPerSegment)
            
            x = segmentLeft + ((feature.coordinate - 1) % ntPerSegment) / ntPerSegment * segmentW
            w = max(1, x1 - min(segmentLeft + segmentW, x1 + feature.length / ntPerSegment * segmentW))
                
            y = chrTop + (iSegment + 1) * segmentHeight + 2
                        
            if feature.type.all().count() > 0:
                type = feature.type.all()[0].name
            else:
                type = ''
                
            if feature.name:
                tip_title = feature.name
            else:
                tip_title = feature.wid
            tip_title = tip_title.replace("'", "\'")
            
            features += '<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);"/></a>' % (
                feature.get_absolute_url(), x, y, w, featureHeight, tip_title, type, )
    
        return '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%s" height="%s" viewport="0 0 %s %s"><style>%s%s%s%s%s%s</style><g class="chr">%s</g><g class="genes">%s</g><g class="promoters">%s</g><g class="tfSites">%s</g><g class="features">%s</g></svg>' % (
            W, H, W, H, style, chrStyle, geneStyle, promoterStyle, tfSiteStyle, featureStyle, chr, genes, promoters, tfSites, features)
 
    def get_as_html_structure_local(self, start_coordinate = None, end_coordinate = None, highlight_wid = None):
        length = end_coordinate - start_coordinate + 1
        
        W = 636        
        geneHeight = 20
        featureHeight = 10
        H = 2 + geneHeight + 2 + 4 + 1 * (2 + featureHeight) + 2
        
        geneY = 2
        chrY = geneY + geneHeight + 4
        promoterY = chrY + 1 + 2
        tfSiteY = promoterY
        featureY = promoterY
        
        #style
        colors = ['3d80b3', '3db34a', 'cd0a0a', 'e78f08']
        style = ''
        for i in range(len(colors)):
            style += '.color-%s{fill:#%s; stroke: #%s;}' % (i, colors[i], colors[i], )
        
        #chromosome
        chrL = 4.5 * len('%s' % start_coordinate) + 4
        chrR = W - 4.5 * len('%s' % end_coordinate) - 2 - 6
        chrW = chrR - chrL        
        
        chrStyle = '\
        .chr text{fill:#222; alignment-baseline:middle; font-size:10px}\
        .chr line{stroke:#666; stroke-width:0.5px;}\
        '
        
        chr = '<text x="%s" y="%s" style="text-anchor:end;">%s</text><line x1="%s" x2="%s" y1="%s" y2="%s" /><text x="%s" y="%s" style="text-anchor:start;">%s</text>' % (
            chrL - 4, chrY, start_coordinate,
            chrL, chrR, chrY, chrY,
            chrR + 2, chrY, end_coordinate)
        
        #genes
        geneStyle = '\
            .genes g polygon{}\
            .genes g text{text-anchor:middle; alignment-baseline:middle; font-size:8px; fill: #222}\
        '            
        genes = ''        
        genesList = self.record.features
        nTus = 0
        iTUs = {}
        tus = []
        for i in range(len(genesList)):
            gene = genesList[i]
            
            if gene.type != "gene":
                continue
            
            if not "gene" in gene.qualifiers:
                continue
            
            gene = Feature(self.record, gene)
            
            #print dir(gene)
            if gene.end > end_coordinate or gene.start + len(gene) - 1 < start_coordinate:
                continue
            
            #tu = gene.transcription_units.all()[0]
            #if iTUs.has_key(tu.wid):
            #    iTu = iTUs[tu.wid]
            #else:
            #    tus.append(tu)
            #    iTu = nTus
            #    iTUs[tu.wid] = iTu
            #    nTus += 1
            
            if gene.strand == 1:
                x1 = chrL + float(gene.start - start_coordinate) / length * chrW
                x3 = chrL + float(gene.start + len(gene) - 1 - start_coordinate) / length * chrW
                x2 = max(x1, x3 - 5)
            else:
                x3 = chrL + float(gene.start - start_coordinate) / length * chrW
                x1 = chrL + float(gene.start + len(gene) - 1 - start_coordinate) / length * chrW
                x2 = min(x1, x3 + 5)
                
            x1 = max(chrL, min(chrR, x1))
            x2 = max(chrL, min(chrR, x2))
            x3 = max(chrL, min(chrR, x3))
                        
            y1 = geneY
            y2 = geneY + geneHeight
            
            #if highlight_wid is None or gene.wid in highlight_wid:
            #    fillOpacity = 0.75
            #    strokeOpacity = 1
            #    strokeWidth = 3
            #else:
            fillOpacity = 0.15
            strokeOpacity = 0.35
            strokeWidth = 1
            
            #if math.fabs(x3 - x1) > len(gene.wid) * 5:
            #    label = gene.wid
            #else:
            #    label = ''
            
            if gene["gene"][0]:
                tip_title = gene["gene"][0]
            else:
                tip_title = gene.wid
            tip_content = 'Transcription unit: %s' % "FIXME" #tu.name
            tip_title = tip_title.replace("'", "\'")
            tip_content = tip_content.replace("'", "\'")
                
            genes += '<g style="">\n\
                <a xlink:href="%s">\n\
                    <polygon class="color-%s" points="%s,%s %s,%s %s,%s %s,%s %s,%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);" style="fill-opacity: %s; stroke-opacity: %s; stroke-width: %spx"/>\n\
                </a>\n\
                <a xlink:href="%s">\n\
                    <text x="%s" y="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);">%s</text>\n\
                </a>\n\
                </g>\n' % (                
                gene.get_absolute_url(),
                10 % len(colors),
                #iTu % len(colors),
                x1, y1,
                x2, y1,
                x3, (y1 + y2) / 2, 
                x2, y2,
                x1, y2,
                tip_title, tip_content,
                fillOpacity, strokeOpacity, strokeWidth,
                gene.get_absolute_url(),
                (x1 + x3) / 2, (y1 + y2) / 2 + 1, 
                tip_title, tip_content,
                gene.name#"FIXMELabel",#label,
                )
                
        #promoters
        promoterStyle = '.promoters rect{fill:#%s; opacity:0.5}' % (colors[0], )
        tfSiteStyle = '.tfSites rect{fill:#%s; opacity:0.5}' % (colors[1], )
        promoters = ''        
        tfSites = ''
        for tu in tus:
            if tu.promoter_35_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_35_coordinate
                tu_length = tu.promoter_35_length
                
                if not (tu_coordinate > end_coordinate or tu_coordinate + tu_length - 1 < start_coordinate):            
                    x1 = chrL + float(tu_coordinate - start_coordinate) / length * chrW
                    x2 = chrL + float(tu_coordinate + tu_length - 1 - start_coordinate) / length * chrW
                    
                    x1 = max(chrL, min(chrR, x1))
                    x2 = max(chrL, min(chrR, x2))
                        
                    if highlight_wid is None or tu.wid in highlight_wid:
                        opacity = 1
                    else:
                        opacity = 0.25
                        
                    if tu.name:
                        tip_title = tu.name
                    else:
                        tip_title = tu.wid
                    tip_title = tip_title.replace("'", "\'")
                        
                    promoters += '<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);" style="opacity: %s;"/></a>' % (
                        tu.get_absolute_url(), x1, promoterY, x2 - x1, featureHeight, tip_title, 'Promoter -35 box', opacity)
                    
            if tu.promoter_10_coordinate is not None:
                tu_coordinate = tu.get_coordinate() + tu.promoter_10_coordinate
                tu_length = tu.promoter_10_length
                
                if not (tu_coordinate > end_coordinate or tu_coordinate + tu_length - 1 < start_coordinate):            
                    x1 = chrL + float(tu_coordinate - start_coordinate) / length * chrW
                    x2 = chrL + float(tu_coordinate + tu_length - 1 - start_coordinate) / length * chrW
                    
                    x1 = max(chrL, min(chrR, x1))
                    x2 = max(chrL, min(chrR, x2))
                        
                    if highlight_wid is None or tu.wid in highlight_wid:
                        opacity = 1
                    else:
                        opacity = 0.25
                        
                    if tu.name:
                        tip_title = tu.name
                    else:
                        tip_title = tu.wid
                    tip_title = tip_title.replace("'", "\'")
                    
                    promoters += '<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);" style="opacity: %s;"/></a>' % (
                        tu.get_absolute_url(), x1, promoterY, x2 - x1, featureHeight, tip_title, 'Promoter -10 box', opacity)

            for tr in tu.transcriptional_regulations.all():
                if tr.binding_site is not None and not (tr.binding_site.coordinate > end_coordinate or tr.binding_site.coordinate + tr.binding_site.length - 1 < start_coordinate):    
                    x1 = chrL + float(tr.binding_site.coordinate - start_coordinate) / length * chrW
                    x2 = chrL + float(tr.binding_site.coordinate + tr.binding_site.length - 1 - start_coordinate) / length * chrW
                    
                    x1 = max(chrL, min(chrR, x1))
                    x2 = max(chrL, min(chrR, x2))
                        
                    if highlight_wid is None or tr.wid in highlight_wid:
                        opacity = 1
                    else:
                        opacity = 0.25
                        
                    if tr.name:
                        tip_title = tr.name
                    else:
                        tip_title = tr.wid
                    tip_title = tip_title.replace("'", "\'")
                    
                    tfSites += '<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);" style="opacity: %s;"/></a>' % (
                        tr.get_absolute_url(), x1, promoterY, x2 - x1, featureHeight, tip_title, 'Transcription factor binding site', opacity)
                
        #features
        featureStyle = '.features rect{fill:#%s;}' % (colors[2], )
        features = ''
        types = {}
        nTypes = 0
        for feature in self.record.features: 
            feature = Feature(self.record, feature)       
            if feature.end > end_coordinate or feature.start + len(feature) - 1 < start_coordinate:
                continue
        
            x1 = chrL + float(feature.start - start_coordinate) / length * chrW
            x2 = chrL + float(feature.start + len(feature) - 1 - start_coordinate) / length * chrW
            
            x1 = max(chrL, min(chrR, x1))
            x2 = max(chrL, min(chrR, x2))
                            
            #if feature.type.all().count() > 0:
            #    type = feature.type.all()[0].name
            #else:
            type = ''
                
            #if highlight_wid is None or feature.wid in highlight_wid:
            #    opacity = 1
            #else:
            opacity = 0.25
                
            #if feature.name:
            #    tip_title = feature.name
            #else:
            #    tip_title = feature.wid
            tip_title = feature.name
            tip_title = tip_title.replace("'", "\'")
            
            features += '<a xlink:href="%s"><rect x="%s" y="%s" width="%s" height="%s" onmousemove="javascript: showToolTip(evt, \'%s\', \'%s\')" onmouseout="javascript: hideToolTip(evt);" style="opacity: %s;"/></a>' % (
                feature.get_absolute_url(),
                x1, featureY, x2 - x1, featureHeight, tip_title, type, opacity)    
        
        return '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%s" height="%s" viewport="0 0 %s %s"><style>%s%s%s%s%s%s</style><g class="chr">%s</g><g class="genes">%s</g><g class="promoters">%s</g><g class="tfSites">%s</g><g class="features">%s</g></svg>' % (
            W, H, W, H, style, chrStyle, geneStyle, promoterStyle, tfSiteStyle, featureStyle, chr, genes, promoters, tfSites, features)      
      
        
        
    class Meta:
        fieldsets = ["test", 
        #    ('Type', {'fields': ['model_type']}),
        #    ('Name', {'fields': ['wid', 'name', 'synonyms', 'cross_references']}), 
        #    ('Classification', {'fields': ['type']}), 
        #    ('Content', {'fields': [
        #        {'verbose_name': 'Metabolites (mM)', 'name': 'biomass_compositions'},
        #        {'verbose_name': 'Protein monomers', 'name': 'protein_monomers'},
        #        {'verbose_name': 'Protein complexes', 'name': 'protein_complexes'},
        #        ]}),
        #    ('Comments', {'fields': ['comments', 'references']}),
        #    ('Metadata', {'fields': [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}]}),
        ]       
        field_list = []
        template = "cyano/detail.html"

class Genes():
    class Meta:
        field_list = [
            ["gene", "Gene"],
            ["locus_tag", "Locus Tag"],
            ["function", "Function"],
        ]
        field_url = "gene"
        template = "cyano/list.html"

# Create your models here.
#class Poll(models.Model):
#    question = models.CharField(max_length = 200)
#    pub_date = models.DateTimeField("data published")
#    
#    def was_published_recently(self):
#        return self.pub_date >= timezone.now() - datetime.timedelta(days = 1)
#    was_published_recently.admin_order_field = "pub_date"
#    was_published_recently.boolean = True
#    was_published_recently.short_description = "Published recently?"
#    
#    def __unicode__(self):
#        return self.question

#class Choice(models.Model):
#    poll = models.ForeignKey(Poll)
#    choice = models.CharField(max_length = 200)
#    votes = models.IntegerField()
#    
#    def __unicode__(self):
#        return self.choice
