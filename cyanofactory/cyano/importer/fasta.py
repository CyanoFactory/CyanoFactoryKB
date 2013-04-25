from importer import Importer
from importer import ParserError
from django.core.exceptions import ObjectDoesNotExist
import cyano.models as models
from cyano.helpers import format_field_detail_view, render_queryset_to_response,\
    objectToQuerySet, format_field_detail_view_diff
from django.db.models.related import RelatedObject
from django.template.defaultfilters import capfirst
import re
import StringIO

def diffgen(di):
    old = ""
    for i, t in di:
        t = re.sub(r"(\s{20})", r"\1&#8203;", t)
            
        if i == -1:
            old += '<span style="background: red">' + t + '</span>'
        elif i == 0:
            old += t
        else:
            pass
    new = ""
    for i, t in di:
        t = re.sub(r"(\s{20})", r"\1&#8203;", t)
        if i == -1:
            pass
        elif i == 0:
            new += t
        else:
            new += '<span style="background: green">' + t + '</span>'

    return old, new

def seqdiffgen(old_seq, new_seq):
    from Bio import pairwise2
    pairwise2.MAX_ALIGNMENTS = 1
    align = pairwise2.align.globalxx(old_seq[:100], new_seq[:100])
    output = StringIO.StringIO()
    
    for o, n in zip(align[0][0], align[0][1]):
        output.write("|" if o == n else " ")
    
    align_graphic = output.getvalue()
    output.close()
    output = StringIO.StringIO()
    
    old_align = re.sub(r'(.{80})', r'\1\n', align[0][0]).split("\n")
    new_align = re.sub(r'(.{80})', r'\1\n', align[0][1]).split("\n")
    align_graphic = re.sub(r'(.{80})', r'\1\n', align_graphic).split("\n")
    output.write("""<div class="sequence">""")
    output.write("Alignment of first 500 bp.")
    if old_seq == new_seq:
        output.write(" Sequences are completly identical.")
    else:
        output.write(" Sequences differ.")
    output.write("<div></div><div>")
    for ol, ag, ne in zip(old_align, align_graphic, new_align):
        output.write(ol + "<br/>")
        output.write(ag.replace(" ", "&nbsp;") + "<br/>")
        output.write(ne + "<br/>")
        output.write("<br/>")
    output.write("</div></div>")
    
    return output.getvalue()

class Fasta(Importer):
    def __init__(self):
        super(Importer, self).__init__()
        self.sequences = []
    
    def __parse_header(self, lno, header):
        """Parses a fasta header and returns name and wid based on the | char.
        Detects with and without gi.
        """
        if len(header) == 0 or header[0] != ">":
            raise ParserError(header, lno, "Invalid FASTA Header")
        
        header = header[1:].split("|")
        if len(header) >= 5:
            if header[0] == "gi":
                wid = header[3].strip()
                name = header[4].strip()
            elif header[0] == "ref":
                wid = header[1].strip()
                
            else:
                raise ParserError(header, lno, "Expected gi| in header")
        elif len(header) >= 2:
            wid = header[0].strip()
            name = header[1].strip()
        else:
            raise ParserError(header, lno, "Header needs wid and name delimited by |")
        
        if re.search(r"\s+", wid):
            raise ParserError(header, lno, "wid must not contain whitespace characters")
        
        return wid, name
    
    def load(self, filename):
        seqcount = -1 # To keep wid ordering
        with open(filename, 'r') as f:
            wid = None
            for lno, line in enumerate(f, start = 1):                
                if line[0] == '>':
                    seqcount += 1
                    wid, name = self.__parse_header(lno, line)
                    self.sequences.append({"wid": wid, "name": name, "sequence": ''})
                else:
                    if wid is None:
                        raise ParserError(line, lno, "No FASTA header found")
                    self.sequences[-1]["sequence"] += line.strip()    
        
        # File IO finished
        #self.messages.append(len(self.sequences) + " sequences read")
        
    #===========================================================================
    #    #retrieve chromosomes
    #    for wid, sequence in sequences.iteritems():
    #        try:
    #            chr = Chromosome.objects.get(species__wid=species_wid, wid=wid)
    #        except ObjectDoesNotExist as error:
    #            error_messages.append(error.message)            
    #            continue        
    #    if len(error_messages) > 0:
    #        return (False, format_list_html(error_messages))
    #    
    #    #validate
    #    for wid, sequence in sequences.iteritems():
    #        chr = Chromosome.objects.get(species__wid=species_wid, wid=wid)
    #        chr.sequence = sequence
    #        try:
    #            chr.full_clean()
    #        except ValidationError as error:
    #            error_messages.append(format_error_as_html_list(error))
    #            continue
    #    if len(error_messages) > 0:
    #        return (False, format_list_html(error_messages))
    # 
    #        
    #    #save
    #    for wid, sequence in sequences.iteritems():    
    #        chr = Chromosome.objects.get(species__wid=species_wid, wid=wid)
    #        chr.sequence = sequence
    #        chr.full_clean()
    #        try:
    #            chr.save()
    #        except ValueError as error:
    #            error_messages.append(error.message)
    #            continue
    #    if len(error_messages) > 0:
    #        return (False, format_list_html(error_messages))
    #                
    #    return (True, 'Sequences successfully saved!')
    #===========================================================================
    
    def preview(self, request, species_wid):
        #super(Importer, self).preview(request, species_wid)
        import diff_match_patch as diff
        d = diff.diff_match_patch()
        
        try:
            species = models.Chromosome.objects.get(wid = species_wid)
        except ObjectDoesNotExist:
            species = models.Species(wid = species_wid)
        
        for sequence in self.sequences:
            wid = sequence["wid"]
            name = sequence["name"]
            seq = sequence["sequence"]
            try:
                old_chr = models.Chromosome.objects.get(species__wid = species_wid, wid = wid)
                chr = models.Chromosome.objects.get(species__wid = species_wid, wid = wid)
            except ObjectDoesNotExist:
                old_chr = None
                chr = models.Chromosome(wid = wid)     
            
            chr.name = name
            chr.sequence = seq

        fieldsets = [
            ("Update " + species.wid, {'fields': [
                {'verbose_name': 'WID', 'name': 'wid'},
                {'verbose_name': 'Name', 'name': 'name'},
                {'verbose_name': 'Sequence', 'name': 'sequence'}]})
            ]
        
        #filter out type, metadata
        fieldset_names = [x[0] for x in fieldsets]
        chr.model_type = models.TableMeta.objects.get(name = "Chromosome")
        if 'Type' in fieldset_names:
            idx = fieldset_names.index('Type')
            del fieldsets[idx]
        
        model = models.Chromosome
        #filter out empty fields
        rmfieldsets = []
        for idx in range(len(fieldsets)):
            rmfields = []
            for idx2 in range(len(fieldsets[idx][1]['fields'])):
                if isinstance(fieldsets[idx][1]['fields'][idx2], dict):
                    field_name = fieldsets[idx][1]['fields'][idx2]['name']
                    verbose_name = fieldsets[idx][1]['fields'][idx2]['verbose_name']
                else:
                    field_name = fieldsets[idx][1]['fields'][idx2]
                    field = model._meta.get_field_by_name(field_name)[0]
                    if isinstance(field, RelatedObject):
                        verbose_name = capfirst(field.get_accessor_name())
                        pass
                    else:
                        verbose_name = field.verbose_name
                    
                old_data, new_data = format_field_detail_view_diff(species, old_chr, chr, field_name, request.user.is_anonymous())
                
                if not old_data and not new_data:
                    rmfields = [idx2] + rmfields
                
                fieldsets[idx][1]['fields'][idx2] = {'verbose_name': verbose_name.replace(" ", '&nbsp;').replace("-", "&#8209;"), 'new_data': new_data, 'old_data': old_data}
            for idx2 in rmfields:
                del fieldsets[idx][1]['fields'][idx2]
            if len(fieldsets[idx][1]['fields']) == 0:
                rmfieldsets = [idx] + rmfieldsets
        for idx in rmfieldsets:
            del fieldsets[idx]
        
        return render_queryset_to_response(
            species = species,        
            request = request, 
            models = [model],
            queryset = objectToQuerySet(chr),
            template = 'cyano/preview.html', 
            data = {
                'model_type': "Chromosome",
                'model': model,
                'fieldsets': fieldsets,
                'message': request.GET.get('message', ''),
                })
    
    def submit(self, request, species):
        pass