"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

# Imports Operon Prediction files from ProOpDB
# Reference:
# Taboada B., Ciria R., Martinez-Guerrer C.E., Merino E., 2012, ProOpDB:
# Prokaryotic Operon DataBase, Nucleic Acids Research, 40(D1), D627-D631

from django.core.exceptions import ObjectDoesNotExist

import cyano.models as cmodels
from .bioparser import BioParser
from django.db.transaction import atomic
from collections import OrderedDict
    
class ProOpDB(BioParser):
    def parse(self, handle):
        if hasattr(self, "notify_progress"):
            self.notify_progress(current = 0, total = 1, message = "Parsing Transcription units...")
        
        self.count = 0
        for line in handle:
            if line.count("\t") == 7:
                # 0: Operon ID
                # 1: GeneName (unused)
                # 2: Locus Name
                # 3: GI (unused)
                # 4: Strand (unused)
                # 5: Left Pos (unused)
                # 6: Right Pos (unused)
                # 7: Gene Position in Operon (unused)
                try:
                    line = line.split("\t")
                    operon = int(line[0])
                    self.count += 1
                except ValueError:
                    continue
        
        handle.seek(0)

        i = 0
        
        self.tu = []
        tus = OrderedDict()
    
        for line in handle:
            if line.count("\t") == 7:
                try:
                    line = line.split("\t")
                    operon = int(line[0])
                    locus = line[2]
                    i += 1
                except ValueError:
                    continue
                
                operon_str = "TU_%04d" % (operon)
                
                try:
                    gene = cmodels.Gene.objects.for_species(self.species).for_wid(locus)
                except ObjectDoesNotExist:
                    #self.stderr.write("Gene with wid {} not found".format(locus))
                    continue
                
                tu = operon_str
                
                if tu in tus:
                    tus[operon_str].append(gene)
                else:
                    tus[operon_str] = [gene]
        
        for k in tus.keys():
            genes = tus[k]
            wid = "TU_" + "-".join(map(lambda x: x.wid, genes[:5]))
            if len(genes) > 5:
                wid += "-et_al"

            tu = cmodels.TranscriptionUnit.objects.for_species(self.species).for_wid(wid, create = True)
            self.tu.append([tu, genes])

    @atomic
    def apply(self):
        
        self.detail.save()

        i = 0
        
        for tu, genes in self.tu:
            i += 1
            
            if hasattr(self, "notify_progress"):
                outstr = "Importing TU %s (%d/%d)" % (tu.wid, i, self.count)
                self.notify_progress(current = i, total = self.count, message = outstr)
            
            #print tu, genes
            
            tu.species = self.species
            tu.save(self.detail)

            for gene in genes:
                tu.genes.add(gene)
