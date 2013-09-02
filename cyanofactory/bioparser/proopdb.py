# Imports Operon Prediction files from ProOpDB
# Reference:
# Taboada B., Ciria R., Martinez-Guerrer C.E., Merino E., 2012, ProOpDB:
# Prokaryotic Operon DataBase, Nucleic Acids Research, 40(D1), D627-D631

from django.core.exceptions import ObjectDoesNotExist

import cyano.models as cmodels
from bioparser import BioParser
from django.db.transaction import commit_on_success
    
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
        tus = []
    
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
                    gene = cmodels.Gene.objects.get(wid = locus)
                except ObjectDoesNotExist:
                    #self.stderr.write("Gene with wid {} not found".format(locus))
                    continue
                
                tu = operon_str
                
                try:
                    tus[self.tu.index(tu)].append(gene)
                except ValueError:
                    tus.append([tu, [gene]])
        
        for tu, genes in tus:
            wid = "TU_" + "-".join(map(lambda x: x.wid, genes))
            try:
                tu = cmodels.TranscriptionUnit.objects.get(wid = wid)
            except ObjectDoesNotExist:
                tu = cmodels.TranscriptionUnit(wid = wid)
            self.tu.append([tu, genes])

    @commit_on_success
    def apply(self):
        self.detail.save()

        i = 0
        
        for tu, genes in self.tu:
            i += 1
            
            if hasattr(self, "notify_progress"):
                outstr = "Importing TU %s (%d/%d)" % (tu.wid, i, self.count)
                self.notify_progress(current = i, total = self.count, message = outstr)
            
            tu.save(self.detail)
            tu.species.add(self.species)

            for gene in genes:
                tu.genes.add(gene)
