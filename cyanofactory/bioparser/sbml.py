"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from .bioparser import BioParser
from django.db.transaction import atomic
from libsbml import SBMLReader
from cyano.helpers import slugify
import cyano.models as cmodels

class SBML(BioParser):
    def parse(self, handle):
        if hasattr(self, "notify_progress"):
            self.notify_progress(current = 0, total = 1, message = "Parsing SBML...")
            
        self.reader = SBMLReader()
        self.document = self.reader.readSBMLFromString("".join(line for line in handle))
        self.model = self.document.getModel()
        self.compartments = list(map(lambda i: self.model.getCompartment(i), range(len(self.model.getListOfCompartments()))))
        self.sbml_species = list(map(lambda i: self.model.getSpecies(i), range(len(self.model.getListOfSpecies()))))
        self.reactions = list(map(lambda i: self.model.getReaction(i), range(len(self.model.getListOfReactions()))))

    @atomic
    def apply(self):
        self.detail.save()
        
        total = len(self.compartments) + len(self.sbml_species) + len(self.reactions)

        # Compartment importer
        for i, compartment in enumerate(self.compartments):
            wid = slugify(compartment.getId())
            
            if compartment.getName():
                name = compartment.getName()
            else:
                name = wid

            if hasattr(self, "notify_progress"):
                out_str = "Importing Compartment %s (%d/%d)" % (wid, i + 1, total)
                self.notify_progress(current = i+1, total = total, message = out_str)
            
            # TODO: compartment.getOutside() not implemented
            cobj = cmodels.Compartment.objects.for_species(self.species).for_wid(wid, create = True)
            
            cobj.name = name
            cobj.species = self.species
            cobj.save(self.detail)

        # Species (= Metabolites) importer
        for i, specie in enumerate(self.sbml_species):
            if not self.model.getCompartment(specie.getCompartment()):
                ##self.stderr.write("WARN: Species {} has invalid compartment {}".format(specie.id, specie.getCompartment()))
                continue
                
            wid = slugify(specie.getId())
            if specie.getName():
                name = specie.getName()
            else:
                name = wid
            
            if hasattr(self, "notify_progress"):
                current = len(self.compartments) + i + 1
                out_str = "Importing Metabolite %s (%d/%d)" % (wid, current, total)
                self.notify_progress(current = current, total = total, message = out_str)
            
            # TODO: specie.getBoundaryCondition() not implemented
            sobj = cmodels.Metabolite.objects.for_species(self.species).for_wid(wid, create = True)
            
            sobj.name = name
            sobj.charge = 0 # TODO
            sobj.is_hydrophobic = False # TODO
            sobj.species = self.species
            sobj.save(self.detail)

        for i, reaction in enumerate(self.reactions):
            wid = slugify(reaction.getId())
            if reaction.getName():
                name = reaction.getName()
            else:
                name = wid
                
            valid = False
            
            if hasattr(self, "notify_progress"):
                current = len(self.compartments) + len(self.sbml_species) + i + 1
                out_str = "Importing Reaction %s (%d/%d)" % (wid, current, total)
                self.notify_progress(current = current, total = total, message = out_str)
            
            # Validation of reactants
            reactants = map(lambda i: reaction.getReactant(i), range(len(reaction.getListOfReactants())))
            products = map(lambda i: reaction.getProduct(i), range(len(reaction.getListOfProducts())))
            
            for reactant in reactants:
                if not self.model.getSpecies(reactant.getSpecies()):
                    ##self.stderr.write("WARN: Reactant {} has invalid species {}".format(reactant.id, reactant.species))
                    break
            else:
                # Validation of products
                for product in products:
                    if not self.model.getSpecies(product.getSpecies()):
                        ##self.stderr.write("WARN: Product {} has invalid species {}".format(product.id, product.species))
                        break
                else:
                    # Validation passed
                    valid = True
            
            if valid:
                reaction_obj = cmodels.Reaction.objects.for_species(self.species).for_wid(wid, create = True)
                
                reaction_obj.name = name
                reaction_obj.direction = 'r' if reaction.getReversible() else 'f'
                reaction_obj.is_spontaneous = False  # TODO
                reaction_obj.species = self.species
                reaction_obj.save(self.detail)

                for reactant in reactants:
                    #try:
                    #    participant_obj = cmodels.ReactionStoichiometryParticipant.objects.get(wid = wid)
                    #except ObjectDoesNotExist:
                    #    participant_obj = cmodels.ReactionStoichiometryParticipant(wid = wid)

                    participant_obj = cmodels.ReactionStoichiometryParticipant()
                    participant_obj.molecule = cmodels.Metabolite.objects.for_species(self.species).for_wid(slugify(reactant.getSpecies()))
                    participant_obj.coefficient = -reactant.getStoichiometry()
                    participant_obj.compartment = cmodels.Compartment.objects.for_species(self.species).for_wid(slugify(self.model.getSpecies(reactant.getSpecies()).getCompartment()))
                    participant_obj.save(self.detail)

                    reaction_obj.stoichiometry.add(participant_obj)

                for product in products:
                    #try:
                    #    participant_obj = cmodels.ReactionStoichiometryParticipant.objects.get(wid = wid)
                    #except ObjectDoesNotExist:
                    #    participant_obj = cmodels.ReactionStoichiometryParticipant(wid = wid)

                    participant_obj = cmodels.ReactionStoichiometryParticipant()
                    participant_obj.molecule = cmodels.Metabolite.objects.for_species(self.species).for_wid(slugify(product.getSpecies()))
                    participant_obj.coefficient = product.getStoichiometry()
                    participant_obj.compartment = cmodels.Compartment.objects.for_species(self.species).for_wid(slugify(self.model.getSpecies(product.getSpecies()).getCompartment()))
                    participant_obj.detail = self.detail
                    participant_obj.save(self.detail)

                    reaction_obj.stoichiometry.add(participant_obj)
