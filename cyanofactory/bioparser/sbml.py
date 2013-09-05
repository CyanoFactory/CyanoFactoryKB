from bioparser import BioParser
from django.db.transaction import commit_on_success
from libsbml import SBMLReader
from cyano.helpers import slugify
import cyano.models as cmodels
from django.core.exceptions import ObjectDoesNotExist

class SBML(BioParser):
    def parse(self, handle):
        if hasattr(self, "notify_progress"):
            self.notify_progress(current = 0, total = 1, message = "Parsing SBML...")
            
        self.reader = SBMLReader()
        self.document = self.reader.readSBMLFromString("".join(line for line in handle))
        self.model = self.document.getModel()
        self.compartments = map(lambda i: self.model.getCompartment(i), range(len(self.model.getListOfCompartments())))
        self.sbml_species = map(lambda i: self.model.getSpecies(i), range(len(self.model.getListOfSpecies())))
        self.reactions = map(lambda i: self.model.getReaction(i), range(len(self.model.getListOfReactions())))

    @commit_on_success
    def apply(self):
        self.detail.save()
        
        total = len(self.compartments) + len(self.sbml_species) + len(self.reactions)

        # Compartment importer
        for i, compartment in enumerate(self.compartments):                
            wid = slugify(compartment.id)
            
            if compartment.name:
                name = compartment.name
            else:
                name = wid

            if hasattr(self, "notify_progress"):
                out_str = "Importing Compartment %s (%d/%d)" % (wid, i + 1, total)
                self.notify_progress(current = i+1, total = total, message = out_str)
            
            # TODO: compartment.getOutside() not implemented
            try:
                cobj = cmodels.Compartment.objects.get(species__pk = self.species.pk, wid = wid)
            except ObjectDoesNotExist:
                cobj = cmodels.Compartment(wid = wid)
            
            cobj.name = name
            cobj.save(self.detail)
            cobj.species.add(self.species)

        # Species (= Metabolites) importer
        for i, specie in enumerate(self.sbml_species):
            if not self.model.getCompartment(specie.getCompartment()):
                ##self.stderr.write("WARN: Species {} has invalid compartment {}".format(specie.id, specie.getCompartment()))
                continue
                
            wid = slugify(specie.id)
            if specie.name:
                name = specie.name
            else:
                name = wid
            
            if hasattr(self, "notify_progress"):
                current = len(self.compartments) + i + 1
                out_str = "Importing Metabolite %s (%d/%d)" % (wid, current, total)
                self.notify_progress(current = current, total = total, message = out_str)
            
            # TODO: specie.getBoundaryCondition() not implemented
            try:
                sobj = cmodels.Metabolite.objects.get(species__pk = self.species.pk, wid = wid)
            except ObjectDoesNotExist:
                sobj = cmodels.Metabolite(wid = wid)
            
            sobj.name = name
            sobj.charge = 0 # TODO
            sobj.is_hydrophobic = False # TODO
            sobj.save(self.detail)
            sobj.species.add(self.species)
        
        for i, reaction in enumerate(self.reactions):
            wid = slugify(reaction.id)
            if reaction.name:
                name = reaction.name
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
                if not self.model.getSpecies(reactant.species):
                    ##self.stderr.write("WARN: Reactant {} has invalid species {}".format(reactant.id, reactant.species))
                    break
            else:
                # Validation of products
                for product in products:
                    if not self.model.getSpecies(product.species):
                        ##self.stderr.write("WARN: Product {} has invalid species {}".format(product.id, product.species))
                        break
                else:
                    # Validation passed
                    valid = True
            
            if valid:
                try:
                    reaction_obj = cmodels.Reaction.objects.get(species__pk = self.species.pk, wid = wid)
                except ObjectDoesNotExist:
                    reaction_obj = cmodels.Reaction(wid = wid)
                
                reaction_obj.name = name
                reaction_obj.direction = 'r' if reaction.reversible else 'f'
                reaction_obj.save(self.detail)
                
                for reactant in reactants:
                    #try:
                    #    participant_obj = cmodels.ReactionStoichiometryParticipant.objects.get(wid = wid)
                    #except ObjectDoesNotExist:
                    #    participant_obj = cmodels.ReactionStoichiometryParticipant(wid = wid)
                
                    participant_obj = cmodels.ReactionStoichiometryParticipant()
                    participant_obj.molecule = cmodels.Metabolite.objects.get(wid = slugify(reactant.species))
                    participant_obj.coefficient = -reactant.stoichiometry
                    participant_obj.compartment = cmodels.Compartment.objects.get(wid = slugify(self.model.getSpecies(reactant.species).compartment))
                    # TODO: EvidencedData needs Revisioning
                    participant_obj.detail = self.detail
                    participant_obj.save()
                    
                    reaction_obj.stoichiometry.add(participant_obj)
                
                for product in products:
                    #try:
                    #    participant_obj = cmodels.ReactionStoichiometryParticipant.objects.get(wid = wid)
                    #except ObjectDoesNotExist:
                    #    participant_obj = cmodels.ReactionStoichiometryParticipant(wid = wid)
                
                    participant_obj = cmodels.ReactionStoichiometryParticipant()
                    participant_obj.molecule = cmodels.Metabolite.objects.get(wid = slugify(product.species))
                    participant_obj.coefficient = product.stoichiometry
                    participant_obj.compartment = cmodels.Compartment.objects.get(wid = slugify(self.model.getSpecies(product.species).compartment))
                    # TODO: EvidencedData needs Revisioning
                    participant_obj.detail = self.detail
                    participant_obj.save()
            
                    reaction_obj.stoichiometry.add(participant_obj)
                
                reaction_obj.save(self.detail)
                reaction_obj.species.add(self.species)
