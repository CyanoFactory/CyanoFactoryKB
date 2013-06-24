from django.core.management.base import BaseCommand
import cyano.models as cmodels
from cyano.helpers import slugify
from django.core.exceptions import ObjectDoesNotExist
from libsbml import SBMLReader
import sys
#from itertools import imap

class Command(BaseCommand):    
    def handle(self, *args, **options):
        for arg in args:
            reader= SBMLReader()
            document = reader.readSBMLFromFile(arg)
            model = document.getModel()
            compartments = map(lambda i: model.getCompartment(i), range(len(model.getListOfCompartments())))
            species = map(lambda i: model.getSpecies(i), range(len(model.getListOfSpecies())))
            reactions = map(lambda i: model.getReaction(i), range(len(model.getListOfReactions())))
            
            species_obj = cmodels.Species.objects.get(wid = "NC_000911")
            
            revdetail = cmodels.RevisionDetail()
            revdetail.user = cmodels.UserProfile.objects.get(user__username__exact = "gabriel")
            revdetail.reason = "Import SBML {}".format(arg)
            
            # Compartment importer
            for i, compartment in enumerate(compartments):                
                wid = slugify(compartment.id)
                
                if compartment.name:
                    name = compartment.name
                else:
                    name = wid
                
                print("Importing Compartment %s (%d/%d)%10s" % (wid, i + 1, len(compartments), " "))
                
                # TODO: compartment.getOutside() not implemented
                try:
                    cobj = cmodels.Compartment.objects.get(wid = wid)
                except ObjectDoesNotExist:
                    cobj = cmodels.Compartment(wid = wid)
                
                cobj.name = name
                cobj.save(revision_detail = revdetail)
                cobj.species.add(species_obj)
                cobj.save(revision_detail = revdetail)

            # Species (= Metabolites) importer
            for i, specie in enumerate(species):
                if not model.getCompartment(specie.getCompartment()):
                    sys.stderr.write("WARN: Species {} has invalid compartment {}".format(specie.id, specie.getCompartment()))
                    continue
                    
                wid = slugify(specie.id)
                if specie.name:
                    name = specie.name
                else:
                    name = wid
                
                print("Importing Metabolite %s (%d/%d)%10s" % (wid, i + 1, len(species), " "))
                
                # TODO: specie.getBoundaryCondition() not implemented
                try:
                    sobj = cmodels.Metabolite.objects.get(wid = wid)
                except ObjectDoesNotExist:
                    sobj = cmodels.Metabolite(wid = wid)
                
                sobj.name = name
                sobj.charge = 0 # TODO
                sobj.is_hydrophobic = False # TODO
                sobj.save(revision_detail = revdetail)
                sobj.species.add(species_obj)
                sobj.save(revision_detail = revdetail)
            
            for i, reaction in enumerate(reactions):
                wid = slugify(reaction.id)
                if reaction.name:
                    name = reaction.name
                else:
                    name = wid
                    
                valid = False
                
                print("Importing Reaction %s (%d/%d)%10s" % (wid, i + 1, len(reactions), " "))
                
                # Validation of reactants
                reactants = map(lambda i: reaction.getReactant(i), range(len(reaction.getListOfReactants())))
                products = map(lambda i: reaction.getProduct(i), range(len(reaction.getListOfProducts())))
                
                for reactant in reactants:
                    if not model.getSpecies(reactant.species):
                        sys.stderr.write("WARN: Reactant {} has invalid species {}".format(reactant.id, reactant.species))
                        break
                else:
                    # Validation of products
                    for product in products:
                        if not model.getSpecies(product.species):
                            sys.stderr.write("WARN: Product {} has invalid species {}".format(product.id, product.species))
                            break
                    else:
                        # Validation passed
                        valid = True
                
                if valid:
                    try:
                        reaction_obj = cmodels.Reaction.objects.get(wid = wid)
                    except ObjectDoesNotExist:
                        reaction_obj = cmodels.Reaction(wid = wid)
                    
                    reaction_obj.name = name
                    reaction_obj.direction = 'r' if reaction.reversible else 'f'
                    reaction_obj.save(revision_detail = revdetail)
                    
                    for reactant in reactants:
                        #try:
                        #    participant_obj = cmodels.ReactionStoichiometryParticipant.objects.get(wid = wid)
                        #except ObjectDoesNotExist:
                        #    participant_obj = cmodels.ReactionStoichiometryParticipant(wid = wid)
                    
                        participant_obj = cmodels.ReactionStoichiometryParticipant()
                        participant_obj.molecule = cmodels.Metabolite.objects.get(wid = slugify(reactant.species))
                        participant_obj.coefficient = -reactant.stoichiometry
                        participant_obj.compartment = cmodels.Compartment.objects.get(wid = slugify(model.getSpecies(reactant.species).compartment))
                        # TODO: EvidencedData needs Revisioning
                        participant_obj.detail = revdetail
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
                        participant_obj.compartment = cmodels.Compartment.objects.get(wid = slugify(model.getSpecies(product.species).compartment))
                        # TODO: EvidencedData needs Revisioning
                        participant_obj.detail = revdetail
                        participant_obj.save()
                
                        reaction_obj.stoichiometry.add(participant_obj)
                    
                    reaction_obj.save(revision_detail = revdetail)
                    reaction_obj.species.add(species_obj)
                    reaction_obj.save(revision_detail = revdetail)
