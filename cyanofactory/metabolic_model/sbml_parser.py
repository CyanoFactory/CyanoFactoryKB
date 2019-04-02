"""
Copyright (c) 2018 Gabriel Kind <kind hs-mittweida de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

# geneProductAssociation JSON
# Helper anpassen auf neue Klasse
# js-code anpassen (Kotz!), vllt. erstmal so belassen, bis Paper durch? Kann man danach noch immer ändern


from xml.sax import make_parser, handler
from io import StringIO

import json
import sys

from .metabolic_model import *
from .sbml_xml_generator import SbmlXMLGenerator
from .sbml_json_generator import SbmlJsonGenerator
from .optgene import OptGeneParser

class SbmlElement(object):
    def __init__(self, name, attrs = None):
        name_split = name.split(":")

        self.ns = name_split[0] if len(name_split) > 1 else ""
        self.name = name_split[0] if len(self.ns) == 0 else name_split[1]
        self.complete_name = name
        self.attrs = attrs

    def __repr__(self):
        return self.complete_name


class SbmlHandlerBase(handler.ContentHandler):
    def __init__(self):
        self.elements = []  # type: List[SbmlElement]

        super().__init__()

        #print("Start of %s" % self)

    def startElement(self, name, attrs):
        self.elements.append(SbmlElement(name, attrs))
        #print("Push: %s" % name)

        self.start(self.elements[-1])

    def start(self, element: SbmlElement):
        pass

    def endElement(self, name):
        if len(self.elements) == 0:
            #print("End of %s | %s" % (self, name))
            self.end(SbmlElement(name))
            pop_handler()
            handlers[-1].endElement(name)
            return

        elem = self.elements.pop()
        #print("Pop: %s" % elem.complete_name)
        if elem.complete_name != name:
            raise ValueError("Tag Mismatch (expected %s, got %s" % (name, elem.complete_name))

        self.end(SbmlElement(name))

    def end(self, element):
        pass

    def common_read_func(self, lst, obj, element):
        lst.append(obj)
        lst[-1].read_attributes(element.attrs)


class SbmlHandler(SbmlHandlerBase):
    def __init__(self):
        self.model = MetabolicModel()
        super().__init__()

    def start(self, element):
        if element.name == "model":
            self.model.read_attributes(element.attrs)
        elif element.name == "listOfCompartments":  # done
            push_handler(CompartmentHandler(self.model.compartments))
        elif element.name == "listOfSpecies":  # done
            push_handler(MetaboliteHandler(self.model.metabolites))
        elif element.name == "listOfObjectives":
            push_handler(ObjectiveHandler(self.model.objectives))
        elif element.name == "listOfGeneProducts":
            push_handler(GeneProductHandler(self.model.genes))
        elif element.name == "listOfGroups":
            push_handler(GroupHandler(self.model.groups))
        elif element.name == "listOfParameters":
            push_handler(ParameterHandler(self.model.parameters))
        elif element.name == "listOfReactions":
            push_handler(ReactionHandler(self.model.reactions))
        elif element.name == "listOfUnitDefinitions":
            push_handler(UnitDefinitionHandler(self.model.unit_definitions))
        elif len(self.elements) >= 2 and self.elements[-2].name == "model" and element.name == "annotation":
            push_handler(AnnotationHandler(self.model))

    def end(self, element):
        if element.name == "sbml":
            a = 0


class CompartmentHandler(SbmlHandlerBase):
    def __init__(self, compartments):
        self.compartments = compartments
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "compartment":
            self.common_read_func(self.compartments, Compartment(), element)
        elif self.elements[-2].name == "compartment" and name == "annotation":
            push_handler(AnnotationHandler(self.compartments[-1]))


class MetaboliteHandler(SbmlHandlerBase):
    def __init__(self, metabolites):
        self.metabolites = metabolites
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "species":
            self.common_read_func(self.metabolites, Metabolite(), element)
        elif self.elements[-2].name == "species" and name == "annotation":
            push_handler(AnnotationHandler(self.metabolites[-1]))


class ObjectiveHandler(SbmlHandlerBase):
    def __init__(self, objectives):
        self.objectives = objectives
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "objective":
            self.common_read_func(self.objectives, Objective(), element)
        elif self.elements[-2].name == "objective" and name == "annotation":
            push_handler(AnnotationHandler(self.objectives[-1]))
        elif name == "listOfFluxObjectives":
            push_handler(FluxObjectiveHandler(self.objectives[-1].flux_objectives))


class FluxObjectiveHandler(SbmlHandlerBase):
    def __init__(self, flux_objectives):
        self.flux_objectives = flux_objectives
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "fluxObjective":
            self.common_read_func(self.flux_objectives, FluxObjective(), element)
        elif self.elements[-2].name == "fluxObjective" and name == "annotation":
            push_handler(AnnotationHandler(self.flux_objectives[-1]))


class GeneProductHandler(SbmlHandlerBase):
    def __init__(self, genes):
        self.genes = genes
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "geneProduct":
            self.common_read_func(self.genes, GeneProduct(), element)
        elif self.elements[-2].name == "geneProduct" and name == "annotation":
            push_handler(AnnotationHandler(self.genes[-1]))


class GroupHandler(SbmlHandlerBase):
    def __init__(self, groups):
        self.groups = groups
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "group":
            self.common_read_func(self.groups, Group(), element)
        elif self.elements[-2].name == "group" and name == "annotation":
            push_handler(AnnotationHandler(self.groups[-1]))
        elif self.elements[-2].name == "group" and name == "listOfMembers":
            push_handler(GroupMemberHandler(self.groups[-1].members))


class GroupMemberHandler(SbmlHandlerBase):
    def __init__(self, members):
        self.members = members
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "member":
            self.common_read_func(self.members, GroupMember(), element)


class ParameterHandler(SbmlHandlerBase):
    def __init__(self, parameters):
        self.parameters = parameters
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "parameter":
            self.common_read_func(self.parameters, Parameter(), element)


class ReactionHandler(SbmlHandlerBase):
    def __init__(self, reactions):
        self.reactions = reactions
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "reaction":
            self.common_read_func(self.reactions, Reaction(), element)
        elif self.elements[-2].name == "reaction" and name == "annotation":
            push_handler(AnnotationHandler(self.reactions[-1]))
        elif self.elements[-2].name == "reaction" and name == "listOfReactants":
            push_handler(MetaboliteReferenceHandler(self.reactions[-1].substrates))
        elif self.elements[-2].name == "reaction" and name == "listOfProducts":
            push_handler(MetaboliteReferenceHandler(self.reactions[-1].products))
        elif self.elements[-2].name == "reaction" and name == "geneProductAssociation":
            push_handler(GeneProductAssociationHandler(self.reactions[-1].gene_products))


class MetaboliteReferenceHandler(SbmlHandlerBase):
    def __init__(self, ref):
        self.ref = ref
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "speciesReference":
            self.common_read_func(self.ref, MetaboliteReference(), element)


class GeneProductAssociationHandler(SbmlHandlerBase):
    def __init__(self, products):
        self.products = products
        super().__init__()

    def start(self, element):
        prod_ref = GeneProductReference()
        name = element.name
        if name in ["and", "or"]:
            prod_ref.logical = name
            self.common_read_func(self.products, prod_ref, element)
            push_handler(GeneProductAssociationHandler(self.products[-1].children))
        elif name == "geneProductRef":
            self.common_read_func(self.products, prod_ref, element)


class UnitDefinitionHandler(SbmlHandlerBase):
    def __init__(self, unit_definitions):
        self.unit_definitions = unit_definitions
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "unitDefinition":
            self.common_read_func(self.unit_definitions, UnitDefinition(), element)
        elif self.elements[-2].name == "unitDefinition" and name == "annotation":
            push_handler(AnnotationHandler(self.unit_definitions[-1]))
        elif self.elements[-2].name == "unitDefinition" and name == "listOfUnits":
            push_handler(UnitHandler(self.unit_definitions[-1].units))


class UnitHandler(SbmlHandlerBase):
    def __init__(self, units):
        self.units = units
        super().__init__()

    def start(self, element):
        name = element.name
        if name == "unit":
            self.common_read_func(self.units, Unit(), element)
        elif self.elements[-2].name == "unit" and name == "annotation":
            push_handler(AnnotationHandler(self.units[-1]))


class AnnotationHandler(SbmlHandlerBase):
    def __init__(self, parent_obj: ElementBase):
        self.parent = parent_obj
        super().__init__()

    def start(self, element):
        if element.ns == "bqbiol" or element.ns == "bqmodel":
            # Verify element stack
            if self.elements[-2].name == "Description" and self.elements[-3].name == "RDF":
                annotation = Annotation()
                annotation.ns = element.ns
                annotation.type = self.elements[-1].name
                if self.parent.description is None:
                    self.parent.description = Description()
                self.parent.description.annotation.append(annotation)
                push_handler(RdfBagHandler(annotation))
            else:
                raise ValueError("Bad bqbiol entry")


class RdfBagHandler(SbmlHandlerBase):
    def __init__(self, annotation: Annotation):
        self.annotation = annotation
        super().__init__()

    def start(self, element):
        if element.name == "li" and "rdf:resource" in element.attrs:
            if self.elements[-2].name != "Bag":
                raise ValueError("Expected bag, got %s" % element.name)
            self.annotation.resources.append(element.attrs["rdf:resource"])


handlers = []

def push_handler(handle):
    handlers.append(handle)
    parser.setContentHandler(handle)


def pop_handler():
    old = handlers.pop()

    if len(handlers) > 0:
        parser.setContentHandler(handlers[-1])
    else:
        a = 0
        pass

    return old

parser = make_parser()
