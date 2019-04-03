"""
Copyright (c) 2018 Gabriel Kind <kind hs-mittweida de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import json
from .sbml_json_generator import SbmlJsonGenerator
from typing import List

class ElementBase(object):
    def __init__(self):
        self.id = ""
        self.metaid = ""
        self.name = ""
        self.sbo_term = ""
        self.description = None

    def __repr__(self):
        def name_and_id():
            if self.name:
                return "%s (%s)" % (self.name, self.id)
            else:
                return self.id

        return "%s: %s" % (type(self).__name__, name_and_id())

    @staticmethod
    def bool(v):
        return not str(v).lower() in ["false", "0"]

    def read_common_attribute(self, k, v):
        if k == "id":
            self.id = v
        elif k == "metaid":
            self.metaid = v
        elif k == "name":
            self.name = v
        elif k == "sboTerm":
            self.sbo_term = v

    def make_common_attributes(self):
        return ["id", self.id,
                "metaid", self.metaid,
                "name", self.name,
                "sboTerm", self.sbo_term]

    def make_common_attributes_ns(self):
        return ["", "id", self.id,
                "", "metaid", self.metaid,
                "", "name", self.name,
                "", "sboTerm", self.sbo_term]

    def write_sbml_description(self, writer):
        if self.description is None:
            return

        has_biol = False
        has_model = False
        for elem in self.description.annotation:
            if elem.ns == "bqbiol":
                has_biol = True
            elif elem.ns == "bqmodel":
                has_model = True

        if has_biol:
            writer.startPrefixMapping("bqbiol", "http://biomodels.net/biology-qualifiers/")
        if has_model:
            writer.startPrefixMapping("bqmodel", "http://biomodels.net/model-qualifiers/")
        writer.startPrefixMapping("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#")

        with writer.element("annotation"):
            with writer.elementNS(ns="rdf", name="RDF"):
                self.description.write_sbml(writer, self.metaid)

        writer.endPrefixMapping("rdf")
        if has_model:
            writer.endPrefixMapping("bqmodel")
        if has_biol:
            writer.endPrefixMapping("bqbiol")

    def to_json(self, **kwargs):
        return json.dumps(self, cls=SbmlJsonGenerator, **kwargs)

    @staticmethod
    def apply_by_id(lst, id, op):
        for i, item in enumerate(lst):
            if item.id == id:
                return op(lst, item, i)
        return None

    @staticmethod
    def apply_by_name(lst, name, op):
        for i, item in enumerate(lst):
            if item.name == name:
                return op(lst, item, i)
        return None

    @staticmethod
    def apply_by_metaid(lst, metaid, op):
        for i, item in enumerate(lst):
            if item.metaid == metaid:
                return op(lst, item, i)
        return None

    @staticmethod
    def list_parse(j, name, lst, cls):
        jlst = j.get(name)
        if jlst is not None:
            if not isinstance(jlst, list):
                raise ValueError(name + " is not a list")
            for item in jlst:
                lst.append(cls.from_json(item))

class LstOp:
    def __init__(self, lst):
        self.lst = lst

    def add(self, obj):
        if self.has(id=obj.id):
            raise ValueError("{} already in list".format(obj))
        self.lst.append(obj)

    def has(self, **kwargs):
        return self.get() is not None

    def get(self, **kwargs):
        return model_getter(self.lst, **kwargs)

    def remove(self, **kwargs):
        if not self.has(**kwargs):
            raise ValueError("{} not in list".format(kwargs))
        return model_deleter(self.lst, **kwargs)

    def change_id(self, new_id, **kwargs):
        obj = self.get(**kwargs)
        if obj is None:
            raise ValueError("{} not in list".format(kwargs))

        if self.has(new_id, **kwargs):
            raise ValueError("Can't change ID to {}, already in list".format(new_id))

        obj.id = new_id

def model_apply(lst_ref, op, **kwargs):
    if "id" in kwargs:
        return ElementBase.apply_by_id(lst_ref, kwargs["id"], op)
    elif "name" in kwargs:
        return ElementBase.apply_by_name(lst_ref, kwargs["name"], op)
    elif "metaid" in kwargs:
        return ElementBase.apply_by_metaid(lst_ref, kwargs["metaid"], op)
    raise ValueError("bad getter")

def model_getter(lst_ref, **kwargs):
    return model_apply(lst_ref, lambda lst, item, idx: item, **kwargs)

def model_deleter(lst_ref, **kwargs):
    return model_apply(lst_ref, lambda lst, item, idx: lst.pop(idx), **kwargs)


class MetabolicModel(ElementBase):
    def __init__(self):
        self.version = ""
        self.metabolites = []  # type: List[Metabolite]
        self.reactions = []  # type: List[Reaction]
        self.genes = []  # ???
        self.compartments = []  # type: List[Compartment]
        self.objectives = []  # type: List[Objective]
        self.parameters = []  # type: List[Parameter]
        self.groups = []  # type: List[Group]
        self.unit_definitions = []  # type: List[UnitDefinition]

        self.extent_units = ""
        self.strict = False
        self.substance_units = ""
        self.time_units = False

        self.lower_bound_ref = None
        self.upper_bound_ref = None
        self.zero_bound_ref = None

        super().__init__()

    @property
    def reaction(self):
        return LstOp(self.reactions)

    @property
    def metabolite(self):
        return LstOp(self.metabolites)

    @property
    def group(self):
        return LstOp(self.groups)

    @property
    def parameter(self):
        return LstOp(self.parameters)

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "extentUnits":
                self.extent_units = v
            elif k == "fbc:strict":
                self.strict = ElementBase.bool(v)
            elif k == "substanceUnits":
                self.substance_units = v
            elif k == "timeUnits":
                self.time_units = v
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        writer.startPrefixMapping("", "http://www.sbml.org/sbml/level3/version1/core")
        writer.startPrefixMapping("fbc", "http://www.sbml.org/sbml/level3/version1/fbc/version2")
        writer.startPrefixMapping("groups", "http://www.sbml.org/sbml/level3/version1/groups/version1")
        if getattr(writer, "enable_wedesign", False):
            writer.startPrefixMapping("wedesign", "http://cyanofactory.hs-mittweida.de/wedesign/version1")

        attrs = ["fbc", "required", False,
                "groups", "required", False,
                "", "level", 3,
                "", "version", 1,
                "", "sboTerm", self.sbo_term]

        with writer.elementNS(name=["http://www.sbml.org/sbml/level3/version1/core", "sbml"], ns="sbml", attrs=attrs):
            with writer.elementNS(ns="", name="model", attrs=self.make_common_attributes_ns() + [
                "", "extentUnits", self.extent_units,
                "fbc", "strict", self.strict,
                "", "substanceUnits", self.substance_units,
                "", "timeUnits", self.time_units
            ]):
                self.write_sbml_description(writer)
                if len(self.objectives) > 0:
                    with writer.elementNS(ns="fbc", name="listOfObjectives", is_list=True):
                        for objective in self.objectives:
                            objective.write_sbml(writer)
                if len(self.genes) > 0:
                    with writer.elementNS(ns="fbc", name="listOfGeneProducts", is_list=True):
                        for product in self.genes:
                            product.write_sbml(writer)
                if len(self.groups) > 0:
                    with writer.elementNS(ns="groups", name="listOfGroups", is_list=True):
                        for group in self.groups:
                            group.write_sbml(writer)
                with writer.element("listOfUnitDefinitions", is_list=True):
                    for unit_def in self.unit_definitions:
                        unit_def.write_sbml(writer)
                with writer.element("listOfCompartments", is_list=True):
                    for compartment in self.compartments:
                        compartment.write_sbml(writer)
                with writer.element("listOfSpecies", is_list=True):
                    for metabolite in self.metabolites:
                        metabolite.write_sbml(writer)
                with writer.element("listOfParameters", is_list=True):
                    for parameter in self.parameters:
                        parameter.write_sbml(writer)
                with writer.element("listOfReactions", is_list=True):
                    for reaction in self.reactions:
                        reaction.write_sbml(writer)

        if getattr(writer, "enable_wedesign", False):
            writer.endPrefixMapping("wedesign")
        writer.endPrefixMapping("groups")
        writer.endPrefixMapping("fbc")
        writer.endPrefixMapping("")

    def create_default_bounds(self):
        self.lower_bound_ref = self.parameter.get(id="cobra_default_lb")
        if self.lower_bound_ref is None:
            p = Parameter()
            p.id = "cobra_default_lb"
            p.name = "cobra default - lb"
            p.sbo_term = "SBO:0000626"
            p.units = "mmol_per_gDW_per_hr"
            p.value = -100000
            self.lower_bound_ref = p
            self.parameters.append(p)

        self.upper_bound_ref = self.parameter.get(id="cobra_default_ub")
        if self.upper_bound_ref is None:
            p = Parameter()
            p.id = "cobra_default_ub"
            p.name = "cobra default - ub"
            p.sbo_term = "SBO:0000626"
            p.units = "mmol_per_gDW_per_hr"
            p.value = 100000
            self.upper_bound_ref = p
            self.parameters.append(p)

        self.zero_bound_ref = self.parameter.get(id="cobra_0_bound")
        if self.zero_bound_ref is None:
            p = Parameter()
            p.id = "cobra_0_bound"
            p.name = "cobra 0 - bound"
            p.sbo_term = "SBO:0000625"
            p.units = "mmol_per_gDW_per_hr"
            p.value = 0
            self.zero_bound_ref = p
            self.parameters.append(p)

    @staticmethod
    def from_json(j):
        if "sbml" in j:
            j = j["sbml"]
        if "model" in j:
            j = j["model"]

        obj = MetabolicModel()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        ElementBase.list_parse(j, "listOfCompartments", obj.compartments, Compartment)
        ElementBase.list_parse(j, "listOfSpecies", obj.metabolites, Metabolite)
        ElementBase.list_parse(j, "fbc:listOfObjectives", obj.objectives, Objective)
        ElementBase.list_parse(j, "fbc:listOfGeneProducts", obj.genes, GeneProduct)
        ElementBase.list_parse(j, "groups:listOfGroups", obj.groups, Group)
        ElementBase.list_parse(j, "listOfParameters", obj.parameters, Parameter)
        ElementBase.list_parse(j, "listOfReactions", obj.reactions, Reaction)
        ElementBase.list_parse(j, "listOfUnitDefinitions", obj.unit_definitions, UnitDefinition)

        obj.create_default_bounds()
        for reac in obj.reactions:
            reac.update_bounds_from_parameters(obj)

        return obj

    def fba(self, objective=None):
        """
        Simplex algorithm through glpk.
        """
        if objective is None:
            raise ValueError("No objective specified")

        # Create fake reactions for external metabolites
        fba_reactions = list(filter(lambda x: x.enabled, self.reactions[:]))
        for metabolite in filter(lambda x: x.is_external(), self.metabolites):
            reaction = Reaction()
            reaction.id = metabolite.id + "_transp"
            reaction.name = metabolite.name + "_transp"
            reaction.reversible = True
            reaction.update_parameters_from_bounds(self)

            ref = MetaboliteReference()
            ref.id = metabolite.id
            ref.stoichiometry = 1.0
            reaction.products.append(ref)

            fba_reactions.append(reaction)

        # Calculate stoic
        stoic = []
        metabolite_ids = list(x.id for x in self.metabolites)

        for i, reac in enumerate(fba_reactions):
            for substrate in reac.substrates:
                mi = metabolite_ids.index(substrate.id)
                stoic.append((mi, i, -substrate.stoichiometry))

            for product in reac.products:
                mi = metabolite_ids.index(product.id)
                stoic.append((mi, i, product.stoichiometry))

        import glpk
        #from scipy.optimize import linprog

        lp = glpk.LPX()
        lp.name = " FBA SOLUTION "
        lp.rows.add(len(metabolite_ids))
        lp.cols.add(len(fba_reactions))

        for i, met in enumerate(metabolite_ids):
            lp.rows[i].name = met

        for i, reac in enumerate(fba_reactions):
            lp.cols[i].name = reac.id

        #constraints

        for i, reac in enumerate(fba_reactions):
            if reac.lower_bound == self.lower_bound_ref.value:
                lb = None
            else:
                lb = reac.lower_bound

            if reac.upper_bound == self.upper_bound_ref.value:
                ub = None
            else:
                ub = reac.upper_bound

            lp.cols[i].bounds = (lb, ub)

        for n in range(len(metabolite_ids)):
            lp.rows[n].bounds = (0., 0.)

        ###### Objective
        lista = [0. for ele in fba_reactions]

        for i, reac in enumerate(fba_reactions):
            if objective is reac:
                lista[i] = 1.0

        lstoic = list(map(lambda x: [0]*len(fba_reactions), [[]]*len(metabolite_ids)))

        for i, reac in enumerate(fba_reactions):
            for j, substrate in enumerate(reac.substrates):
                mi = metabolite_ids.index(substrate.id)
                lstoic[mi][i] = -substrate.stoichiometry

            for j, product in enumerate(reac.products):
                mi = metabolite_ids.index(product.id)
                lstoic[mi][i] = product.stoichiometry

        #res = linprog(c=lista, A_ub=lstoic, bounds=[reac.constraint for reac in fba_reactions], options={"disp": True})
        #print res

        lp.obj[:] = lista[:]
        lp.obj.maximize = True
        ###### Matrix
        lp.matrix = stoic
        lp.simplex()

        res = SimulationResult()
        for flux, reac in zip(lp.cols, fba_reactions):
            res.results.append(SimulationResult.SingleResult(reac.id, flux.value))

        res.solution = lp.status

        return res


# Single compartment in listOfCompartments
class Compartment(ElementBase):
    def __init__(self):
        self.constant = False
        self.units = ""

        super().__init__()

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "constant":
                self.constant = ElementBase.bool(v)
            elif k == "units":
                self.units = v
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = self.make_common_attributes() + [
                "constant", self.constant,
                "units", self.units]

        with writer.element("compartment", attrs):
            self.write_sbml_description(writer)

    @staticmethod
    def from_json(j):
        obj = Compartment()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        return obj


# Single metabolite in listOfSpecies
class Metabolite(ElementBase):
    def __init__(self):
        self.compartment = ""
        self.charge = 0
        self.formula = ""
        self.constant = False
        self.boundary_condition = False
        self.has_only_substance_units = False

        super().__init__()

    def is_external(self):
        return self.compartment == "e"

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "compartment":
                self.compartment = v
            elif k == "fbc:charge":
                self.charge = int(v)
            elif k == "fbc:chemicalFormula":
                self.formula = v
            elif k == "constant":
                self.constant = ElementBase.bool(v)
            elif k == "boundaryCondition":
                self.boundary_condition = ElementBase.bool(v)
            elif k == "hasOnlySubstanceUnits":
                self.has_only_substance_units = ElementBase.bool(v)
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = self.make_common_attributes_ns() + [
            "", "compartment", self.compartment,
            "fbc", "charge", self.charge,
            "fbc", "chemicalFormula", self.formula,
            "", "constant", self.constant,
            "", "boundaryCondition", self.boundary_condition,
            "", "hasOnlySubstanceUnits", self.has_only_substance_units]

        with writer.elementNS(ns="", name="species", attrs=attrs):
            self.write_sbml_description(writer)

    @staticmethod
    def from_json(j):
        obj = Metabolite()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        return obj


# Single reaction in listOfReactions
class Reaction(ElementBase):
    def __init__(self):
        self.metabolites = {}  # JSON?
        self.lower_bound = None  # only JSON
        self.upper_bound = None  # only JSON
        self.lower_bound_name = ""
        self.upper_bound_name = ""
        self.gene_reaction_rule = ""  # only JSON
        self.subsystem = ""  # only JSON == group
        self.reversible = False
        self.fast = False
        self.substrates = []  # type: List[MetaboliteReference]
        self.products = []  # type: List[MetaboliteReference]
        self.gene_products = []  # type: List[GeneProductReference]
        self.enabled = True

        super().__init__()

    def update_parameters_from_bounds(self, model: MetabolicModel):
        if self.lower_bound is None or self.lower_bound == model.lower_bound_ref.value:
            self.lower_bound_name = model.lower_bound_ref.id
        elif self.lower_bound == 0:
            self.lower_bound_name = model.zero_bound_ref.id
        else:
            p = model.parameter.get(id=self.id + "_lower_bound")
            if p is None:
                p = Parameter()
                p.id = self.id + "_lower_bound"
                p.name = self.name + " lower bound"
                p.sbo_term = "SBO:0000625"
                p.units = "mmol_per_gDW_per_hr"
                model.parameters.append(p)
            p.value = self.lower_bound
            self.lower_bound_name = p.id

        if self.upper_bound is None or self.upper_bound_name == model.upper_bound_ref.value:
            self.upper_bound_name = model.upper_bound_ref.id
        elif self.upper_bound == 0:
            self.upper_bound_name = model.zero_bound_ref.id
        else:
            p = model.parameter.get(id=self.id + "_upper_bound")
            if p is None:
                p = Parameter()
                p.id = self.id + "_upper_bound"
                p.name = self.name + " upper bound"
                p.sbo_term = "SBO:0000625"
                p.units = "mmol_per_gDW_per_hr"
                model.parameters.append(p)
            p.value = self.upper_bound
            self.upper_bound_name = p.id

    def update_bounds_from_parameters(self, model):
        if self.lower_bound_name == model.lower_bound_ref.id:
            self.lower_bound = model.lower_bound_ref.value
        elif self.lower_bound_name == model.zero_bound_ref.id:
            self.lower_bound = 0
        else:
            p = model.parameter.get(id=self.id + "_lower_bound")
            self.lower_bound = p.value

        if self.upper_bound_name == model.upper_bound_ref.id:
            self.upper_bound = model.upper_bound_ref.value
        elif self.upper_bound_name == model.zero_bound_ref.id:
            self.upper_bound = 0
        else:
            p = model.parameter.get(id=self.id + "_upper_bound")
            self.upper_bound = p.value

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "fbc:lowerFluxBound":
                self.lower_bound_name = v
            elif k == "fbc:upperFluxBound":
                self.upper_bound_name = v
            elif k == "reversible":
                self.reversible = ElementBase.bool(v)
            elif k == "fast":
                self.fast = ElementBase.bool(v)
            elif k == "wedesign:enabled":
                self.enabled = ElementBase.bool(v)
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = self.make_common_attributes_ns() + [
            "fbc", "lowerFluxBound", self.lower_bound_name,
            "fbc", "upperFluxBound", self.upper_bound_name,
            "", "reversible", self.reversible,
            "", "fast", self.fast
        ]
        if getattr(writer, "enable_wedesign", False):
            attrs += [
                "wedesign", "enabled", self.enabled
            ]

        with writer.elementNS(ns="", name="reaction", attrs=attrs):
            self.write_sbml_description(writer)
            if len(self.gene_products) > 0:
                with writer.elementNS(ns="fbc", name="geneProductAssociation", is_list=True):
                    for gene_product in self.gene_products:
                        gene_product.write_sbml(writer)
            with writer.element("listOfReactants", is_list=True):
                for substrate in self.substrates:
                    substrate.write_sbml(writer)
            with writer.element("listOfProducts", is_list=True):
                for product in self.products:
                    product.write_sbml(writer)

    @staticmethod
    def from_json(j):
        obj = Reaction()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        ElementBase.list_parse(j, "listOfReactants", obj.substrates, MetaboliteReference)
        ElementBase.list_parse(j, "listOfProducts", obj.products, MetaboliteReference)
        ElementBase.list_parse(j, "fbc:geneProductAssociation", obj.gene_products, GeneProductReference)
        return obj


# Single reference to a metabolite in reaction -> listOfReactants/listOfProducts
class MetaboliteReference(ElementBase):
    def __init__(self):
        self.constant = False
        self.stoichiometry = 1.0

        super().__init__()

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "constant":
                self.constant = ElementBase.bool(v)
            elif k == "stoichiometry":
                self.stoichiometry = float(v)
            elif k == "species":
                self.id = v
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = [
            "species", self.id,
            "metaid", self.metaid,
            "name", self.name,
            "sboTerm", self.sbo_term,
            "constant", self.constant,
            "stoichiometry", self.stoichiometry,
        ]

        with writer.element("speciesReference", attrs):
            self.write_sbml_description(writer)

    @staticmethod
    def from_json(j):
        obj = MetaboliteReference()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        return obj


# Single GeneProduct in fbc:listOfGeneProducts
class GeneProduct(ElementBase):
    def __init__(self):
        self.label = ""

        super().__init__()

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "fbc:id":
                self.id = v
            elif k == "fbc:label":
                self.label = v
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = [
            "", "metaid", self.metaid,
            "", "sboTerm", self.sbo_term,
            "fbc", "id", self.id,
            "fbc", "label", self.label
        ]

        with writer.elementNS(ns="fbc", name="geneProduct", attrs=attrs):
            self.write_sbml_description(writer)

    @staticmethod
    def from_json(j):
        obj = GeneProduct()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        return obj


# Single reference to a GeneProductRef in fbc:geneProductAssociation of reaction
class GeneProductReference(ElementBase):
    def __init__(self):
        self.logical = None
        self.gene_product = ""
        self.children = []  # Filled when "and" or "or"

        super().__init__()

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "fbc:geneProduct":
                self.gene_product = v
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        return

        attrs = self.make_common_attributes_ns() + [
            "fbc", "geneProduct", self.gene_product
        ]

        if self.logical is not None:
            with writer.elementNS(ns="fbc", name=self.logical, attrs=attrs):
                self.write_sbml_description(writer)
                for child in self.children:
                    child.write_sbml(writer)
            a = 0
        else:
            with writer.elementNS(ns="fbc", name="geneProductRef", attrs=attrs):
                self.write_sbml_description(writer)

    @staticmethod
    def from_json(j):
        obj = GeneProductReference()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        return obj


# The rdf:Description tag
class Description(ElementBase):
    def __init__(self):
        self.annotation = []  # type: List[Annotation]

        super().__init__()

    def write_sbml(self, writer, src_id):
        with writer.elementNS(ns="rdf", name="Description", attrs=[
            "rdf", "about", "#" + src_id
        ]):
            for annotation in self.annotation:
                annotation.write_sbml(writer)

    @staticmethod
    def from_json(j):
        if j is None:
            return None

        obj = Description()
        for item in j.items():
            obj.annotation.append(Annotation.from_json(item))

        return obj

# Single bqbiol/bqmodel entry of a rdf:Description
class Annotation(ElementBase):
    def __init__(self):
        self.ns = ""
        self.type = ""
        self.resources = []

        super().__init__()

    def write_sbml(self, writer):
        # Technically not a list, but "Bag" is skipped in JSON writer
        with writer.elementNS(name=self.type, ns=self.ns, is_list=True):
            with writer.elementNS(ns="rdf", name="Bag"):
                for resource in self.resources:
                    with writer.elementNS(ns="rdf", name="li", attrs=[
                        "rdf", "resource", resource
                    ]):
                        pass

    @staticmethod
    def from_json(j):
        annotation = Annotation()
        annotation.ns, annotation.type = j[0].split(":")
        annotation.resources = j[1][:]
        return annotation


# Single fbc:objective in fbc:listOfObjectives
class Objective(ElementBase):
    def __init__(self):
        self.type = ""
        self.flux_objectives = []  # type: List[FluxObjective]

        super().__init__()

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "fbc:id":
                self.id = v
            elif k == "fbc:type":
                self.type = v
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = [
            "", "metaid", self.metaid,
            "", "name", self.name,
            "", "sboTerm", self.sbo_term,
            "fbc", "id", self.id,
            "fbc", "type", self.type
        ]

        with writer.elementNS(ns="fbc", name="objective", attrs=attrs):
            self.write_sbml_description(writer)
            for obj in self.flux_objectives:
                with writer.elementNS(ns="fbc", name="listOfFluxObjectives", is_list=True):
                    obj.write_sbml(writer)

    @staticmethod
    def from_json(j):
        obj = Objective()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        ElementBase.list_parse(j, "fbc:listOfFluxObjectives", obj.flux_objectives, FluxObjective)
        return obj


# Single fbc:fluxObjective in fbc:listOfFluxObjectives
class FluxObjective(ElementBase):
    def __init__(self):
        self.coefficient = 0
        self.reaction = ""

        super().__init__()

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "fbc:coefficient":
                self.coefficient = float(v)
            elif k == "fbc:reaction":
                self.reaction = v
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = self.make_common_attributes_ns() + [
            "fbc", "coefficient", self.coefficient,
            "fbc", "reaction", self.reaction
        ]

        with writer.elementNS(ns="fbc", name="fluxObjective", attrs=attrs):
            self.write_sbml_description(writer)

    @staticmethod
    def from_json(j):
        obj = FluxObjective()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        return obj


# Single parameter in listOfParameters
class Parameter(ElementBase):
    def __init__(self):
        self.constant = False
        self.units = ""
        self.value = 0

        super().__init__()

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "units":
                self.units = v
            elif k == "constant":
                self.constant = ElementBase.bool(v)
            elif k == "value":
                self.value = float(v)
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = self.make_common_attributes() + [
            "units", self.units,
            "constant", self.constant,
            "value", self.value
        ]

        with writer.element("parameter", attrs):
            self.write_sbml_description(writer)

    @staticmethod
    def from_json(j):
        obj = Parameter()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        return obj


# Single groups:group in groups:listOfGroups
class Group(ElementBase):
    def __init__(self):
        self.kind = ""
        self.members = []  # type: List[GroupMember]

        super().__init__()

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "groups:id":
                self.id = v
            elif k == "groups:kind":
                self.kind = v
            elif k == "groups:name":
                self.name = v
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = [
            "groups", "id", self.id,
            "groups", "kind", self.kind,
            "groups", "name", self.name,
            "", "sboTerm", self.sbo_term
        ]

        with writer.elementNS(ns="groups", name="group", attrs=attrs):
            self.write_sbml_description(writer)
            with writer.elementNS(ns="groups", name="listOfMembers", is_list=True):
                for member in self.members:
                    member.write_sbml(writer)

    @staticmethod
    def from_json(j):
        obj = Group()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        ElementBase.list_parse(j, "groups:listOfMembers", obj.members, GroupMember)
        return obj


# Single groups:member in groups:listOfMembers
class GroupMember(ElementBase):
    def __init__(self):
        self.id_ref = ""

        super().__init__()

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "groups:idRef":
                self.id_ref = v
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = self.make_common_attributes_ns() + [
            "groups", "idRef", self.id_ref
        ]

        with writer.elementNS(ns="groups", name="member", attrs=attrs):
            self.write_sbml_description(writer)

    @staticmethod
    def from_json(j):
        obj = GroupMember()
        obj.id_ref = j
        return obj


# Single UnitDefinition in listOfUnitDefinitions
class UnitDefinition(ElementBase):
    def __init__(self):
        self.units = []  # type: List[Unit]

        super().__init__()

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = self.make_common_attributes()

        with writer.element("unitDefinition", attrs):
            self.write_sbml_description(writer)
            with writer.element("listOfUnits", is_list=True):
                for unit in self.units:
                    unit.write_sbml(writer)

    @staticmethod
    def from_json(j):
        obj = UnitDefinition()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        ElementBase.list_parse(j, "listOfUnits", obj.units, Unit)
        return obj


# Single Unit in listOfUnits
class Unit(ElementBase):
    def __init__(self):
        self.exponent = 0
        self.kind = ""
        self.multiplier = 0
        self.scale = 0

        super().__init__()

    def read_attributes(self, attrs):
        for k, v in attrs.items():
            if k == "exponent":
                self.exponent = int(v)
            elif k == "kind":
                self.kind = v
            elif k == "multiplier":
                self.multiplier = int(v)
            elif k == "scale":
                self.scale = int(v)
            else:
                self.read_common_attribute(k, v)

    def write_sbml(self, writer):
        attrs = self.make_common_attributes() + [
            "exponent", self.exponent,
            "kind", self.kind,
            "multiplier", self.multiplier,
            "scale", self.scale
        ]

        with writer.element("unit", attrs):
            self.write_sbml_description(writer)

    @staticmethod
    def from_json(j):
        obj = Unit()
        obj.read_attributes(j)
        obj.description = Description.from_json(j.get("annotation"))
        return obj

class SimulationResult(object):
    class SingleResult(object):
        def __init__(self, reaction, flux):
            self.reaction = reaction
            self.flux = flux

    def __init__(self):
        self.solution = ""
        self.results: List[SimulationResult.SingleResult] = []

    def solution_text(self):
        if self.solution == "opt":
            return "Optimal"
        elif self.solution == "undef":
            return "Undefined"
        elif self.solution == "feas":
            return "Maybe not optimal"
        elif self.solution == "infeas" or self.solution == "nofeas":
            return "Unfeasible"
        elif self.solution == "unbnd":
            return "Unbound (check constraints)"
        else:
            return "Unknown Error"
