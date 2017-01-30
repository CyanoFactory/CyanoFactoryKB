"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from collections import OrderedDict

from PyNetMet2.metabolism import Metabolism
from PyNetMet2.enzyme import Enzyme

class JsonModel(object):
    class Compound(object):
        def __init__(self, **kwargs):
            self.name = None
            self.stoichiometry = 1

            self.__dict__.update(**kwargs)

            if "stoichiometry" in kwargs:
                num = float(kwargs["stoichiometry"])

                if num.is_integer():
                    self.stoichiometry = int(num)

        def to_json(self):
            out = OrderedDict()
            order = [
                "name",
                "stoichiometry",
            ]
            for i in order:
                val = self.__dict__[i]
                if val:
                    out.update({i: self.__dict__[i]})
            return out

        def __str__(self):
            return u"{} {}".format(self.stoichiometry, self.name)

    class Reaction(object):
        def __init__(self, **kwargs):
            self.name = None
            self.substrates = []
            self.products = []
            self.constraints = []
            self.reversible = False
            self.pathway = None
            self.disabled = False
            self.favourite = False

            self.__dict__.update(**kwargs)
            if "substrates" in kwargs:
                self.substrates = list(map(lambda x: JsonModel.Compound(**x), kwargs["substrates"]))
            if "products" in kwargs:
                self.products = list(map(lambda x: JsonModel.Compound(**x), kwargs["products"]))

            if len(self.constraints) not in [0, 2]:
                raise ValueError("Bad constraints")

        def to_json(self):
            out = OrderedDict()
            order = [
                "name",
                "substrates",
                "products",
                "constraints",
                "reversible",
                "pathway",
                "disabled",
                "favourite"
            ]
            for i in order:
                val = self.__dict__[i]
                if val:
                    if i in ["substrates", "products"]:
                        out.update({i: [x.to_json() for x in self.__dict__[i]]})
                    else:
                        out.update({i: self.__dict__[i]})
            return out

        def __str__(self):
            return u"{} : {} {} {}".format(
                self.name,
                u" + ".join([str(x) for x in self.substrates]),
                u"<->" if self.reversible else "->",
                u" + ".join([str(x) for x in self.products])
            )

    class Metabolite(object):
        def __init__(self, **kwargs):
            self.name = None
            self.external = False

            self.__dict__.update(**kwargs)

        def to_json(self):
            out = OrderedDict()
            order = [
                "name",
                "external",
            ]
            for i in order:
                val = self.__dict__[i]
                if val:
                    out.update({i: self.__dict__[i]})
            return out

        def __str__(self):
            return self.name

    class Objective(object):
        def __init__(self, **kwargs):
            self.name = None
            self.maximize = False

            self.__dict__.update(**kwargs)

        def to_json(self):
            out = OrderedDict()
            order = [
                "name",
                "maximize",
            ]
            for i in order:
                val = self.__dict__[i]
                out.update({i: self.__dict__[i]})
            return out

        def __str__(self):
            return u"{} {}".format(self.name, 1 if self.maximize else -1)

    def __init__(self):
        self.reactions = []
        self.metabolites = []
        self.objectives = []
        self.design_objectives = []
        self.target_reactions = []

    @staticmethod
    def __index_predicate(iterator, pred):
        for i, it in enumerate(iterator):
            if pred(it):
                return i

        raise ValueError("Predicate not matching any in list")

    @staticmethod
    def from_model(model):
        """

        :type model: PyNetMet2.metabolism.Metabolism
        """
        this = JsonModel()

        for enz in model.enzymes:
            reac = JsonModel.Reaction()
            if enz.pathway == "_TRANSPORT_":
                continue
            reac.name = enz.name
            reac.substrates = []
            for s, m in zip(enz.stoic[0], enz.substrates):
                reac.substrates.append(
                    JsonModel.Compound(name=m, stoichiometry=s)
                )
            for s, m in zip(enz.stoic[1], enz.products):
                reac.products.append(
                    JsonModel.Compound(name=m, stoichiometry=s)
                )
            if None not in enz.constraint:
                reac.constraints = enz.constraint
            reac.reversible = enz.reversible
            reac.pathway = enz.pathway

            this.reactions.append(reac)

        for metabolite in model.mets:
            this.metabolites.append(
                JsonModel.Metabolite(name=metabolite,
                                     external=metabolite in model.external)
            )

        for obj in model.obj:
            try:
                this.objectives.append(
                    JsonModel.Objective(name=obj[0],
                                        maximize=int(float(obj[1])) == 1)
                )
            except ValueError:
                raise ValueError("Bad objective")

        for obj in model.design_obj:
            try:
                this.design_objectives.append(
                    JsonModel.Objective(name=obj[0],
                                        maximize=int(float(obj[1])) == 1)
                )
            except ValueError:
                raise ValueError("Bad objective")

        return this

    def to_model(self):
        reactions = []
        external = []
        objective = []
        design_objective = []

        for reac in self.reactions:
            if reac.disabled:
                continue

            enz = Enzyme(reac)
            enz.pathway = reac.pathway

            reactions.append(enz)

            if reac.constraints:
                enz.constraint = (float(reac.constraints[0]), float(reac.constraints[1]))

        for metabolite in self.metabolites:
            if metabolite.external:
                external.append(metabolite.name)

        for obj in self.objectives:
            objective.append(u"{} 1 {}".format(
                obj.name,
                "1" if obj.maximize else "-1"
            ))

        for obj in self.design_objectives:
            design_objective.append(u"{} 1 {}".format(
                obj.name,
                "1" if obj.maximize else "-1"
            ))

        m = Metabolism(
            filein="",
            fromfile=False,
            reactions=[],
            constraints=[],
            external=external,
            objective=objective,
            design_objective=design_objective
        )
        m.enzymes = reactions
        m.calcs()
        return m

    @staticmethod
    def from_json(jstr):
        jm = JsonModel()

        if "reactions" in jstr:
            jm.reactions = list(JsonModel.Reaction(**x) for x in jstr["reactions"])

        if "metabolites" in jstr:
            jm.metabolites = list(JsonModel.Metabolite(**x) for x in jstr["metabolites"])

        if "objectives" in jstr:
            jm.objectives = list(JsonModel.Objective(**x) for x in jstr["objectives"])

        if "design_objectives" in jstr:
            jm.design_objectives = list(JsonModel.Objective(**x) for x in jstr["design_objectives"])

        if "target_reactions" in jstr:
            jm.target_reactions = list(JsonModel.Objective(**x) for x in jstr["target_reactions"])

        return jm

    def to_json(self):
        out = OrderedDict()
        order = [
            "reactions",
            "metabolites",
            "objectives",
            "design_objectives",
            "target_reactions"
        ]
        for i in order:
            val = self.__dict__[i]
            if val:
                out.update({i: [x.to_json() for x in self.__dict__[i]]})
        return out

    def add_reaction(self, reaction):
        """
        Add reaction string to model
        """
        if self.has_reaction(reaction.name):
            raise ValueError("Reaction already in model: " + reaction.name)

        self.reactions.append(reaction)

    def get_reaction(self, name):
        """
        Get reaction with name
        """
        try:
            return self.reactions[JsonModel.__index_predicate(self.reactions, lambda x: x.name == name)]
        except ValueError:
            return None

    def rename_reaction(self, reac_from, reac_to):
        reac = self.get_reaction(reac_from)

        if reac is None:
            raise ValueError("Reaction not in model: " + reac_from)

        if self.has_reaction(reac_to):
            raise ValueError("Reaction already in model: " + reac_to)

        reac.name = reac_to

    def remove_reaction(self, name):
        """
        Remove reaction with passed name from model
        """
        reac = self.get_reaction(name)

        if reac is None:
            raise ValueError("Reaction not in model: " + name)

        del self.reactions[self.reactions.index(reac)]

    def has_reaction(self, name):
        return self.get_reaction(name) is not None

    def add_metabolite(self, name, external):
        JsonModel.check_metabolite_name(name)

        if self.has_metabolite(name):
            raise ValueError("Metabolite already in model: " + name)

        self.metabolites.append(JsonModel.Metabolite(name=name, external=external))

    def get_metabolite(self, name):
        """
        Get reaction with name
        """
        try:
            return self.metabolites[JsonModel.__index_predicate(self.metabolites, lambda x: x.name == name)]
        except ValueError:
            return None

    def has_metabolite(self, name):
        return self.get_metabolite(name) is not None

    def remove_metabolite(self, name):
        from itertools import chain

        for reaction in self.reactions:
            for metabolite in chain(reaction.substrates, reaction.products):
                if metabolite.name == name:
                    raise ValueError("Metabolite in use: " + name)

        self.metabolites = list(filter(lambda x: x.name != name, self.metabolites))

    @staticmethod
    def check_metabolite_name(name):
        import re

        if re.match(r"^[0-9]", name):
            raise ValueError("Name must not begin with a number")

        if re.match(r"(^<?\-> | <?\-> | <?\->$)", name):
            raise ValueError("Name must not contain lonely <-> or ->")

        if re.match(r"(^\+ | \+ )", name):
            raise ValueError("Lonely + only allowed at end of name")

    def rename_metabolite(self, met_from, met_to):
        from itertools import chain

        JsonModel.check_metabolite_name(met_to)

        met = self.get_metabolite(met_from)

        if met is None:
            raise ValueError("Metabolite not in model: " + met_from)

        if self.has_metabolite(met_to):
            raise ValueError("Metabolite already in model: " + met_to)

        met.name = met_to

        for reaction in self.reactions:
            for metabolite in chain(reaction.substrates, reaction.products):
                if metabolite.name == met_from:
                    metabolite.name = met_to

    def set_metabolite_external(self, name, external):
        met = self.get_metabolite(name)

        if met is None:
            raise ValueError("Metabolite not in model: " + name)

        met.external = external
