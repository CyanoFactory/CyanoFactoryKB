#    Copyright (C) 2012 Daniel Gamermann <daniel.gamermann@ucv.es>
#    Copyright (C) 2015 Gabriel Kind <gkind@hs-mittweida.de>
#
#    This file is part of PyNetMet
#
#    PyNetMet is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PyNetMet is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with PyNetMet.  If not, see <http://www.gnu.org/licenses/>.
#
#
#    Please, cite us in your research!
#
from .util import python_2_unicode_compatible, create_sid

import libsbml

@python_2_unicode_compatible
class Enzyme(object):
    '''
    This class defines a chemical reaction object. The input should be an
    string containing a reaction in OptGene format.
    '''
    def __init__(self, model, linea=None, reac=None, pathway="GP"):
        from itertools import takewhile
        import re

        self.model = model
        self.transport = pathway == "_TRANSPORT_"

        if reac is None:
            self._sbml = model.createReaction()
        else:
            self._sbml = reac
            if not self.meta_id:
                self.meta_id = self.id
            return

        if linea is None:
            # Default ctor
            return

        # Tokenize
        line = re.split("\s+", str(linea))

        # Entry name
        name = list(takewhile(lambda x: not x.endswith(":"), line))
        line = line[len(name):]
        # takewhile did not take the item with the colon
        if len(line) == 0:
            raise ValueError("Invalid reaction: " + linea)
        name.append(line[0])
        line = line[1:]
        # Remove last char which is the :
        self.name = " ".join(name)[:-1].strip()
        self.id = self.name
        self.meta_id = self.name
        self.issues_info = []
        self.substrates = []
        self.products = []
        self.stoic = [[], []]
        self.constraint = (0,0)
        reversible = None
        self.pathway = pathway

        # 0 reading substrates, 1 reading products
        state = 0

        while len(line) > 0:
            value = 1
            # Parse frac
            try:
                if "/" in line[0]:
                    fl = line[0].split("/")
                    value = float(fl[0]) / float(fl[1])
                else:
                    value = float(line[0])
                line = line[1:]
            except ValueError:
                pass  # identifier

            # Identifier, take until operator is found
            metabolite = list(takewhile(lambda x: x != "+" and x != "->" and x != "<->", line))
            line = line[len(metabolite):]
            # Take all following operators
            # If it's a + followed by an arrow it belongs to the name
            ops = list(takewhile(lambda x: x == "+" or x == "->" or x == "<->", line))
            line = line[len(ops):]
            if len(ops) > 0:
                if "->" in ops[:-1] or "<->" in ops[:-1]:
                    raise ValueError("Lonely -> or <-> not allowed in reaction names: " + " ".join(ops[:-1]))
                metabolite += ops[:-1]
                # Push real operator back in
                line.insert(0, ops[-1])

            metabolite = " ".join(metabolite).strip()

            if len(metabolite) > 0:
                if state == 0:
                    target = self.substrate_list
                else:
                    target = self.product_list

                # reference to reactants or products
                # Check for duplicate metabolites first
                for i, t in enumerate(target):
                    if t.getSpecies() == metabolite:
                        target[i].setStoichiometry(target[i].getStoichiometry() + value)
                        if target[i].getStoichiometry() == 0:
                            self.issues_info.append("Eliminated metabolite: " + t)
                            target[i].removeFromParentAndDelete()
                        break
                else:
                    # not found
                    sr = libsbml.SpeciesReference(3,1)
                    sid = create_sid(metabolite)
                    sr.setSpecies(sid)
                    sr.setMetaId(sid)
                    sr.setStoichiometry(value)

                    if not model.getSpecies(sid):
                        species = model.createSpecies()
                        species.setMetaId(sid)
                        species.setId(sid)
                        species.setName(metabolite)

                    target.append(sr)

            if len(line) == 0:
                break

            # First line item is now on an operator
            if line[0] == "->" or line[0] == "<->":
                if reversible is not None:
                    raise ValueError("More then one reaction arrow: " + linea)

                self.reversible = line[0] == "<->"
                reversible = self.reversible
                # Products
                state = 1

            # Skip the operator
            line = line[1:]

        if len(self.substrates) == 0:
            self.issues_info.append("No substrates")

        if len(self.products) == 0:
            self.issues_info.append("No products")

        for su in self.substrates:
            for pr in self.products:
                if su == pr:
                    self.issues_info.append("Metabolite " + su +
                                            " is both substrate and product")

        if reversible is None:
            raise ValueError("No reaction arrow found: " + linea)

        if len(self.substrates) == 0 and self.reversible:
            self.constraint = (None, None)
        elif len(self.substrates) == 0 and not self.reversible:
            self.constraint = (0., None)
        elif len(self.products) == 0 and self.reversible:
            self.constraint = (None, None)
        elif len(self.products) == 0 and not self.reversible:
            self.constraint = (0., None)
        elif self.reversible:
            self.constraint = (None, None)
        else:
            self.constraint = (0., None)

        if pathway:
            self.pathway = pathway

    def get_name_or_id(self):
        return self.name or self.meta_id

    @property
    def meta_id(self):
        return self._sbml.getMetaId()

    @meta_id.setter
    def meta_id(self, value):
        self._sbml.setMetaId(create_sid(value))

    @property
    def id(self):
        return self._sbml.getId()

    @id.setter
    def id(self, value):
        self._sbml.setId(create_sid(value))

    @property
    def name(self):
        return self._sbml.getName()

    @name.setter
    def name(self, value):
        if len(value) > 0 and value.strip()[-1] == ":":
            raise ValueError(": not allowed at end of name")
        self._sbml.setName(value)

    @property
    def reversible(self):
        return self._sbml.getReversible()

    @reversible.setter
    def reversible(self, value):
        self._sbml.setReversible(value)

    @property
    def substrates(self):
        return list(map(lambda x: x.getSpecies(), self._sbml.getListOfReactants()))

    @substrates.setter
    def substrates(self, value):
        self._sbml.getListOfReactants().clear(True)
        for v in value:
            r = self._sbml.createReactant()
            r.setSpecies(v)

    @property
    def products(self):
        return list(map(lambda x: x.getSpecies(), self._sbml.getListOfProducts()))

    @products.setter
    def products(self, value):
        self._sbml.getListOfProducts().clear(True)
        for v in value:
            r = self._sbml.createProduct()
            r.setSpecies(v)

    @property
    def stoic(self):
        return [
            list(map(lambda x: x.getStoichiometry(), self._sbml.getListOfReactants())),
            list(map(lambda x: x.getStoichiometry(), self._sbml.getListOfProducts()))
        ]

    @property
    def substrate_list(self):
        return self._sbml.getListOfReactants()

    @property
    def product_list(self):
        return self._sbml.getListOfProducts()

    @property
    def sbml(self):
        return self._sbml

    @stoic.setter
    def stoic(self, value):
        for s, v in zip(self._sbml.getListOfReactants(), value[0]):
            s.setStoichiometry(v)
        for s, v in zip(self._sbml.getListOfProducts(), value[1]):
            s.setStoichiometry(v)

    @property
    def metabolites(self):
        return self.substrates + self.products

    @property
    def issues(self):
        return len(self.issues_info) > 0

    @property
    def constraint(self):
        lb = self._sbml.getPlugin(0).getLowerFluxBound()
        ub = self._sbml.getPlugin(0).getUpperFluxBound()

        lb_changed = False
        ub_changed = False

        if lb == "cobra_default_lb":
            lb = None
            lb_changed = True
        elif lb:
            lb = self.model.getParameter(lb).getValue()
            lb_changed = True

        if ub == "cobra_default_ub":
            ub = None
            ub_changed = True
        elif ub:
            ub = self.model.getParameter(ub).getValue()
            ub_changed = True

        if not (lb_changed or ub_changed):
            print(lb_changed, ub_changed, self.id, self.reversible)
            if len(self.substrates) == 0 and self.reversible:
                return (None, None)
            elif len(self.substrates) == 0 and not self.reversible:
                return (0., None)
            elif len(self.products) == 0 and self.reversible:
                return (None, None)
            elif len(self.products) == 0 and not self.reversible:
                return (0., None)
            elif self.reversible:
                return (None, None)
            else:
                return (0., None)

        return [lb, ub]

    @constraint.setter
    def constraint(self, val):
        if len(val) < 2:
            if len(self.substrates) == 0 and self.reversible:
                self.constraint = (None, None)
            elif len(self.substrates) == 0 and not self.reversible:
                self.constraint = (0., None)
            elif len(self.products) == 0 and self.reversible:
                self.constraint = (None, None)
            elif len(self.products) == 0 and not self.reversible:
                self.constraint = (0., None)
            elif self.reversible:
                self.constraint = (None, None)
            else:
                self.constraint = (0., None)
            return

        lb = val[0]
        ub = val[1]

        if lb is None:
            self._sbml.getPlugin(0).setLowerFluxBound(self.model.getParameter("cobra_default_lb").getId())
        elif lb == 0:
            self._sbml.getPlugin(0).setLowerFluxBound(self.model.getParameter("cobra_0_bound").getId())
        else:
            param = self.model.getParameter(self.meta_id + "_lower_bound")
            if not param:
                param = self.model.createParameter()
                param.setId(self.meta_id + "_lower_bound")
                param.setName(self.name + " lower bound")
                param.setConstant(True)
                param.setUnits("mmol_per_gDW_per_hr")
                param.setValue(lb)
            self._sbml.getPlugin(0).setLowerFluxBound(self.meta_id + "_lower_bound")

        if ub is None:
            self._sbml.getPlugin(0).setUpperFluxBound(self.model.getParameter("cobra_default_ub").getId())
        elif ub == 0:
            self._sbml.getPlugin(0).setUpperFluxBound(self.model.getParameter("cobra_0_bound").getId())
        else:
            param = self.model.getParameter(self.meta_id + "_upper_bound")
            if not param:
                param = self.model.createParameter()
                param.setId(self.meta_id + "_upper_bound")
                param.setName(self.name + " upper bound")
                param.setConstant(True)
                param.setUnits("mmol_per_gDW_per_hr")
                param.setValue(ub)

            self._sbml.getPlugin(0).setUpperFluxBound(self.meta_id + "_upper_bound")

    @property
    def pathway(self):
        if self.transport:
            return "_TRANSPORT_"

        for g in self.model.getPlugin("groups").getListOfGroups():
            m = g.getMemberByIdRef(self.meta_id)
            if m:
                return m.getIdRef()

        return ""

    @pathway.setter
    def pathway(self, val):
        if val == "_TRANSPORT_":
            self.transport = True
            return

        for g in self.model.getPlugin("groups").getListOfGroups():
            m = g.getMemberByIdRef(self.meta_id)
            if m:
                m.removeFromParentAndDelete()
                break

        g = self.model.getPlugin("groups").getGroup(create_sid(val))

        if not g:
            g = self.model.getPlugin("groups").createGroup()
            g.setId(create_sid(val))
            g.setName(val)
            g.setKind("partonomy")

        m = g.createMember()
        m.setIdRef(self.meta_id)

    @property
    def disabled(self):
        return self.model.getParameter("celldesign_" + self.id + "_disabled") is not None

    @disabled.setter
    def disabled(self, val):
        arg = "celldesign_" + self.id + "_disabled"

        param = self.model.getParameter(arg)

        if val:
            if not param:
                param = self.model.createParameter()
            param.setId(arg)
            param.setValue(1)

            param.setUnits(("_" + str(self.constraint[0])+"X"+str(self.constraint[1])).replace(".", "_"))

            self.constraint = [0, 0]
        else:
            if param:
                const = param.getUnits()[1:].replace("_", ".").split("X")
                self.constraint = [float(const[0]), float(const[1])]
                param.removeFromParentAndDelete()

    def __str__(self):
        """
        Returns the reaction in OptGene format.
        """
        def reac_printer(arg):
            stoic, metabolite = arg

            if stoic == 1:
                return metabolite
            else:
                if stoic == int(stoic):
                    stoic = int(stoic)

                return u"%s %s" % (stoic, metabolite)

        return u"{} : {} {} {}".format(
            self.name,
            " + ".join(map(reac_printer, zip(self.stoic[0], self.substrates))),
            "<->" if self.reversible else "->",
            " + ".join(map(reac_printer, zip(self.stoic[1], self.products))))

    def __rmul__(self, num, new_name=""):
        """
        Multiplication by number: Multiplies all stoichiometric coefficients 
        by the number.
        """
        if not nname:
            nname = self.name + "*"+num
        if num > 0:
            nre = self.copy()
        elif num < 0:
            nre = -self.copy()
            num = -num
        elif num == 0:
            raise ValueError("One should not multiply a reaction by zero!")
        nre.name = nname
        nre.stoic[0] = [num*ele for ele in nre.stoic[0]]
        nre.stoic[1] = [num*ele for ele in nre.stoic[1]]
        return nre

    def __add__(self, other, nname=""):
        """
        Sum of reactions. Sums substrates and products. Also eliminates
        superfulous metabolites and identifies catalyzers.
        """
        # Reaction name
        if not nname:
            nname = self.name + "+" + other.name
        stri_reac = nname + " : "
        # Metabolites
        mets_el = list(set(self.products)&set(other.substrates))
        mets_cata = list(set(self.substrates)&set(other.products))
        all_subs = list(set(self.substrates)|set(other.substrates))
        all_prods = list(set(self.products)|set(other.products))
        # Substrates
        to_add = []
        for sub in all_subs:
            stoic_s = 0
            stoic_o = 0
            stoic_el = 0
            stoic_cata = 0
            if sub in self.substrates:
                stoic_s = self.stoic_n(sub)
            if sub in other.substrates:
                stoic_o = other.stoic_n(sub)
            if sub in mets_el:
                stoic_el = self.stoic_n(sub)
            if sub in mets_cata:
                stoic_cata = other.stoic_n(sub)
            ncoef = stoic_s + stoic_o - stoic_el - stoic_cata
            if ncoef > 0 :
                to_add.append(ncoef + " " + sub)
        stri_reac += " + ".join(to_add)
        # Reversibility
        if not self.reversible or not other.reversible:
            rev = " -> "
        else:
            rev = " <-> "
        stri_reac += rev
        # Products
        to_add = []
        for prod in all_prods:
            stoic_s = 0
            stoic_o = 0
            stoic_el = 0
            stoic_cata = 0
            if prod in self.products:
                stoic_s = self.stoic_n(prod)
            if prod in other.products:
                stoic_o = other.stoic_n(prod)
            if prod in mets_el:
                stoic_el = other.stoic_n(prod)
            if prod in mets_cata:
                stoic_cata = self.stoic_n(prod)
            ncoef = stoic_s + stoic_o - stoic_el - stoic_cata
            if ncoef > 0 :
                to_add.append(ncoef + " " + prod)
        stri_reac += " + ".join(to_add)
        enz = Enzyme(stri_reac)
        if len(mets_el):
            enz.issues_info.append("Possibly eliminated: "+", ".join(mets_el))
        if len(mets_cata):
            enz.issues_info.append("Possible catalyzers: "+", ".join(mets_cata))
        return enz

    def __neg__(self, nname=""):
        """
        Return the reaction reversed.
        """
        return self.rev_reac(nname)

    def __sub__(self, other, nname=""):
        """
        Sums the first reaction with the second one reversed.
        """
        return self.__add__(-other, nname)

    def __eq__(self,other):
        """
        Checks if two enzymes do the same reaction.
        (if they connect the same metabolites with the same stoichiometry.)
        """
        if len(self.metabolites) != len(other.metabolites):
            return False
        if self.reversible ^ other.reversible:
            return False

        mets1 = self.metabolites[:]
        mets2 = other.metabolites[:]
        mets1.sort()
        mets2.sort()
        if mets1 != mets2:
            return False
        subs1 = self.substrates[:]
        pros1 = self.products[:]
        subs2 = other.substrates[:]
        pros2 = other.products[:]
        subs1.sort()
        subs2.sort()
        pros1.sort()
        pros2.sort()
        if not self.reversible:
            if subs1 != subs2 or pros1 != pros2:
                return False
            else:
                blac = True
                for ii in range(len(subs1)):
                    blac = blac and (self.stoic_n(subs1[ii]) == \
                                    other.stoic_n(subs2[ii]))
                for ii in range(len(pros1)):
                    blac = blac and (self.stoic_n(pros1[ii]) == \
                                    other.stoic_n(pros2[ii]))
                return blac
        else:
            if subs1 != subs2:
                if pros1 != subs2 or pros2 != subs1:
                    return False
                else:
                    blac = True
                    for ii in range(len(subs1)):
                        blac = blac and (self.stoic_n(subs1[ii]) == \
                                       other.stoic_n(pros2[ii]))
                    for ii in range(len(pros1)):
                        blac = blac and (self.stoic_n(pros1[ii]) == \
                                       other.stoic_n(subs2[ii]))
                    return blac
            else:
                if pros1 != pros2:
                    return False
                else:
                    blac = True
                    for ii in range(len(subs1)):
                        blac = blac and (self.stoic_n(subs1[ii]) == \
                                       other.stoic_n(subs2[ii]))
                    for ii in range(len(pros1)):
                        blac = blac and (self.stoic_n(pros1[ii]) == \
                                       other.stoic_n(pros2[ii]))
                    return blac

    def has_metabolite(self,metabol):
        ''' Checks if the reaction has the metabolite metabol '''
        if metabol in self.metabolites:
            return True
        else:
            return False

    def has_substrate(self,metabol):
        ''' Checks if the reaction has the substrate metabol '''
        if metabol in self.substrates:
            return True
        else:
            return False

    def has_product(self,metabol):
        ''' Checks if the reaction has the product metabol '''
        if metabol in self.products:
            return True
        else:
            return False

    def has_substrate_rev(self,metabol):
        ''' Checks if the reaction has the substrate metabol, taking into
            account that the reaction may be reversible.'''
        if self.reversible and metabol in self.metabolites:
            return True
        elif metabol in self.substrates:
            return True
        else:
            return False

    def has_product_rev(self,metabol):
        ''' Checks if the reaction has the product metabol, taking into account
            that the reaction may be reversible.'''
        if self.reversible and metabol in self.metabolites:
            return True
        elif metabol in self.products:
            return True
        else:
            return False

    def connects(self,mets):
        ''' Checks if the two metabolites given as a list or tuple in mets are
            connected by the reaction. '''
        if self.has_metabolite(mets[0]) and self.has_metabolite(mets[1]):
            if self.reversible:
                if self.has_substrate(mets[0]) and self.has_product(mets[1]):
                    return True
                elif self.has_substrate(mets[1]) and self.has_product(mets[0]):
                    return True
                else:
                    return False
            else:
                if self.has_substrate(mets[0]) and self.has_product(mets[1]):
                    return True
            return False
        else:
            return False

    def make_irr(self):
        ''' Returns the irreversible version of the reaction. '''
        from copy import copy
        e = copy(self)
        e.reversible = False
        return e

    def make_rev(self):
        ''' Returns the reversible version of the reaction. '''
        from copy import copy
        e = copy(self)
        e.reversible = True
        return e

    def rev_reac(self, nname=""):
        ''' Returns the reaction in the reversed order. '''
        if not nname:
            name = self.name + "_R"
        else:
            name = nname
        sustr = [ele for ele in self.products]
        prods = [ele for ele in self.substrates]
        stoic = [ele for ele in self.stoic]
        stoic.reverse()
        rev = self.reversible
        stri = name + " : "
        for i in range(len(sustr)):
            stri += stoic[0][i] + " " + sustr[i] + " + "
        if len(sustr) > 0:
            stri = stri[:len(stri)-2]
        if rev:
            stri += " <-> "
        else:
            stri += " -> "
        for i in range(len(prods)):
            stri += stoic[1][i] + " " + prods[i] + " + "
        if len(prods) > 0:
            stri = stri[:len(stri)-2]
        return Enzyme(stri)

    def pop(self,met):
        ''' Takes out the metabolite met from the reaction.
        (Careful: a real chemical reaction would have its stoichiometry 
        broken by this operation!!!)'''
        if met in self.substrates:
            ii = self.substrates.index(met)
            self.substrates.pop(ii)
            self.stoic[0].pop(ii)
        if met in self.products:
            ii = self.products.index(met)
            self.products.pop(ii)
            self.stoic[1].pop(ii)
        if len(self.substrates) == 0 or len(self.products) == 0:
            self.issues_info.append("lost all sustrates or products after \
                                      a pop")

    def __copy__(self):
        ''' Returns a copy of the reaction. '''
        e = Enzyme()
        e.name = self.name
        e.issues_info = self.issues_info
        e.substrates = self.substrates[:]
        e.products = self.products[:]
        e.stoic = [[], []]
        e.stoic[0] = self.stoic[0][:]
        e.stoic[1] = self.stoic[1][:]
        e.reversible = self.reversible
        e.pathway = self.pathway
        return e

    def stoic_n(self, met):
        ''' Returns the stoichiometric coefficient of the metabolite met. '''
        if self.has_substrate(met):
            return self.stoic[0][self.substrates.index(met)]
        elif self.has_product(met):
            return self.stoic[1][self.products.index(met)]
        else:
            raise ValueError("Metabolite " + met + " not in the reaction!")

    __repr__ = __str__