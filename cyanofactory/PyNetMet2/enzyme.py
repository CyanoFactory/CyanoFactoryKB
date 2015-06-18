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

from __future__ import unicode_literals

class Enzyme(object):
    '''
    This class defines a chemical reaction object. The input should be an
    string containing a reaction in OptGene format.
    '''
    def __init__(self, linea=None, pathway="GP"):
        from itertools import takewhile
        import re

        if linea is None:
            return

        # Tokenize
        line = re.split("\s+", linea)

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
        self.issues_info = []
        self.substrates = []
        self.products = []
        self.stoic = [[], []]
        self.constraint = []
        self.reversible = None
        self.pathway = pathway
        self.disabled = False

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
                    target = self.substrates
                    target_s = self.stoic[0]
                else:
                    target = self.products
                    target_s = self.stoic[1]

                # reference to reactants or products
                # Check for duplicate metabolites first
                for i, t in enumerate(target):
                    if t == metabolite:
                        target_s[i] += value
                        if target_s[i] == 0:
                            self.issues_info.append("Eliminated metabolite: " + t)
                            target.pop(i)
                            target_s.pop(i)
                        break
                else:
                    # not found
                    target.append(metabolite)
                    target_s.append(value)

            if len(line) == 0:
                break

            # First line item is now on an operator
            if line[0] == "->" or line[0] == "<->":
                if self.reversible is not None:
                    raise ValueError("More then one reaction arrow: " + linea)

                self.reversible = line[0] == "<->"
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

        if self.reversible is None:
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

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if value.strip()[-1] == ":":
            raise ValueError(": not allowed at end of name")
        self._name = value

    @property
    def metabolites(self):
        return self.substrates + self.products

    @property
    def issues(self):
        return len(self.issues_info) > 0

    def __repr__(self):
        """
        Returns the reaction in OptGene format.
        """
        return unicode(self).encode("utf-8")

    def __unicode__(self):
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

                return "%s %s" % (stoic, metabolite)

        return "{} : {} {} {}".format(
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
            nname = self.name + "*"+unicode(num)
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
                to_add.append(unicode(ncoef) + " " + sub)
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
                to_add.append(unicode(ncoef) + " " + prod)
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
                for ii in xrange(len(subs1)):
                    blac = blac and (self.stoic_n(subs1[ii]) == \
                                    other.stoic_n(subs2[ii]))
                for ii in xrange(len(pros1)):
                    blac = blac and (self.stoic_n(pros1[ii]) == \
                                    other.stoic_n(pros2[ii]))
                return blac
        else:
            if subs1 != subs2:
                if pros1 != subs2 or pros2 != subs1:
                    return False
                else:
                    blac = True
                    for ii in xrange(len(subs1)):
                        blac = blac and (self.stoic_n(subs1[ii]) == \
                                       other.stoic_n(pros2[ii]))
                    for ii in xrange(len(pros1)):
                        blac = blac and (self.stoic_n(pros1[ii]) == \
                                       other.stoic_n(subs2[ii]))
                    return blac
            else:
                if pros1 != pros2:
                    return False
                else:
                    blac = True
                    for ii in xrange(len(subs1)):
                        blac = blac and (self.stoic_n(subs1[ii]) == \
                                       other.stoic_n(subs2[ii]))
                    for ii in xrange(len(pros1)):
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
            stri += unicode(stoic[0][i]) + " " + sustr[i] + " + "
        if len(sustr) > 0:
            stri = stri[:len(stri)-2]
        if rev:
            stri += " <-> "
        else:
            stri += " -> "
        for i in range(len(prods)):
            stri += unicode(stoic[1][i]) + " " + prods[i] + " + "
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
