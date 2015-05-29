#    Copyright (C) 2012 Daniel Gamermann <daniel.gamermann@ucv.es>
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
#    Please, cite us in your reasearch!
#


class EnzError(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return self.value




class Enzyme:
    '''This class defines a chemical reaction object. The input should be an   
       string containing a reaction in OptGene format.
       '''

    def __init__(self,linea, pathway="GP"):
        parts = linea.split(" : ")
        name = parts[0].replace(" ", "")
        #print linea
        self.pathway = pathway
        self.name = name
        self.issues = False
        self.issues_info = []
        try:
            reac = parts[1]
        except IndexError:
            raise EnzError("Badly written reaction: %s"%linea)
        reac = reac.replace('\r', '')
        reac = reac.replace('\n', '')
        if "<->" in reac and (" <-> " not in reac):
            reac = reac.replace("<->", " <->")
        if " <-> " not in reac:
            reacion = reac.split(" -> ")
            self.reversible = False
            sus = reacion[0]
            try:
                pro = reacion[1]
            except IndexError:
                raise EnzError("Badly written reaction: %s"%linea)
        else:
            reacion = reac.split(" <-> ")
            self.reversible = True
            sus = reacion[0]
            try:
                pro = reacion[1]
            except IndexError:
                raise EnzError("Badly written reaction: %s"%linea)
        suss = sus.split(" + ")
        pros = pro.split(" + ")
        substrates = []
        products = []
        stoicsust = []
        stoicprod = []
        for ele in suss:
            mets = ele.split(" ")
            topop = []
            for i, met in enumerate(mets):
                if met == "": topop.append(i)
            topop.sort(reverse=True)
            for ii in topop:
                mets.pop(ii)
            if len(mets) == 0:
                self.issues = True
                self.issues_info.append("No substrate")
            elif len(mets) == 1:
                stoicsust.append(1.)
                substrates.append(mets[0])
            else:
                try:
                    stoic = float(mets[0])
                    na = "".join(mets[1:])
                    stoicsust.append(stoic)
                    substrates.append(na)
                except:
                    na = "".join(mets)
                    stoicsust.append(1.)
                    substrates.append(na)
        for ele in pros:
            mets = ele.split(" ")
            topop = []
            for i, met in enumerate(mets):
                if met == "": topop.append(i)
            topop.sort(reverse=True)
            for ii in topop:
                mets.pop(ii)
            if len(mets) == 0:
                self.issues = True
                self.issues_info.append("No product")
            elif len(mets) == 1:
                stoicprod.append(1.)
                products.append(mets[0])
            else:
                try:
                    stoic = float(mets[0])
                    na = "".join(mets[1:])
                    stoicprod.append(stoic)
                    products.append(na)
                except:
                    na = "".join(mets)
                    stoicprod.append(1.)
                    products.append(na)
        for su in substrates:
            for pr in products:
                if su == pr:
                    self.issues = True
                    self.issues_info.append("Metabolite " + su + 
                                            " connected to itself")
        # checks repeated metabolites
        if len(substrates) != len(list(set(substrates))):
            subs = list(set(substrates))
            stoic_sus = []
            for ii, su in enumerate(subs):
                inds = [kk for kk in xrange(len(substrates)) if substrates[kk] == su]
                coefs = [stoicsust[ele] for ele in inds]
                ncoef = sum(coefs)
                stoic_sus.append(ncoef)
        else:
            subs = substrates
            stoic_sus = stoicsust
        if len(products) != len(list(set(products))):
            pros = list(set(products))
            stoic_pro = []
            for ii, pr in enumerate(pros):
                inds = [kk for kk in xrange(len(products)) if products[kk] == pr]
                coefs = [stoicprod[ele] for ele in inds]
                ncoef = sum(coefs)
                stoic_pro.append(ncoef)
        else:
            pros = products
            stoic_pro = stoicprod
        self.substrates = subs
        self.products = pros
        self.metabolites = subs + pros
        nsusu = len(subs)
        nprod = len(pros)
        self.Nsubstrates = nsusu
        self.Nproducts = nprod
        self.stoic = [stoic_sus,stoic_pro]
        if self.reversible:
            if nsusu > nprod:
                self.tup = (nprod,nsusu)
            else:
                self.tup = (nsusu,nprod)
        else:
            self.tup = (nsusu,nprod)

    def __repr__(self):
        """
        Returns the reaction in OptGene format.
        """
        return str(self)

    def __str__(self):
        """
        Returns the reaction in OptGene format.
        """
        name = self.name
        sustr = self.substrates
        prods = self.products
        stoic = self.stoic
        rev = self.reversible
        stri = name + " : "
        for i in range(len(sustr)):
            stri += str(stoic[0][i]) + " " + sustr[i] + " + "
        if len(sustr) > 0:
            stri = stri[:len(stri)-2]
        if rev:
            stri += " <-> "
        else:
            stri += " -> "
        for i in range(len(prods)):
            stri += str(stoic[1][i]) + " " + prods[i] + " + "
        if self.Nproducts > 0:
            stri = stri[:len(stri)-2]
        return stri

    def __rmul__(self, num, nname=""):
        """
        Multiplication by number: Multiplies all stoichiometric coefficients 
        by the number.
        """
        if not nname:
            nname = self.name + "*"+str(num)
        if num > 0:
            nre = self.copy()
        elif num < 0:
            nre = -self.copy()
            num = -num
        elif num == 0:
            raise EnzError("One should not multiply an reaction by zero!")
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
                to_add.append(str(ncoef) + " " + sub)
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
                to_add.append(str(ncoef) + " " + prod)
        stri_reac += " + ".join(to_add)
        enz = Enzyme(stri_reac)
        if len(mets_el):
            enz.issues = True
            enz.issues_info.append("Possibly eliminated: "+", ".join(mets_el))
        if len(mets_cata):
            enz.issues = True
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
        if len(self.metabolites) != len (other.metabolites):
            return False
        if (self.reversible ^ other.reversible):
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
        raise EnzError(" Comparison went very wrong!")

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
        ''' Checks if the two metabolites given as a list or tupple in mets are
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
        kk=str(self)
        return Enzyme(kk.replace("<->", "->"))

    def make_rev(self):
        ''' Returns the reversible version of the reaction. '''
        kk=str(self)
        return Enzyme(kk.replace("->", "<->"))

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
            stri += str(stoic[0][i]) + " " + sustr[i] + " + "
        if len(sustr) > 0:
            stri = stri[:len(stri)-2]
        if rev:
            stri += " <-> "
        else:
            stri += " -> "
        for i in range(len(prods)):
            stri += str(stoic[1][i]) + " " + prods[i] + " + "
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
            self.metabolites.pop(self.metabolites.index(met))
            self.stoic[0].pop(ii)
            self.Nsubstrates -= 1
            nsusu = self.Nsubstrates
            nprod = self.Nproducts
            if self.reversible:
                if nsusu > nprod:
                    self.tup = (nprod, nsusu)
                else:
                    self.tup = (nsusu, nprod)
            else:
                self.tup = (nsusu, nprod)
        elif met in self.products:
            ii = self.products.index(met)
            self.products.pop(ii)
            self.metabolites.pop(self.metabolites.index(met))
            self.stoic[1].pop(ii)
            self.Nproducts-=1
            nsusu=self.Nsubstrates
            nprod=self.Nproducts
            if self.reversible:
                if nsusu > nprod:
                    self.tup = (nprod, nsusu)
                else:
                    self.tup = (nsusu, nprod)
            else:
                self.tup = (nsusu, nprod)
        if self.Nsubstrates == 0 or self.Nproducts == 0:
            self.issues = True
            self.issues_info.append("lost all sustrates or products after \
                                      a pop")

    def copy(self):
        ''' Returns a copy of the reaction. '''
        bbb = str(self)
        return Enzyme(bbb)

    def stoic_n(self, met):
        ''' Returns the stoichiometric coefficient of the metabolite met. '''
        if self.has_substrate(met):
            return self.stoic[0][self.substrates.index(met)]
        elif self.has_product(met):
            return self.stoic[1][self.products.index(met)]
        else:
            raise EnzError("Metabolite " + met + " not in the reaction!")

