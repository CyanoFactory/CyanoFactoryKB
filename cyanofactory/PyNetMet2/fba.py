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

#import glpk

from .util import python_2_unicode_compatible

@python_2_unicode_compatible
class FBA:
    """
    This class perform FBA analysis for a metabolism object.
    """
    def __init__(self, org, eps=1.e-10, maximize = True):
        """
        Prepares the FBA object.
        """
        # External Mets
        self.max = maximize
        self.eps = eps
        external_in = org.external_in
        external_out = org.external_out
        reacs = org.enzymes
        mets = org.mets
        reac_names = [ele.id for ele in reacs]
        # Stoichiometric matrix
        # Matrix
        lista = []
        for j, reac in enumerate(reacs):
            for met in reac.substrates:
                i = org.dic_mets[met]
                MM = -reac.stoic_n(met)
                if MM != 0.:
                   lista.append((i, j, MM)) # (metabolite, reaction, M_coef)
            for met in reac.products:
                i = org.dic_mets[met]
                MM = reac.stoic_n(met)
                if MM != 0.:
                   lista.append((i, j, MM)) # (metabolite, reaction, M_coef)
        # checks and tries to correct repeated indexes in Mstoic
        indexes = [(ele[0],ele[1]) for ele in lista]
        bla = [indexes.count(ele) for ele in indexes]
        topop = []
        for ii, ele in enumerate(bla):
            if ele != 1:
                if ii in topop:
                    pass
                else:
                    new = ii
                    kk = 1
                    ind = indexes[ii]
                    nv = lista[ii][2]
                    while kk < ele:
                        new = indexes.index(ind,new+1)
                        topop.append(new)
                        nv += lista[new][2]
                        kk += 1
                    lista[ii] = (indexes[ii][0], indexes[ii][1], nv)
                    print("Warning: repeated index in Stoichiometric " + \
                         "matrix. Look into reaction %i."%indexes[ii][1] + \
                         " FBA might not be correct.")
        topop.sort(reverse=True)
        for ele in topop:
            lista.pop(ele)
        MM = lista[:]
        #constraints
        constr = list(map(lambda x: x.constraint, reacs))
        # FBA
        self.ext_in = external_in
        self.ext_out = external_out
        self.mets = mets
        self.reac_names = reac_names
        self.reacs = reacs
        self.constr = constr
        self.Mstoic = MM
        self.obj = org.obj
        self.design_obj = org.design_obj
        self.fba()

    def __str__(self, keyz=lambda x:x[1], rev=False):
        fluxes =[(self.lp.cols[ii].name, self.lp.cols[ii].value) for ii in \
                                                           range(len(self.reacs))]
        fluxes = sorted(fluxes, key=keyz, reverse=rev)
        stri = u""
        for ele in fluxes:
            stri += u"%30s ---->  %9.9f  \n" % ele
        stri += u"\nFlux on objective                   : %9.9f \n" % self.lp.obj.value
        fluxes = [ele[1] for ele in fluxes]
        no_flux = len(list(filter(lambda x:abs(x)<self.eps,fluxes)))
        wi_flux = len(fluxes)-no_flux
        stri += u"Reactions with flux    (flux>eps)   : %i \n"%wi_flux
        stri += u"Reactions without flux (flux<eps)   : %i \n"%no_flux
        stri += u"Solution status: " + self.get_status()
        return stri

    def __sub__(self, other, keyz=lambda x:x[3], rev=False):
        """
        Subtracting two fba's ones obtains the difference in flux between the
        two. This should be used for the same metabolism under different
        conditions. The function defined in keyz is the way the results will 
        be ordered. x =(reaction name, flux in self, flux in other, difference)
        """
        common_reacs = list(set(self.reac_names) & set(other.reac_names))
        percent_list = []
        for ele in common_reacs:
            pos_s = self.reac_names.index(ele)
            pos_o = other.reac_names.index(ele)
            flux_s = self.flux[pos_s]
            flux_o = other.flux[pos_o]
            if abs(flux_s) < self.eps:
                if abs(flux_o) < self.eps:
                    diff = 0.
                else:
                    diff = "NA"
            else:
                diff = 100.*abs(flux_s-flux_o)/abs(flux_s)
            percent_list.append([ele, flux_s, flux_o, diff])
        lista = sorted(percent_list, key=keyz, reverse=rev)
        stri = ""
        for ele in lista:
            blac = "%40s  %20s  %20s   |   dif=%25s \n" % (ele[0], str(ele[1]),
                                                      str(ele[2]),str(ele[3]))
            stri += blac
        return stri

    def print_flux(self, ii):
        """
        Prints the flux for one reaction.
        """
        return "%30s ---->  %9.9f  " % (self.lp.cols[ii].name,
                                        self.lp.cols[ii].value)

    def fba(self):
        """
        Simplex algorithm through glpk.
        """
        import glpk

        lp = glpk.LPX()
        lp.name = " FBA SOLUTION "
        lp.rows.add(len(self.mets))
        lp.cols.add(len(self.reacs))
        for ii, ele in enumerate(lp.rows):
            ele.name = self.mets[ii]
        for ii, ele in enumerate(lp.cols):
            ele.name = self.reac_names[ii]
        #constraints
        ii = 0
        for n in range(len(self.reacs)):
            lp.cols[n].bounds = (self.constr[ii][0], self.constr[ii][1])
            ii += 1
        for n in range(len(self.mets)):
            lp.rows[n].bounds = (0., 0.)
        ###### Objective
        lista = [0. for ele in self.reac_names]
        for obj in self.obj:
            iobj = self.reac_names.index(obj[0])
            lista[iobj] = float(obj[1])
        lp.obj[:] = lista[:]
        lp.obj.maximize = self.max
        ###### Matrix
        lp.matrix = self.Mstoic[:]
        lp.simplex()
        flux = [lp.cols[ii].value for ii in range(len(self.reacs))]
        self.lp = lp
        self.flux = flux
        self.Z = lp.obj.value

    def shadow(self, nrea, ele=0.9, relat = True):
        """
        Calculates the derivative of the flux on the objective reaction
        with respect to one flux. If the flux was not zero, it reduces the flux
        by a factor ele. If the flux was zero it tries to force flux in the
        reaction.        
        """
        ZZ = self.Z
        flu = self.flux[nrea]
        c_old = self.lp.cols[nrea].bounds[:]
        if abs(flu) < self.eps:
            c_new = (ele*.0001,ele*.0001)
        else:
            c_new = (ele*flu,ele*flu)
        self.lp.cols[nrea].bounds = c_new[0], c_new[1]
        self.lp.simplex()
        Znew = self.lp.obj.value
        flu_new = self.lp.cols[nrea].value
        if abs(flu_new-flu)<self.eps or abs(Znew-ZZ)<self.eps:
            shdw = 0.
        else:
            if not relat:
                shdw = (Znew-ZZ)/(flu_new-flu)
            else:
                shdw = abs(flu/ZZ)*(Znew-ZZ)/(flu_new-flu)
        if self.lp.status=="nofeas":
            shdw = "NA"
        self.lp.cols[nrea].bounds = c_old[0], c_old[1]
        self.lp.simplex()
        self.flux = [self.lp.cols[ii].value for ii in range(len(self.reacs))]
        self.Z = self.lp.obj.value
        return shdw

    def essential(self, nreac):
        """
        Tests if a reaction is essential for generating flux at the objective.
        """
        ZZ = self.Z
        fluc = self.flux[nreac]
        if abs(fluc) < self.eps:
            return [False, 0.]
        c_old = self.lp.cols[nreac].bounds[:]
        c_new = (0, 0)
        self.lp.cols[nreac].bounds = c_new[0], c_new[1]
        self.lp.simplex()
        Znew = self.lp.obj.value
        self.lp.cols[nreac].bounds = c_old[0], c_old[1]
        self.lp.simplex()
        self.flux = [self.lp.cols[ii].value for ii in range(len(self.reacs))]
        self.Z = self.lp.obj.value
        if abs(Znew) < self.eps:
            return [True, 1.0]
        else:
            ess = False
        red = (ZZ-Znew)/ZZ
        if abs(red) < self.eps:
            red = 0.
        return [ess, red]

    def KO_fba(self):
        """
        Tests production of the main objective and the design objective in single Knock-Outs of all the reactions.
        """
        # Extract the objective and the design objective
        obj = self.obj
        obj_name = obj[0][0]
        design_obj = self.design_obj
        design_obj_name = design_obj[0][0]
        # Start KO test
        results=[]
        for i in range(len(self.reacs)):
            reac=self.reacs[i]
            if reac.pathway == "_TRANSPORT_": # if it's a transport reaction -> pass
                continue

            if reac.name == obj_name: # if it's the objective function -> calculate only the design objective
                Z = 0
                stat1 = 'main objective deleted'
                constraints = self.constr[i]
                self.constr[i]=(0, 0)
                self.obj = design_obj
                self.fba()
                des_Z = self.Z
                stat2 = self.lp.status
            else:
                constraints = self.constr[i]
                self.constr[i]=(0, 0)
                self.fba()
                Z = self.Z
                stat1 = self.lp.status
                if reac.name == design_obj_name: # if it's the design objective function -> pass the second part
                    des_Z = 0
                    stat2 = 'design objective deleted'
                else:
                    constr = self.constr[self.reac_names.index(obj_name)]
                    self.obj = design_obj
                    self.constr[self.reac_names.index(obj_name)] = (Z, Z)
                    self.fba()
                    des_Z= self.Z
                    stat2 = self.lp.status
                    self.constr[self.reac_names.index(obj_name)] = constr
            results.append([reac.name, round(Z,9), stat1, round(des_Z,9), stat2])
            self.constr[i] = constraints
            self.obj = obj
        return results

    def design_fba(self, units=None):
        saved_constr = []
        saved_obj = self.obj[:]
        all_flux = []

        for obj in self.obj:
            iobj = self.reac_names.index(obj[0])
            saved_constr.append(self.constr[iobj][:])

        if units is None:
            # Determine max and min of obj by constraining obj from 0 to nothing
            self.fba()

            max_flux = self.flux[self.reac_names.index(self.obj[0][0])]

            self.obj = self.design_obj[:]

            units = []

            for i in range(0, 9):
                units.append(round(max_flux * (i / 9), 2))

            units.append(round(max_flux, 2))

        self.obj = self.design_obj[:]

        for unit in units:
            for obj in saved_obj:
                iobj = self.reac_names.index(obj[0])
                self.constr[iobj] = [unit, unit]

            self.fba()

            if self.lp.status != "opt":
                continue

            all_flux.append([unit, self.flux])

        # Restore previous settings
        self.obj = saved_obj
        for i, obj in enumerate(self.obj):
            iobj = self.reac_names.index(obj[0])
            self.constr[iobj] = saved_constr[i]

        return all_flux

    def max_min(self, ii, fixobj=0.5):
        """
        Checks the maximum and minimum flux in each reaction if set as objective
        for a fixed value of the original flux in the objective function.
        """
        # Fixes flux in objective
        obj_old = [ele for ele in self.lp.obj[:]]
        max_old = self.lp.obj.maximize
        cs_old = []
        for obj in self.obj:
            iobj = self.reac_names.index(obj[0])
            cs_old.append(self.constr[iobj])
            fluc = self.flux[iobj]
            nflux = fluc*fixobj
            self.lp.cols[iobj].bounds = nflux, nflux
        nobj = [0.0 for ele in self.reac_names]
        nobj[ii] = 1.0
        self.lp.obj[:] = nobj[:]
        # In direct direction
        self.lp.cols[ii].bounds = (0.0, None)
        self.lp.obj.maximize = False
        self.lp.simplex()
        s1 = self.lp.status
        mini1 = self.lp.obj.value
        self.lp.obj.maximize = True
        self.lp.simplex()
        s2 = self.lp.status
        maxi1 = self.lp.obj.value
        # In reverse direction
        if self.reacs[ii].reversible:
            self.lp.cols[ii].bounds = (None, 0.0)
            self.lp.obj.maximize = False
            self.lp.simplex()
            s3 = self.lp.status
            mini2 = self.lp.obj.value
            self.lp.obj.maximize = True
            self.lp.simplex()
            s4 = self.lp.status
            maxi2 = self.lp.obj.value
        else:
            s3, s4 = 'na', 'na'
            mini2 = 0.
            maxi2 = 0.
        # everything back:
        self.lp.cols[ii].bounds = self.constr[ii]
        self.lp.obj.maximize = max_old
        for jj, obj in enumerate(self.obj):
            self.lp.cols[iobj].bounds = cs_old[jj]
        self.lp.obj[:] = obj_old[:]
        self.lp.simplex()
        # returning results:
        if mini1 < self.eps: mini1 = 0.
        if maxi1 < self.eps: maxi1 = 0.
        if abs(mini2) < self.eps: mini2 = 0.
        if abs(maxi2) < self.eps: maxi2 = 0.
        if s1 != "opt":
            mini1 = "X"
        if s2 != "opt":
            maxi1 = "X"
        if s3 != "opt":
            mini2 = "X"
        if s4 != "opt":
            maxi2 = "X"
        return [(mini1,maxi1),(mini2,maxi2)]

    def check_flux(self):
        """
        Returns the sum of all fluxes over all metabolites.
        This result should be zero!
        """
        reacs_by_met = [[] for ele in self.mets]
        ii = 0
        for ele in self.mets:
            jj = 0
            for reac in self.reacs:
                if reac.has_metabolite(ele):
                    reacs_by_met[ii].append(jj)
                jj += 1
            ii += 1
        self.reacs_by_met = reacs_by_met
        summ = 0.
        for ii in range(len(self.mets)):
            summ += abs(self.sum_flux_met(ii))
        return summ

    def sum_flux_met(self,ii):
        """
        Sums all fluxes in one metabolite.
        This result should be zero!
        """
        met = self.mets[ii]
        reacs = [self.reacs[jj] for jj in self.reacs_by_met[ii]]
        nreacs = len(reacs)
        cond = []
        stoic = []
        for ele in reacs:
            if ele.has_substrate(met):
                cond.append(1.)
                pos = ele.substrates.index(met)
                stoic.append(ele.stoic[0][pos])
            elif ele.has_product(met):
                cond.append(-1.)
                pos = ele.products.index(met)
                stoic.append(ele.stoic[1][pos])
            else:
                print("HELP HELP HELP!!! SOMETHING WENT WRONG")
        summ = 0.
        for jj in range(nreacs):
            summ += cond[jj]*stoic[jj]*self.flux[self.reacs_by_met[ii][jj]]
        return summ

    def get_status(self):
        if self.lp.status == "opt":
            return "Optimal"
        elif self.lp.status == "undef":
            return "Undefined"
        elif self.lp.status == "feas":
            return "Maybe not optimal"
        elif self.lp.status == "infeas" or self.lp.status == "nofeas":
            return "Unfeasible"
        elif self.lp.status == "unbnd":
            return "Unbound (check constraints)"
        else:
            return "Unknown Error"

    __repr__ = __str__