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

from __future__ import print_function

from codecs import open
import os
from .network import *
from .enzyme import *
from xml.dom.minidom import parseString
import re
import six

class Metabolism(object):
    """
        This class defines a metabolism object (set of chemical reactions,
        metabolites and a network).
    """

    def __init__(self, filein, filetype="auto", reactions=list(), constraints=list(),
                                     external=list(), objective=list(), design_objective=list(), fromfile=True):
        if fromfile:
            if filetype == "auto":
                filetype = self.autodetect_format(filein)

            if filetype == "opt":
                [reactions, constraints, external, objective, design_objective] = self.prepare_opt(filein)
            elif filetype == "sbml":
                [reactions, constraints, external, objective] = self.prepare_sbml(filein)

        self.file_name = filein
        self.constraints = constraints
        self.external = external
        self.objective = objective
        self.design_objective = design_objective

        # Get reactions and pathways
        pathnames = []
        enzymes = []
        pathna = ""
        for line in reactions:
            if "#" in line:
                parts = line.split("#")
                comments = parts[1].strip()
                if parts[0] == "":
                    pathnames.append(comments)
                    pathna = comments
                elif ":" in parts[0] and "->" in parts[0] and line[0] != "#":
                    enzi = Enzyme(parts[0],pathna)
                    enzymes.append(enzi)
            elif ":" in line and "->" in line:
                enzi = Enzyme(line,pathna)
                enzymes.append(enzi)

        self.pathnames = pathnames
        self.enzymes = enzymes

        for line in constraints:
            if not line.startswith("#"):
                line = line.strip()
                # rightmost [ defines interval
                pos = line.rfind("[")
                name = line[:pos].strip()
                m = re.match(r"\[\s?(\S+)\s?,\s?(\S+)\s?\]$", line[pos:])
                if m:
                    try:
                        begin = float(m.group(1))
                        end = float(m.group(2))
                    except ValueError:
                        continue

                    for enz in self.enzymes:
                        if enz.name == name:
                            enz.constraint = (begin, end)

        # Objective
        objs = []
        for ele in self.objective:
            if ele[0] != "#":
                if "#" in ele:
                    bla = ele.split("#")
                    bla = bla[0]
                else:
                    bla = ele
                line = bla.rsplit(" ", 2)
                if len(line) == 3:
                    objs.append([line[0],line[2]])
                if len(line) == 2:
                    objs.append([line[0],line[1]])
                if len(line) == 1:
                    objs.append([line[0],"1"])

        # Design objective
        design_objs = []
        for ele in self.design_objective:
            if ele[0] != "#":
                if "#" in ele:
                    bla = ele.split("#")
                    bla = bla[0]
                else:
                    bla = ele
                line = bla.rsplit(" ", 2)
                if len(line) == 3:
                    design_objs.append([line[0],line[2]])
                if len(line) == 2:
                    design_objs.append([line[0],line[1]])
                if len(line) == 1:
                    design_objs.append([line[0],"1"])

        self.design_obj = design_objs
        self.obj = objs

        # Most attributes calculated in the calcs routine.
        self.calcs()

    def __repr__(self):
        """
        Show number of reactions and metabolites in the metabolism.
        """
        aa = "# Reactions:" + str(len(self.enzymes)) + "\n"
        aa += "# Metabolites:" + str(len(self.mets)) + "\n"
        return aa

    def calcs(self):
        """
        Creates the metabolism attributes based on the enzymes list.
        """
        # Reads metabolites
        subs = {}
        metabol = []
        reac_per_met = []
        ii = 0
        for reac in self.enzymes:
            for sus in reac.substrates:
                if sus not in metabol:
                    metabol.append(sus)
                    subs[sus] = ii
                    ii += 1
                    reac_per_met.append([reac.name])
                else:
                    isus = metabol.index(sus)
                    reac_per_met[isus].append(reac.name)
            for prod in reac.products:
                if prod not in metabol:
                    metabol.append(prod)
                    subs[prod] = ii
                    ii += 1
                    reac_per_met.append([reac.name])
                else:
                    iprod = metabol.index(prod)
                    reac_per_met[iprod].append(reac.name)
        reacs_per_met = [len(ele) for ele in reac_per_met]
        # Prepares external mets and M matrix
        list_rev = []
        list_irr = []
        external_in = []
        external_out = []
        transport = []

        # Throw away transport reactions, they are regenerated
        self.enzymes = [item for item in self.enzymes if item.pathway != "_TRANSPORT_"]

        enzinv = {}
        for ii, enzy in enumerate(self.enzymes):
            enzinv[enzy.name] = ii

            if enzy.reversible:
                list_rev.append(ii)
            else:
                list_irr.append(ii)
        self.external = list(set(self.external)|set(external_in))
        self.external = list(set(self.external)|set(external_out))

        # Apends reactions for external mets
        enzisinv = ii+1
        self.pathnames = []
        self.pathnames.append("_TRANSPORT_")
        for line in self.external:
            if line[0] != "#":
                if (line not in external_in) and (line not in external_out) and (line in metabol):
                    self.enzymes.append(Enzyme(line + "_transp : <-> " +
                                               line, "_TRANSPORT_"))
                    enzinv[line + "_transp"]=enzisinv
                    list_rev.append(enzisinv)
                    transport.append(enzisinv)
                    external_in.append(line)
                    external_out.append(line)
                    enzisinv += 1

        # writes reactions in pathways
        pathways = [[] for ele in self.pathnames]
        for ii in range(len(self.enzymes)):
            pathway = self.enzymes[ii].pathway
            if pathway in self.pathnames:
                ipath = self.pathnames.index(pathway)
            else:
                self.pathnames.append(pathway)
                ipath = self.pathnames.index(pathway)
                pathways.append([])
            pathways[ipath].append(ii)

        # Constructs the M matrix (metabolites connections)
        nsubs = len(subs)
        M = [[0 for iii in range(nsubs)] for jjj in range(nsubs)]
        for enzy in self.enzymes:
            for su in enzy.substrates:
                for pr in enzy.products:
                    M[subs[su]][subs[pr]] = 1
                    if enzy.reversible:
                        M[subs[pr]][subs[su]] = 1

        self.dic_enzs = enzinv
        self.pathways = pathways
        self.transport = transport
        self.dic_mets = subs
        self.mets = metabol
        self.external_in = external_in
        self.external_out = external_out
        self.M = M
        self.reac_rev = list_rev
        self.reac_irr = list_irr    
        self.reac_per_met = reac_per_met
        self.reacs_per_met = reacs_per_met
        self._net = None

        # Check for duplicated enzyme names and rename
        non_uniq_run = False
        non_uniq = filter(lambda x: self.has_reaction(x), self.mets)
        for x in non_uniq:
            self.rename_reaction(x, x + "_r", no_update=True)
            non_uniq_run = True
        if non_uniq_run:
            self.calcs()

    @property
    def net(self):
        if self._net is None:
            self.net = Network(self.M, self.mets)

        return self._net

    @net.setter
    def net(self, value):
        self._net = value

    def add_reacs(self, reacs):
        """
        Adds new reactions to the metabolism.
        """
        # Reads metabolites
        for reac in reacs:
            r = Enzyme(reac,"GP")
            if r.name in self.dic_enzs:
                raise ValueError("Reaction with same name already in model: " + r.name)

        subs = self.dic_mets
        metabol = self.mets
        reac_per_met = self.reac_per_met
        ii = len(self.mets)
        nreacs_old = len(self.enzymes)
        for jj, rea in enumerate(reacs):
            reac = Enzyme(rea,"GP")
            for sus in reac.substrates:
                if sus not in metabol:
                    metabol.append(sus)
                    subs[sus] = ii
                    ii += 1
                    reac_per_met.append([reac.name])
                else:
                    isus = metabol.index(sus)
                    reac_per_met[isus].append(reac.name)
            for prod in reac.products:
                if prod not in metabol:
                    metabol.append(prod)
                    subs[prod] = ii
                    ii += 1
                    reac_per_met.append([reac.name])
                else:
                    iprod = metabol.index(prod)
                    reac_per_met[iprod].append(reac.name)
            self.enzymes.append(reac)
            if reac.reversible:
                self.reac_rev.append(jj+nreacs_old)
            else:
                self.reac_irr.append(jj+nreacs_old)
        # Prepares external mets and M matrix
        for ii, enzy in enumerate(self.enzymes[nreacs_old:]):
            self.dic_enzs[enzy.name] = ii+nreacs_old
            if enzy.issues:
                if len(enzy.products) == 0:
                    self.external_out.extend(enzy.metabolites)
                    self.transport.append(ii+nreacs_old)
                elif len(enzy.substrates) == 0:
                    self.external_in.extend(enzy.metabolites)
                    self.transport.append(ii+nreacs_old)
        # writes reactions in pathways
        for ii in range(len(self.enzymes[nreacs_old:])):
            pathway = self.enzymes[ii+nreacs_old].pathway
            if pathway in self.pathnames:
                ipath = self.pathnames.index(pathway)
            else:
                self.pathnames.append(pathway)
                ipath = self.pathnames.index(pathway)
                self.pathways.append([])
            self.pathways[ipath].append(ii+nreacs_old)
        self.dic_mets = subs
        self.mets = metabol
        M = self.M_matrix()
        self.M = M
        self.calcs()

    def pop(self, iname):
        """
        Remove reaction number iname from the metabolism.
        """
        self.enzymes.pop(iname)
        self.calcs()

    def add_reaction(self, reaction):
        """
        Add reaction string to model
        """
        self.add_reacs([reaction])
        self.calcs()

    def get_reaction(self, name):
        """
        Get reaction with name
        """
        idx = self.dic_enzs.get(name)
        if idx is None:
            return None

        return self.enzymes[idx]

    def get_reaction_with_pathway(self, pathway_name):
        return list(filter(lambda x: x.pathway == pathway_name, self.enzymes))

    def rename_reaction(self, reac_from, reac_to, no_update=False):
        reac = self.get_reaction(reac_from)

        if reac is None:
            raise ValueError("Reaction not in model: " + reac_from)

        if self.has_reaction(reac_to):
            raise ValueError("Reaction already in model: " + reac_to)

        for obj in self.obj:
            if obj[0] == reac_from:
                obj[0] = reac_to

        for obj in self.design_obj:
            if obj[0] == reac_from:
                obj[0] = reac_to

        reac.name = reac_to

        if not no_update:
            self.calcs()

    def remove_reaction(self, name):
        """
        Remove reaction with passed name from model
        """
        idx = self.dic_enzs.get(name)
        if idx is None:
            raise ValueError("Reaction not in model: " + name)

        self.obj = list(filter(lambda x: x[0] != name, self.obj))
        self.design_obj = list(filter(lambda x: x[0] != name, self.obj))

        self.pop(idx)

    def has_reaction(self, name):
        return self.get_reaction(name) is not None

    def rename_metabolite(self, met_from, met_to):
        if re.match(r"^[0-9]", met_to):
            raise ValueError("Name must not begin with a number")

        if re.match(r"(^<?\-> | <?\-> | <?\->$)", met_to):
            raise ValueError("Name must not contain lonely <-> or ->")

        if re.match(r"(^\+ | \+ )", met_to):
            raise ValueError("Lonely + only allowed at end of name")

        if not self.has_metabolite(met_from):
            raise ValueError("Metabolite not in model: " + met_from)

        if self.has_metabolite(met_to):
            raise ValueError("Metabolite already in model: " + met_to)

        for reaction in self.enzymes:
            for i, metabolite in enumerate(reaction.substrates):
                if metabolite == met_from:
                    reaction.substrates[i] = met_to

            for i, metabolite in enumerate(reaction.products):
                if metabolite == met_from:
                    reaction.products[i] = met_to

            for i, line in enumerate(self.external):
                if line[0] != "#":
                    if line == met_from:
                        self.external[i] = met_to

        self.calcs()

    def has_metabolite(self, name):
        return name in self.dic_mets

    def make_metabolite_external(self, name):
        if not self.has_metabolite(name):
            raise ValueError("Metabolite not in model: " + name)

        for i, line in enumerate(self.external):
            if line[0] != "#":
                if line == name:
                    # Already external
                    return

        self.external.append(name)

        self.calcs()

    def make_metabolite_internal(self, name):
        if not self.has_metabolite(name):
            raise ValueError("Metabolite not in model: " + name)

        self.external = list(filter(lambda x: x[0] == "#" or x != name, self.external))

        self.calcs()

    def autodetect_format(self, infile):
        if isinstance(infile, six.string_types):
            with open(infile, encoding='utf-8') as f:
                if f.readline().startswith("<?xml"):
                    return "sbml"
                return "opt"
        else:
            try:
                if infile.readline().startswith("<?xml"):
                    return "sbml"
                return "opt"
            finally:
                infile.seek(0)

    def prepare_opt(self, infile):
        if isinstance(infile, six.string_types):
            with open(infile, encoding='utf-8') as f:
                return self.parse_opt(f)
        else:
            return self.parse_opt(infile)

    def parse_opt(self, infile):
        """
        Reads the OptGene file and gathers the reaction information,
        contrains, external metabolites and objective.     
        """
        # reads the file
        reactions = []
        constraints = []
        external = []
        objective = []
        design_objective = []

        # Find key positions in the file
        keyy = ""
        on = True

        for data in infile:
            line = data.strip()

            # Skip empty lines
            if len(line) == 0:
                continue

            if line[0] == "-":
                keyy = line
            elif line[0] == "%" and on:
                on = False
            elif line[0] == "%" and not on:
                on = True
            elif keyy == "-REACTIONS" and line != keyy and on:
                reactions.append(line)
            elif keyy == "-CONSTRAINTS" and line != keyy and on:
                constraints.append(line)
            elif keyy == "-EXTERNAL METABOLITES" and line != keyy and on:
                external.append(line)
            elif keyy == "-OBJ" and line != keyy and on:
                objective.append(line)
            elif keyy == "-DESIGNOBJ" and line != keyy and on:
                design_objective.append(line)
            else:
                pass

        # Return results
        return [reactions, constraints, external, objective, design_objective]

    def prepare_sbml(self, infile):
        if isinstance(infile, six.string_types):
            with open(infile, encoding='utf-8') as f:
                return self.parse_sbml(f)
        else:
            return self.parse_sbml(infile)

    def parse_sbml(self, filein):
        """
            Reads information from sbml file.
        """
        # Reads sbml
        data = filein.read()
        # Parse the file
        dom = parseString(data)
        reactions = dom.getElementsByTagName('reaction')
        species = dom.getElementsByTagName('species')
        nreacs = len(reactions)
        # Dictionary with species
        specie = {}
        externals = []
        for i in range(len(species)):
           name = species[i]._attrs["name"].value.replace(":", ";")
           try:
               iddd = species[i]._attrs["id"].value
           except:
               iddd = name
           try:
               bcon = species[i]._attrs["boundaryCondition"].value
           except:
               bcon = "false"
           try:
              compart = species[i]._attrs["compartment"].value
           except:
              compart = ""
           specie[iddd] = (name, compart)
           if bcon == "true":
               externals.append(iddd)
        # Prepares reactions in optgene format
        reacts = ["# GPW"]
        constrs = []
        objs = []
        for i in range(nreacs):
           reacs = reactions[i].getElementsByTagName('listOfReactants')
           prods = reactions[i].getElementsByTagName('listOfProducts')
           params = reactions[i].getElementsByTagName('parameter')
           try:
              reacname = reactions[i]._attrs["id"].value
           except:
              reacname = reactions[i]._attrs["name"].value
           try:
              rever = reactions[i]._attrs["reversible"].value
           except:
              rever = "true"
           rever = rever.replace(" ","")
           sustrate = ""
           for j in range(len(reacs)):
              spe = reacs[j].getElementsByTagName('speciesReference')
              for k in range(len(spe)):
                 kk = spe[k]._attrs["species"].value
                 try:
                    stoic = spe[k]._attrs["stoichiometry"].value
                 except:
                    stoic = " 1 "
                 sustrate += " + "+str(float(stoic))+" "+kk.replace(" ", "")
           product = ""
           for j in range(len(prods)):
              spe=prods[j].getElementsByTagName('speciesReference')
              for k in range(len(spe)):
                 kk = spe[k]._attrs["species"].value
                 try:
                    stoic = spe[k]._attrs["stoichiometry"].value
                 except:
                    stoic = " 1 "
                 product += " + "+str(float(stoic))+" "+kk.replace(" ", "")
           stri = ""
           sustrate = sustrate[3:len(sustrate)]
           product = product[3:len(product)]
           if rever == "false":
              reacao = sustrate + " -> " + product
           else:
              reacao = sustrate + " <-> " + product
           reacao = reacname.replace(" ", "") + " : " + reacao
           reacts.append(reacao)
           try:
              lbound = [ele._attrs["value"].value for ele in params \
                                    if ele._attrs["id"].value == "LOWER_BOUND"]
              ubound = [ele._attrs["value"].value for ele in params \
                                    if ele._attrs["id"].value == "UPPER_BOUND"]
              objcoe = [ele._attrs["value"].value for ele in params \
                          if ele._attrs["id"].value == "OBJECTIVE_COEFFICIENT"]
           except:
              lbound = [ele._attrs["value"].value for ele in params \
                                  if ele._attrs["name"].value == "LOWER_BOUND"]
              ubound = [ele._attrs["value"].value for ele in params \
                                  if ele._attrs["name"].value == "UPPER_BOUND"]
              objcoe = [ele._attrs["value"].value for ele in params \
                         if ele._attrs["name"].value == "OBJECTIVE_COEFFICIENT"]
           try:
               lbound = float(lbound[0])
           except:
               lbound = None
           try:
               ubound = float(ubound[0])
           except:
               ubound = None
           try:
               objcoe = float(objcoe[0])
           except:
               objcoe = None
           if lbound != None and ubound != None:
               constr = reacname+" ["+str(lbound)+", "+str(ubound)+"]" 
               constrs.append(constr)
           if objcoe != None and objcoe != 0.:
               obj = reacname+" 1 "+str(objcoe)
               objs.append(obj)
        self.dic_specie = specie
        return [reacts, constrs, externals, objs]

    def dump(self, fileout="model.txt", filetype="opt", printgenes = False, \
             dic_genes = {}, alllimits=False, allobjs=False):
        """
        Creates an OptGene or SBML file with the metabolism.
        (Reactions, Constraints, External Metabolites and Objective).
        """
        if isinstance(fileout, six.string_types):
            fil = open(fileout, "w", encoding='utf-8')
        else:
            fil = fileout

        if filetype == "opt":
            print("-REACTIONS", file=fil)
            for ii, pathname in enumerate(self.pathnames):
                if pathname == "_TRANSPORT_":
                    continue
                print("", file=fil)
                if pathname is not None:
                    print("# " + pathname, file=fil)
                for reac in self.pathways[ii]:
                    fil.write(self.enzymes[reac])
                    fil.write("\n")
            print("", file=fil)
            print("-CONSTRAINTS", file=fil)
            print("", file=fil)
            for constr in self.enzymes:
                if constr.constraint[0] is None or constr.constraint[1] is None:
                    continue
                print(constr.name + " [" + str(constr.constraint[0]) + ", " + str(constr.constraint[1]) + "]", file=fil)
            print("", file=fil)
            print("-EXTERNAL METABOLITES", file=fil)
            print("", file=fil)
            exts = set(self.external) | set(self.external_in)
            exts = exts|set(self.external_out)
            exts = list(exts)
            exts.sort()
            for ext in exts:
                print(ext, file=fil)

            if self.obj:
                print("", file=fil)
                print("-OBJ", file=fil)
                print("", file=fil)
                for obj in self.obj:
                    print(obj[0] + " 1 " + str(obj[1]), file=fil)

            if self.design_obj:
                print("", file=fil)
                print("-DESIGNOBJ", file=fil)
                print("", file=fil)
                for obj in self.design_obj:
                    print(obj[0] + " 1 " + str(obj[1]), file=fil)

        elif filetype == "sbml":
            if printgenes:
                enzs = dic_genes.keys()
            print('<?xml version="1.0" encoding="utf-8"?>', file=fil)
            file_name_only = os.path.splitext(os.path.basename(self.file_name))[0]
            print('<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" version="1" level="3" xmlns:xhtml="http://www.w3.org/1999/xhtml">', file=fil)
            print('<model id="%s" name="%s" >' % (file_name_only, file_name_only), file=fil)
            tab = "    "
            print("<notes>", file=fil)
            print(tab+"<xhtml:body>", file=fil)
            print(tab*2+"<xhtml:p>Generated by PyNetMet tools. INTERTECH@UPV</xhtml:p>", file=fil)
            print(tab+"</xhtml:body>", file=fil)
            print("</notes>", file=fil)
            print("<listOfCompartments>", file=fil)
            print(tab+'<compartment id="Extra_cellular" constant="false"/>', file=fil)
            print(tab+'<compartment id="Cytosol" constant="false"/>', file=fil)
            print("</listOfCompartments>", file=fil)
            print("<listOfSpecies>", file=fil)
            exts = set(self.external) | set(self.external_in)
            exts = exts|set(self.external_out)
            exts = list(exts)
            for met in self.mets:
                if met in exts:
                    bcon = "true"
                    comp = "Extra_cellular"
                else:
                    bcon = "false"
                    comp = "Cytosol"
                print(tab+'<species id="%s" name="%s" boundaryCondition="%s" compartment="%s" constant="false" hasOnlySubstanceUnits="false"/>'%(met,met,bcon,comp), file=fil)
            print("</listOfSpecies>", file=fil)
            print("<listOfReactions>", file=fil)
            objs = [ele[0] for ele in self.obj]
            for ii,enzy in enumerate(self.enzymes):
                if enzy.pathway == "_TRANSPORT_":
                    continue
                if enzy.reversible:
                    rev = "true"
                else:
                    rev = "false"
                print(tab+'<reaction id="%s" name="%s" reversible="%s" fast="false">'%(enzy.name,enzy.name,rev), file=fil)
                if printgenes:
                    print(2*tab+'<notes>', file=fil)
                    print(3*tab+'</xhtml:body>', file=fil)
                    if enzy.name in enzs:
                        assoc = "Genes : "+", ".join(dic_genes[enzy.name])
                    else:
                        assoc = "No known genes."
                    print(4*tab+'<xhtml:p>'+assoc+'</xhtml:p>', file=fil)
                    print(3*tab+'</xhtml:body>', file=fil)
                    print(2*tab+'</notes>', file=fil)
                print(2*tab+"<listOfReactants>", file=fil)
                for met in enzy.substrates:
                    print(3*tab+'<speciesReference species="%s" stoichiometry="%s" constant="false"/>'%(met,str(enzy.stoic_n(met))), file=fil)
                print(2*tab+"</listOfReactants>", file=fil)
                print(2*tab+"<listOfProducts>", file=fil)
                for met in enzy.products:
                    print(3*tab+'<speciesReference species="%s" stoichiometry="%s" constant="false"/>'%(met,str(enzy.stoic_n(met))), file=fil)
                print(2*tab+"</listOfProducts>", file=fil)

                [lb, ub] = enzy.constraint
                obj = False
                if enzy.name in objs:
                    obj = True
                    iobj = objs.index(enzy.name)
                if (ub != None and lb != None) or obj or alllimits or allobjs:
                    print(2*tab+"<kineticLaw>", file=fil)
                    print(3*tab+"<listOfParameters>", file=fil)
                    if lb != None and ub != None:
                        print(4*tab+'<parameter id="LOWER_BOUND" value="%s" constant="false"/>'%str(lb), file=fil)
                        print(4*tab+'<parameter id="UPPER_BOUND" value="%s" constant="false"/>'%str(ub), file=fil)
                    elif alllimits:
                        if lb == None:
                            lb = -99999999.
                        if ub == None:
                            ub = 99999999.
                        print(4*tab+'<parameter id="LOWER_BOUND" value="%s" constant="false"/>'%str(lb), file=fil)
                        print(4*tab+'<parameter id="UPPER_BOUND" value="%s" constant="false"/>'%str(ub), file=fil)
                    if obj:
                        print(4*tab+'<parameter id= "OBJECTIVE_COEFFICIENT" value="%s" constant="false"/>'%self.obj[iobj][1], file=fil)
                    elif allobjs:
                        print(4*tab+'<parameter id= "OBJECTIVE_COEFFICIENT" value="0.0" constant="false"/>', file=fil)
                    print(3*tab+"</listOfParameters>", file=fil)
                    print(2*tab+"</kineticLaw>", file=fil)
                print(tab+'</reaction>', file=fil)
            print("</listOfReactions>", file=fil)
            print("</model>", file=fil)
            print("</sbml>", file=fil)

        if isinstance(fileout, six.string_types):
            fil.close()

    def M_matrix(self,symetric=False):
        """
        Returns the adjacent matrix for the metabolites. It can be directed or
        indirected if symmetric is True or False.
        """
        nsubs = len(self.dic_mets)
        subs=self.dic_mets
        M = [[0 for iii in range(nsubs)] for jjj in range(nsubs)]
        for enzy in self.enzymes:
            for su in enzy.substrates:
                for pr in enzy.products:
                    M[subs[su]][subs[pr]] = 1
                    if enzy.reversible or symetric:
                        M[subs[pr]][subs[su]] = 1
        return M

    def fba(self):
        from .fba import FBA
        return FBA(self)

    def M_matrix_reacs(self,symmetric=False):
        """
        Returns the bipartite adjacent matrix for the reactions and 
        metabolites. It can be directed or indirected if symmetric 
        is True or False.
        """
        nsubs = len(self.dic_mets)
        nreacs = len(self.enzymes)
        subs = self.dic_mets
        enzs = self.dic_enzs
        reacs = self.enzymes
        nnodes = nsubs + nreacs
        nodes = self.mets + [ele.name for ele in reacs]
        M = [[0 for iii in range(nnodes)] for jjj in range(nnodes)]
        for ii, reac in enumerate(reacs):
            for su in reac.substrates:
                M[subs[su]][enzs[reac.name]+nsubs] = 1
                if reac.reversible or symmetric:
                    M[enzs[reac.name]+nsubs][subs[su]] = 1
            for pr in reac.products:
                M[enzs[reac.name]+nsubs][subs[pr]] = 1
                if reac.reversible or symmetric:
                    M[subs[pr]][enzs[reac.name]+nsubs] = 1
        return [M,nodes]

    def bad_reacs(self):
        """
        Filters reactions that are for sure disconnected from the giant
        component or that can not participate in the FBA.
        """
        # Selects biggest connected component
        self.net.components()
        self.net.disc_comps.sort(key=lambda x:len(x), reverse=True)
        topop = []
        mets_in_reacs = [enzy.metabolites for enzy in self.enzymes]
        for comp in self.net.disc_comps[1:]:
            for met in comp:
                bla = [self.mets[met] in reac for reac in mets_in_reacs]
                reacs = [ii for ii,ele in enumerate(bla) if ele]
                topop = list(set(topop)|set(reacs))
        # Searches for reactions where one substrate and one product appear
        # only in this reaction.
        for ii,enzy in enumerate(self.enzymes):
            bla_sus = [self.reacs_per_met[self.dic_mets[ele]] == 1
                     for ele in enzy.substrates]
            bla_pro = [self.reacs_per_met[self.dic_mets[ele]] == 1
                     for ele in enzy.products]
            if len(bla_sus) and len(bla_pro):
                sus = reduce(lambda x, y:x or y, bla_sus)
                pro = reduce(lambda x, y:x or y, bla_pro)
            else:
                sus = False
                pro = False
            if (sus and pro) and ii not in topop:
                topop.append(ii)
        for ele in self.transport:
            if "_transp" in self.enzymes[ele].name:
                topop.append(ele)
        topop=list(set(topop))
        self.pathnames.pop(self.pathnames.index("_TRANSPORT_"))
        topop.sort(reverse=True)
        poped = [self.enzymes.pop(ele) for ele in topop]
        self.poped = poped
        self.calcs()

    def write_log(self, fileout="Log_model.txt", nshow=50):
        """
        Writes a log file with the metabolism information.
        """
        ## Prepares LOG file
        probs = open(fileout, "w", encoding='utf-8')
        print("-Reactions with issues:", file=probs)
        print("----------------", file=probs)
        print(" ", file=probs)
        print(" ", file=probs)
        print(" ", file=probs)
        print(" General:", file=probs)
        print(" ", file=probs)
        for enzy in self.enzymes:
            if enzy.issues:
                print("____________________________________", file=probs)
                print("Issue in reaction " + enzy.name, file=probs)
                print(enzy.issues_info, file=probs)
                print("____________________________________", file=probs)
        # Searches for reactions where one substrate and one product appear
        #   only in this reaction
        print(" ", file=probs)
        print(" ", file=probs)
        for ii, enzy in enumerate(self.enzymes):
            bla_sus = [self.reacs_per_met[self.dic_mets[ele]] == 1
                       for ele in enzy.substrates]
            bla_pro = [self.reacs_per_met[self.dic_mets[ele]] == 1
                       for ele in enzy.products]
            if len(bla_sus) and len(bla_pro):
                sus = reduce(lambda x, y:x or y, bla_sus)
                pro = reduce(lambda x, y:x or y, bla_pro)
            else:
                sus = False
                pro = False
            if (sus and pro):
                print(" Possible disconnected reaction: " + enzy.name, file=probs)
        # General information about network
        print("      ", file=probs)
        print("      ", file=probs)
        print("      ", file=probs)
        print(" Information About Network     ", file=probs)
        print(" ----------- ----- -------     ", file=probs)
        print("      ", file=probs)
        print("Number of reactions: " + str(len(self.enzymes)), file=probs)
        print("Number of pathways: " + str(len(self.pathnames)), file=probs)
        print("Number of metabolites: " + str(len(self.mets)), file=probs)
        print("Reversible: " + self.nreac_rev, file=probs)
        print("Irreversible:" + self.nreac_irr, file=probs)
        print("      ", file=probs)
        print("      ", file=probs)
        # Information about enzyme-metabolite connections
        print("      ", file=probs)
        print("   ENZYME-METABOLITE", file=probs)
        print("      ", file=probs)
        print(" Most connected Metabolites to enzymes: ", file=probs)
        print("      ", file=probs)
        con_met = sorted(range(len(self.mets)), key=lambda x:self.reacs_per_met[x],
                         reverse=True)
        for ele in con_met[0:nshow+1]:
            print("%30s : %5i" % (self.mets[ele], self.reacs_per_met[ele]), file=probs)

        print("      ", file=probs)
        print(" Most connected Metabolites as substrates: ", file=probs)
        print("      ", file=probs)
        substr = [sum([met in reac.substrates for reac in self.enzymes]+
                  [met in reac.products for reac in self.enzymes if
                  reac.reversible]) for met in self.mets]
        con_met = sorted(range(len(self.mets)), key=lambda x:substr[x],
                         reverse=True)
        for ele in con_met[0:nshow+1]:
            print("%30s : %5i" % (self.mets[ele], substr[ele]), file=probs)
        print("      ", file=probs)
        print(" Most connected Metabolites as products: ", file=probs)
        print("      ", file=probs)
        substr = [sum([met in reac.products for reac in self.enzymes]+
                  [met in reac.substrates for reac in self.enzymes if 
                  reac.reversible]) for met in self.mets]
        con_met = sorted(range(len(self.mets)), key=lambda x:substr[x],
                         reverse=True)
        for ele in con_met[0:nshow+1]:
            print("%30s : %5i" % (self.mets[ele], substr[ele]), file=probs)
        print("      ", file=probs)
        print("   METABOLITE-METABOLITE", file=probs)
        print("      ", file=probs)
        print(" Most connected Metabolites: ", file=probs)
        print("      ", file=probs)
        nlinks = [len(ele) for ele in self.net.neigbs]
        con_met = sorted(range(len(self.nmets)), key=lambda x:nlinks[x],
                         reverse=True)
        for ele in con_met[0:nshow+1]:
            print("%30s : %5i" % (self.mets[ele], nlinks[ele]), file=probs)
        print("      ", file=probs)
        print(" Most out-connected Metabolites: ", file=probs)
        print("      ", file=probs)
        nlinks = [len(ele) for ele in self.net.linksout]
        con_met = sorted(range(len(self.nmets)), key=lambda x:nlinks[x],
                         reverse=True)
        for ele in con_met[0:nshow+1]:
            print("%30s : %5i" % (self.mets[ele], nlinks[ele]), file=probs)
        print("      ", file=probs)
        print(" Most in-connected Metabolites: ", file=probs)
        print("      ", file=probs)
        nlinks = [len(ele) for ele in self.net.linksin]
        con_met = sorted(range(len(self.nmets)), key=lambda x:nlinks[x],
                         reverse=True)
        for ele in con_met[0:nshow+1]:
            print("%30s : %5i" % (self.mets[ele], nlinks[ele]), file=probs)
        self.log_file = "logs/Log_" + self.file_name + ".txt"
        print(" ", file=probs)
        probs.close()

