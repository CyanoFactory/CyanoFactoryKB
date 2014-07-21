"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from collections import OrderedDict

import re

class OptGeneParser:
    """Parses the specified OptGene file.
    
    OptGene (or BioOpt) is a text data format to describe interaction
    between different metabolites.
    
    The original specification (http://129.16.106.142/tools.php?c=bioopt)
    is incompatible to some sample files, because of that this parser is more
    tolerant about the file structure. The format basicly consists of
    identifiers delimited by operators. Identifiers are allowed to contain
    spaces (spaces to the left and right are removed)
    
    An OptGene file consists of multiple sections:
    
    - ``-REACTIONS``
    - ``-CONSTRAINTS``
    - ``-EXTERNAL METABOLITES`` (or ``-UNCONSTRAINED METABOLITES`` and ``-EXTRACELLULAR METABOLITES (UNCONSTRAINED)``)
    - ``-OBJECTIVE``
    - ``-MAXIMIZE``
    - ``-MINIMIZE``
    
    ``-REACTIONS`` contains a list of chemical reactions, e.g.
    
      GLK1_1 : alpha-D-Glucose + 2.31 ATP -> alpha-D-Glucose 6-phosphate + 3/4 ADP
    
    The identifier is everything to the left until the parser encounters a colon.
    
    The identifier is everything to the left until the parser encounters a colon
    followed by a whitespace, no whitespace before the colon is needed. Then the
    formula parsing follows. In a formula the identifiers are delimited by ``+``
    and the reactants and products are delimited by an arrow, ``<->`` means
    reversible, ``->`` means not reversible. Only one arrow is allowed in the
    formula. Before the identifier can be a    float or a fraction number of two
    floats, when this is missing 1 is assumed.
    
    ``-CONSTRAINTS`` contain values for minflux and maxflux.
    
      CO2xtI [min, max]

    Everything up to the last ``[`` is parsed as an identifier (must be a known
    identifier from a ``REACTIONS`` item, otherwise an error is reported. ``min``
    and ``max`` are two float values. If a reaction is not listed the values
    ``[-1000, +1000]`` for reversible and ``[0, +max]`` for irreversible are
    assumed.
    
    ``-EXTERNAL METABOLITES`` contains external metabolites, is a list with
    one identifier per line, the whole line is used as the identifier name.
    An error is reported when a metabolite has the same name as a reaction.
    
    Other sections are not supported yet and skipped. A note is added to the
    error list for every unknown section.
    
    ``Variables``:
    
    - ``errors``: Contains parsing errors
    - ``reactions``: All reactions
    - ``constraints``: All constraints
    - ``external_metabolites``: All external metabolites
    """
    def __init__(self, stream):
        """Parses the specified OptGene stream"""
        self._file = stream
        self.header = []
        self.errors = []
        self.line_no = 0
        self.current_section = self._parse_header
        self.metabolites = OrderedDict()
        self.reactions = OrderedDict()
        self.temp_comment_storage = []
        self.multi_line_comment = False

        self.parse()

    def _parse_header(self, line):
        """
        Reads "file header"
        Everything until first section, is usually empty
        """
        for item in self.temp_comment_storage:
            self.header.append(item)
        self.temp_comment_storage = []
        self.header.append(line)

    class Metabolite(object):
        def __init__(self, name, model=None):
            self._name = name
            self.external = False
            self.comments = []
            self.model = model

        def __new__(cls, name, model=None):
            if model and name in model.metabolites:
                raise TypeError("Metabolite already in model")

            instance = super(OptGeneParser.Metabolite, cls).__new__(cls)
            model.metabolites[name] = instance

            return instance

        def __repr__(self):
            return self.name

        @property
        def name(self):
            return self._name

        @name.setter
        def name(self, value):
            if value in self.model.metabolites and value is not self.name:
                raise ValueError("Duplicate metabolite name")

            self.model.metabolites = OrderedDict([(value, v) if k == self.name else (k, v) for k, v in self.model.metabolites.items()])

            self._name = value

        #@property
        #def model(self):
        #    return self._model

        #@model.setter
        #def model(self, value):
        #    # self assign
        #    if self.model is value:
        #        return
        #
        #    self.model.metabolites

    class Reaction(object):
        def __init__(self, name, model=None, reactants=[], reactants_stoic=[], products=[], products_stoic=[], reversible=False):
            self._name = name
            self.reactants = reactants
            self.reactants_stoic = reactants_stoic
            self.products = products
            self.products_stoic = products_stoic
            self.reversible = reversible
            self.comments = {}
            self._constraint = None
            self.flux = None
            self.model = model

            optgene.reactions[name] = self

        def __new__(cls, optgene, name, *args, **kwargs):
            if name in optgene.reactions:
                raise TypeError("Reaction already in model")

            instance = super(OptGeneParser.Reaction, cls).__new__(cls)
            optgene.reactions[name] = instance

            return instance

        @staticmethod
        def from_string(self, string):
            line = re.split("\s+", string)  # Tokenize

            # Entry name
            name, line = self._takewhile(lambda x: not x.endswith(":"), line, inclusive=True)
            name = " ".join(name)[:-1].strip()

            if len(line) == 0:
                raise ValueError("Invalid reaction: " + string)

            reactants = []
            arrow = ""
            products = []
            state = 0

            while len(line) > 0:
                value = 1
                try:
                    if "/" in line[0]:
                        fl = line[0].split("/")
                        value = float(fl[0]) / float(fl[1])
                    else:
                        value = float(line[0])
                    line = line[1:]
                except ValueError:
                    pass  # identifier

                # Identifier...
                iname, line = self._takewhile(lambda x: x != "+" and x != "->" and x != "<->", line)
                # Test if there a multiple operators, all except the last belong to the name...
                ops, line = self._takewhile(lambda x: x == "+" or x == "->" or x == "<->", line)
                if len(ops) > 0:
                    iname += ops[:-1]
                    line.insert(0, ops[-1])

                iname = " ".join(iname)
                #print "XXX", iname, line

                if len(iname) == 0:
                    raise ValueError("Invalid reaction: " + string)

                # Must be +, <- or ->
                if state == 0:
                    target = reactants
                else:
                    target = products

                if not iname in self.metabolites:
                    OptGeneParser.Metabolite(self, iname)
                metabolite = self.metabolites[iname]

                # reference to reactants or products
                target += [[value, metabolite]]

                if len(line) == 0:
                    break

                if line[0] == "->" or line[0] == "<->":
                    if arrow != "":
                        raise ValueError("More then one reaction arrow: " + string)

                    arrow = line[0]
                    state = 1

                line = line[1:]

            return OptGeneParser.Reaction(
                name,
                reactants=map(lambda x: x[1], reactants),
                reactants_stoic=map(lambda x: x[0], reactants),
                products=map(lambda x: x[1], products),
                products_stoic=map(lambda x: x[0], products),
                reversible=arrow == "<->"
            )

        def is_constrained(self):
            return self._constraint is not None

        @property
        def name(self):
            return self._name

        @name.setter
        def name(self, value):
            if value in self.model.reactions and value is not self.name:
                raise ValueError("Duplicate reaction name")

            self.model.reactions = OrderedDict([(value, v) if k == self.name else (k, v) for k, v in self.model.reactions.items()])

            self._name = value

        @property
        def constraint(self):
            if self._constraint is None:
                if self.reversible:
                    return [None, None]
                else:
                    return [0, None]
            return self._constraint

        @constraint.setter
        def constraint(self, value):
            self._constraint = value

        def __repr__(self):
            return "{} : {} {} {}".format(
                self.name,
                " + ".join(map(lambda x: "%s %s" % (x[0], x[1]), zip(self.reactants_stoic, self.reactants))),
                "<->" if self.reversible else "->",
                " + ".join(map(lambda x: "%s %s" % (x[0], x[1]), zip(self.products_stoic, self.products))))

    def has_metabolite(self, name):
        return name in self.metabolites

    def has_reaction(self, name):
        return name in self.reactions

    def _parse_reactions(self, line):
        reaction = OptGeneParser.Reaction.from_string(line)

        if name in self.reactions:
            raise ValueError("Reaction " + name + " already in model")

        reaction = OptGeneParser.Reaction(
            self,
            name,
            reactants=map(lambda x: x[1], reactants),
            reactants_stoic=map(lambda x: x[0], reactants),
            products=map(lambda x: x[1], products),
            products_stoic=map(lambda x: x[0], products),
            reversible=arrow == "<->"
        )
        reaction.comments["reactions"] = self.temp_comment_storage
        self.temp_comment_storage = []

    def _parse_constraints(self, line):
        # rightmost [ defines interval
        pos = line.rfind("[")
        if pos == -1:
            raise ValueError("Invalid entry in constraints: " + line)

        name = line[:pos].strip()
        m = re.match(r"\[\s?(\S+)\s?,\s?(\S+)\s?\]$", line[pos:])
        if m:
            try:
                begin = float(m.group(1))
                end = float(m.group(2))

                if begin > end:
                    raise ValueError("Invalid interval (begin > end): " + str(begin) + " " + str(end))
            except ValueError:
                raise ValueError("Interval values must be numbers: " + m.group(1) + " " + m.group(2))

            if not self.has_reaction(name):
                raise ValueError("Constraint for unknown Reaction: " + name)

            reaction = self.reactions[name]
            reaction.comments["constraints"] = self.temp_comment_storage
            reaction.constraint = [begin, end]
        else:
            raise ValueError("No interval found")

    def _parse_external_metabolites(self, line):
        # Singleton
        if not self.has_metabolite(line):
            raise ValueError("Unused external metabolite: " + line)

        metabolite = self.metabolites[line]
        metabolite.comments = self.temp_comment_storage
        metabolite.external = True

        # Create reaction for external meta
        OptGeneParser.Reaction(
            self,
            metabolite.name + "_transp",
            products=[metabolite],
            products_stoic=[1.0],
            reversible=True
        )

        self.temp_comment_storage = []

    def _parse_objective(self, line):
        m = re.match(r"^(.+)\s+(-?[0-9]+)\s+(-?[0-9]+)$", line)
        if m:
            pass
        else:
            raise ValueError("Invalid Objective " + line)

    def _parse_design_objective(self, line):
        m = re.match(r"^(.+)\s+(-?[0-9]+)\s+(-?[0-9]+)$", line)
        if m:
            pass
        else:
            raise ValueError("Invalid Design Objective " + line)

    def _parse_unknown(self, line):
        print "Unknown", line

    Sections = {
        "reactions": _parse_reactions,
        "constraints": _parse_constraints,
        "external metabolites": _parse_external_metabolites,
        "unconstrained metabolites": _parse_external_metabolites,
        "extracellular metabolites": _parse_external_metabolites,
        "extracellular metabolites (unconstrained)": _parse_external_metabolites,
        "objective": _parse_objective,
        "designobj": _parse_design_objective
    }

    def _add_error(self, line_no, msg):
        self.errors.append("L" + str(line_no) + ": " + msg)

    def _takewhile(self, predicate, iterable, inclusive=False):
        res = []
        for x in iterable:
            if predicate(x):
                res.append(x)
            else:
                if inclusive:
                    res.append(x)
                break

        return res, iterable[len(res):]

    def parse(self):
        for line_no, line in enumerate(self._file, start=1):
            line = line.strip()
            if len(line) == 0:
                continue

            if line.startswith("%") or self.multi_line_comment:
                if line.startswith("%"):
                    self.multi_line_comment = not self.multi_line_comment
                self.temp_comment_storage.append(line)
                continue
            elif line.startswith("#"):
                self.temp_comment_storage.append(line)
                continue
            elif line.startswith("-"):
                if len(self.temp_comment_storage) > 0:
                    self._add_error(line_no, "Lonely comment omitted: " + "; ".join(self.temp_comment_storage))
                    self.temp_comment_storage = []
                self.current_section = self.Sections.get(line[1:].lower(), OptGeneParser._parse_unknown)
                continue

            try:
                self.current_section(self, line)
            except ValueError as e:
                self._add_error(line_no, str(e))
                self.temp_comment_storage = []

    def write(self, handle):
        for item in self.header:
            handle.write(item + "\n")

        handle.write("\n-REACTIONS\n\n")
        for item in self.reactions.keys():
            reaction = self.reactions[item]

            for comment in reaction.comments.get("reactions", []):
                handle.write(comment + "\n")

            handle.write(str(reaction) + "\n")

        handle.write("\n-CONSTRAINTS\n\n")
        for item in self.reactions.keys():
            reaction = self.reactions[item]

            if not reaction.is_constrained():
                continue

            for comment in reaction.comments.get("constraints", []):
                handle.write(comment + "\n")

            handle.write(reaction.name + " " + str(reaction.constraint) + "\n")

        handle.write("\n-EXTERNAL METABOLITES\n\n")
        for item in self.metabolites.keys():
            metabolite = self.metabolites[item]

            if metabolite.external:
                for comment in metabolite.comments:
                    handle.write(comment + "\n")

                handle.write(metabolite.name + "\n")

    def fba(self):
        """
        Simplex algorithm through glpk.
        """

        # Calculate stoic
        stoic = []
        for i, reac in enumerate(self.reactions):
            reac = self.reactions[reac]
            for j, substrate in enumerate(zip(reac.reactants, reac.reactants_stoic)):
                subst, subst_stoic = substrate
                mi = self.metabolites.values().index(subst)
                stoic.append((mi, i, -subst_stoic))

            for j, product in enumerate(zip(reac.products, reac.products_stoic)):
                prod, prod_stoic = product
                mi = self.metabolites.values().index(prod)
                stoic.append((mi, i, prod_stoic))

        import glpk

        lp = glpk.LPX()
        lp.name = " FBA SOLUTION "
        lp.rows.add(len(self.metabolites))
        lp.cols.add(len(self.reactions))

        for i, met in enumerate(self.metabolites):
            met = self.metabolites[met]

            lp.rows[i].name = met.name

        for i, reac in enumerate(self.reactions):
            reac = self.reactions[reac]

            lp.cols[i].name = reac.name

        #constraints

        for i, reac in enumerate(self.reactions):
            reac = self.reactions[reac]

            lp.cols[i].bounds = tuple(reac.constraint)

        for n in xrange(len(self.metabolites)):
            lp.rows[n].bounds = (0., 0.)

        ###### Objective
        lista = [0. for ele in self.reactions]
        lista[4] = 1.0

        lp.obj[:] = lista[:]
        lp.obj.maximize = True  # self.max
        ###### Matrix
        lp.matrix = stoic   #self.Mstoic[:]
        lp.simplex()

        for flux, reac in zip(lp.cols, self.reactions):
            flux = flux.value
            reac = self.reactions[reac]

            reac.flux = flux
