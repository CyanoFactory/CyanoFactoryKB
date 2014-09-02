"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from collections import OrderedDict

import re

class OptGeneParser:
    """
    Parses the specified OptGene file.
    
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
        self.current_section_name = ""
        self.current_section = OptGeneParser._parse_header
        self.reactions = []
        self.temp_comment_storage = []
        self.multi_line_comment = False
        self.unknown_sections = OrderedDict()
        self._objective = None
        self._design_objective = None

        self.parse()

    @property
    def objective(self):
        return self._objective

    @objective.setter
    def objective(self, value):
        for reac in self.reactions:
            if reac is value:
                self._objective = value

        raise ValueError("Objective not in model: " + value.name)

    @property
    def design_objective(self):
        return self._objective

    @design_objective.setter
    def design_objective(self, value):
        for reac in self.reactions:
            if reac is value:
                self._design_objective = value

        raise ValueError("Design objective not in model: " + value.name)

    def add_reaction(self, reaction):
        for reac in self.reactions:
            if reaction.name == reac.name:
                raise ValueError("Reaction already in model: " + reac.name)

        self.reactions.append(reaction)

    def remove_reaction(self, reaction):
        for i, reac in enumerate(self.reactions):
            if reaction is reac:
                del self.reactions[i]

        raise ValueError("Reaction not in model: " + reac.name)

    def get_reaction(self, name):
        for reac in self.reactions:
            if reac.name == name:
                return reac

        return None

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
        __metabolites = []

        @staticmethod
        def instance(name):
            for met in OptGeneParser.Metabolite.__metabolites:
                if met.name == name:
                    return met

            met = OptGeneParser.Metabolite()
            met.name = name
            OptGeneParser.Metabolite.__metabolites.append(met)
            return met

        def __init__(self):
            self.name = ""
            self.external = False
            self.comments = []

        def __repr__(self):
            return self.name

    class Reaction(object):
        def __init__(self, name, reactants=[], reactants_stoic=[], products=[], products_stoic=[], reversible=False, disabled=False):
            self.name = name
            self.reactants = reactants
            self.reactants_stoic = reactants_stoic
            self.products = products
            self.products_stoic = products_stoic
            self.reversible = reversible
            self.comments = {}
            self._constraint = None
            self.flux = None
            self.disabled = disabled

        def replace_with(self, with_reaction):
            self.name = with_reaction.name
            self.reactants = with_reaction.reactants
            self.reactants_stoic = with_reaction.reactants_stoic
            self.products = with_reaction.products
            self.products_stoic = with_reaction.products_stoic
            self.reversible = with_reaction.reversible
            self.constraint = with_reaction.constraint
            self.disabled = with_reaction.disabled

        @staticmethod
        def from_string(string):
            line = re.split("\s+", string)  # Tokenize

            # Entry name
            name, line = OptGeneParser._takewhile(lambda x: not x.endswith(":"), line, inclusive=True)
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
                iname, line = OptGeneParser._takewhile(lambda x: x != "+" and x != "->" and x != "<->", line)
                # Test if there a multiple operators, all except the last belong to the name...
                ops, line = OptGeneParser._takewhile(lambda x: x == "+" or x == "->" or x == "<->", line)
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

                metabolite = OptGeneParser.Metabolite.instance(iname)

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

            disabled = False
            if name[:1] == "!":
                name = name[1:]
                disabled = True

            return OptGeneParser.Reaction(
                name,
                reactants=map(lambda x: x[1], reactants),
                reactants_stoic=map(lambda x: x[0], reactants),
                products=map(lambda x: x[1], products),
                products_stoic=map(lambda x: x[0], products),
                reversible=arrow == "<->",
                disabled=disabled
            )

        def is_constrained(self):
            return self._constraint is not None

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
            def reac_printer(arg):
                stoic, metabolite = arg

                if stoic == 1:
                    return str(metabolite)
                else:
                    if stoic == int(stoic):
                        stoic = int(stoic)

                    return "%s %s" % (stoic, str(metabolite))

            return "{}{} : {} {} {}".format(
                "!" if self.disabled else "",
                self.name,
                " + ".join(map(reac_printer, zip(self.reactants_stoic, self.reactants))),
                "<->" if self.reversible else "->",
                " + ".join(map(reac_printer, zip(self.products_stoic, self.products))))

    def has_metabolite(self, name):
        from itertools import chain

        for reac in self.reactions:
            for met in chain(reac.reactants, reac.products):
                if met.name == name:
                    return True
        return False

    def has_reaction(self, name):
        return self.get_reaction(name) is not None

    def get_metabolites(self):
        from itertools import chain

        metabolites = {}
        met_list = []

        for reac in self.reactions:
            for met in chain(reac.reactants, reac.products):
                if not met.name in metabolites:
                    metabolites[met.name] = 1
                    met_list.append(met)
        return met_list

    def _parse_reactions(self, line):
        reaction = OptGeneParser.Reaction.from_string(line)
        reaction.comments["reactions"] = self.temp_comment_storage

        self.add_reaction(reaction)

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

            reaction = self.get_reaction(name)
            reaction.comments["constraints"] = self.temp_comment_storage
            self.temp_comment_storage = []
            reaction.constraint = [begin, end]
        else:
            raise ValueError("No interval found")

    def _parse_external_metabolites(self, line):
        # Singleton
        if not self.has_metabolite(line):
            raise ValueError("Unused external metabolite: " + line)

        metabolite = OptGeneParser.Metabolite.instance(line)
        metabolite.external = True

        self.temp_comment_storage = []

    def _parse_objective(self, line):
        m = re.match(r"^(.+)\s+(-?[0-9]+)\s+(-?[0-9]+)$", line)
        if m:
            if self.objective is not None:
                raise ValueError("Only one objective supported: " + line)

            self.objective = self.get_reaction(m.group(1))
            self.objective.comments["objective"] = self.temp_comment_storage
            self.temp_comment_storage = []
        else:
            raise ValueError("Invalid Objective " + line)

    def _parse_design_objective(self, line):
        m = re.match(r"^(.+)\s+(-?[0-9]+)\s+(-?[0-9]+)$", line)
        if m:
            if self.design_objective is not None:
                raise ValueError("Only one design objective supported: " + line)

            self.design_objective = self.get_reaction(m.group(1))
            self.design_objective.comments["design objective"] = self.temp_comment_storage
            self.temp_comment_storage = []
        else:
            raise ValueError("Invalid Design Objective " + line)

    def _parse_unknown(self, line):
        if not self.current_section_name in self.unknown_sections:
            self.unknown_sections[self.current_section_name] = []

        for item in self.temp_comment_storage:
            self.unknown_sections[self.current_section_name].append(item)
        self.temp_comment_storage = []

        self.unknown_sections[self.current_section_name].append(line)

    Sections = {
        "reactions": _parse_reactions,
        "constraints": _parse_constraints,
        "external metabolites": _parse_external_metabolites,
        "unconstrained metabolites": _parse_external_metabolites,
        "extracellular metabolites": _parse_external_metabolites,
        "extracellular metabolites (unconstrained)": _parse_external_metabolites,
        "obj": _parse_objective,
        "objective": _parse_objective,
        "designobj": _parse_design_objective
    }

    def _add_error(self, line_no, msg):
        self.errors.append("L" + str(line_no) + ": " + msg)

    @staticmethod
    def _takewhile(predicate, iterable, inclusive=False):
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
                    self._add_error(line_no, "End of section comment omitted: " + "; ".join(self.temp_comment_storage))
                    self.temp_comment_storage = []
                self.current_section_name = line[1:]
                self.current_section = self.Sections.get(self.current_section_name.lower(), OptGeneParser._parse_unknown)
                continue

            try:
                self.current_section(self, line)
            except ValueError as e:
                self._add_error(line_no, str(e))
                self.temp_comment_storage = []

        # Check for reactions that have same name as metabolites
        metabolites = self.get_metabolites()
        for reaction in self.reactions:
            for metabolite in metabolites:
                if reaction.name == metabolite.name:
                    self._add_error(0, "Reaction %s conflicted with metabolite %s and was renamed" % (reaction.name, metabolite))
                    reaction.name += "_reac"

    def write(self, handle):
        for item in self.header:
            handle.write(item + "\n")

        handle.write("\n-REACTIONS\n\n")
        for reaction in self.reactions:
            for comment in reaction.comments.get("reactions", []):
                handle.write(comment + "\n")

            handle.write(str(reaction) + "\n")

        handle.write("\n-CONSTRAINTS\n\n")
        for reaction in self.reactions:

            if not reaction.is_constrained():
                continue

            for comment in reaction.comments.get("constraints", []):
                handle.write(comment + "\n")

            handle.write(reaction.name + " " + str(reaction.constraint) + "\n")

        handle.write("\n-EXTERNAL METABOLITES\n\n")
        for metabolite in filter(lambda m: m.external, self.get_metabolites()):
            for comment in metabolite.comments:
                handle.write(comment + "\n")

            handle.write(metabolite.name + "\n")

        if self.objective:
            handle.write("\n-OBJECTIVE\n\n")

            for comment in self.objective.comments.get("objective", []):
                handle.write(comment + "\n")

            handle.write(self.objective.name + " 1 1\n")

        if self.design_objective:
            handle.write("\n-DESIGNOBJ\n\n")

            for comment in self.design_objective.comments.get("design objective", []):
                handle.write(comment + "\n")

            handle.write(self.design_objective.name + " 1 1\n")

        for unknown_section in self.unknown_sections:
            handle.write("\n-" + unknown_section + "\n\n")

            for line in self.unknown_sections[unknown_section]:
                handle.write(line + "\n")

    def fba(self):
        """
        Simplex algorithm through glpk.
        """
        if self.objective is None:
            raise ValueError("No objective specified")

        # Create fake reactions for external metabolites
        fba_reactions = filter(lambda x: not x.disabled, self.reactions[:])
        for metabolite in filter(lambda x: x.external, self.get_metabolites()):
            reaction = OptGeneParser.Reaction(
                metabolite.name + "_transp",
                products=[metabolite],
                products_stoic=[1.0],
                reversible=True
            )

            fba_reactions.append(reaction)

        # Calculate stoic
        stoic = []
        metabolites = self.get_metabolites()
        for i, reac in enumerate(fba_reactions):
            for j, substrate in enumerate(zip(reac.reactants, reac.reactants_stoic)):
                subst, subst_stoic = substrate
                mi = metabolites.index(subst)
                stoic.append((mi, i, -subst_stoic))

            for j, product in enumerate(zip(reac.products, reac.products_stoic)):
                prod, prod_stoic = product
                mi = metabolites.index(prod)
                stoic.append((mi, i, prod_stoic))

        import glpk

        lp = glpk.LPX()
        lp.name = " FBA SOLUTION "
        lp.rows.add(len(metabolites))
        lp.cols.add(len(fba_reactions))

        for i, met in enumerate(metabolites):
            lp.rows[i].name = met.name

        for i, reac in enumerate(fba_reactions):
            lp.cols[i].name = reac.name

        #constraints

        for i, reac in enumerate(fba_reactions):
            lp.cols[i].bounds = tuple(reac.constraint)
            #print tuple(reac.constraint)

        for n in xrange(len(metabolites)):
            lp.rows[n].bounds = (0., 0.)

        ###### Objective
        lista = [0. for ele in fba_reactions]

        for i, reac in enumerate(fba_reactions):
            if self.objective is reac:
                lista[i] = 1.0

        lp.obj[:] = lista[:]
        lp.obj.maximize = True
        ###### Matrix
        lp.matrix = stoic
        lp.simplex()

        for flux, reac in zip(lp.cols, fba_reactions):
            flux = flux.value

            reac.flux = round(flux, 9)
