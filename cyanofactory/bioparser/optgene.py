import re

class OptGeneParser:
	def __init__(self, filename):
		self._file = filename
		self.errors = []
		self.reactions = {}
		self.constraints = {}
		self.external_metabolites = {}
		self.parse()
		
	# via http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
	@staticmethod
	def _enum(*sequential, **named):
		enums = dict(zip(sequential, range(len(sequential))), **named)
		return type('Enum', (), enums)

	Sections = _enum("none",
					"unknown",
					"reactions",
					"constraints",
					"external metabolites",
					"unconstrained metabolites",
					"extracellular metabolites",
					"extracellular metabolites (unconstrained)")
					#"objective",
					#"maximize",
					#"minimize")

	def _add_error(self, line_no, msg):
		self.errors.append("L" + str(line_no) + ": " + msg)
		
	def _takewhile(self, predicate, iterable, inclusive = False):
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
		error = False
		multi_comment = False
		current_section = OptGeneParser.Sections.none
		with open(self._file) as f:
			line_no = 0
			for line in f:
				line_no += 1
				line = re.sub(r"(\s+)", r" ", line)
				line = line.strip()
				
				if "#" in line: # Remove comments
					line = line[:line.find("#")]
				
				if multi_comment: # In comment check for end %
					if not "%" in line:
						continue
				
				if "%" in line: # Multiline comment
					multi_comment = not multi_comment
					line = line[:line.find("%")]
				
				if len(line) == 0:
					continue;
				
				if line[0] == "-": # Section start?
					sec = line[1:].lower()
					if hasattr(OptGeneParser.Sections, sec):
						current_section = getattr(OptGeneParser.Sections, sec)
					else:
						self._add_error(line_no, "Skipping unknown section " + line[1:])
						current_section = OptGeneParser.Sections.unknown
					continue
				
				if current_section == OptGeneParser.Sections.none:
					self._add_error(line_no, "Invalid identifier outside of section: " + line)
					continue
				if current_section == OptGeneParser.Sections.unknown:
					continue
				
				elif current_section == OptGeneParser.Sections.reactions:
					line = re.split("\s+", line) # Tokenize
					
					# Entry name
					name, line = self._takewhile(lambda x: x[-1] != ":", line, inclusive = True)
					name = " ".join(name)[:-1].strip()
					if len(line) == 0:
						self._add_error(line_no, "Invalid entry in reactions")
						continue
					
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
							pass # identifier
						
						# Identifier...
						iname, line = self._takewhile(lambda x: x != "+" and x != "->" and x != "<->", line)
						
						if len(iname) == 0:
							self._add_error(line_no, "Invalid entry in reactions")
							error = True
							break
						
						# Must be +, <- or ->
						if state == 0:
							target = reactants
						else:
							target = products
						target += [[value, " ".join(iname)]]
						
						if len(line) == 0:
							break
						
						if line[0] == "->" or line[0] == "<->":
							if arrow != "":
								self._add_error(line_no, "More then one reaction arrow: " + line[0])
								error = True
								break
							arrow = line[0]		
							state = 1
						
						line = line[1:]

					if not error:
						self.reactions[name] = {\
							"name": name,
							"index": len(self.reactions),
							"reactants": reactants,
							"arrow": arrow,
							"products": products}
					error = False
				
				elif current_section == OptGeneParser.Sections.constraints:
					# rightmost [ defines interval
					pos = line.rfind("[")
					if pos == -1:
						self._add_error(line_no, "Invalid entry in constraints: " + line[0])
						continue
					
					name = line[:pos].strip()
					m = re.match(r"\[\s?(\S+)\s?,\s?(\S+)\s?\]$", line[pos:])
					if m:
						try:
							begin = float(m.group(1))
							end = float(m.group(2))
						
							if begin > end:
								self._add_error(line_no, "Invalid interval (begin > end): " + str(begin) + " " + str(end))
								continue
						except ValueError:
							self._add_error(line_no, "Interval values must be numbers: " + m.group(1) + " " + m.group(2))
							continue
								
						self.constraints[name] = {\
							"name": name,
							"index": len(self.constraints),
							"interval": [[begin, end]]}
						
						if not name in self.reactions:
							self._add_error(line_no, "Constraint for unknown Reaction: " + name)
					else:
						self._add_error(line_no, "No interval found")
				
				elif current_section == getattr(OptGeneParser.Sections, "external metabolites") or\
						current_section == getattr(OptGeneParser.Sections, "unconstrained metabolites") or\
						current_section == getattr(OptGeneParser.Sections, "extracellular metabolites") or\
						current_section == getattr(OptGeneParser.Sections, "extracellular metabolites (unconstrained)"):
					# trivial
					self.external_metabolites[line] = {\
							"name": line,
							"index": len(self.external_metabolites)}
					
					if line in self.reactions:
						self._add_error(line_no, "External Metabolite conflicts with Reaction name: " + line)