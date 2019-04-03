from .util import python_2_unicode_compatible, create_sid

import libsbml

@python_2_unicode_compatible
class Objective(object):
    '''
    This class defines a chemical reaction object. The input should be an
    string containing a reaction in OptGene format.
    '''
    def __init__(self, model, obj=None):
        if obj is None:
            self._sbml = libsbml.Objective(3, 1)
            model.getPlugin("fbc").getListOfObjectives().append(self._sbml)
        else:
            self._sbml = obj

        if len(self._sbml.getListOfFluxObjectives()) == 0:
            self._sbml.createFluxObjective()

    @property
    def id(self):
        return self._sbml.getId()

    @id.setter
    def id(self, val):
        self._sbml.setId(create_sid(val))

    @property
    def type(self):
        return self._sbml.getType()

    @type.setter
    def type(self, val):
        self._sbml.setType(val)

    @property
    def coefficient(self):
        return self._sbml.getFluxObjective(0).getCoefficient()

    @coefficient.setter
    def coefficient(self, val):
        self._sbml.getFluxObjective(0).setCoefficient(val)

    @property
    def reaction(self):
        return self._sbml.getFluxObjective(0).getReaction()

    @reaction.setter
    def reaction(self, val):
        self._sbml.getFluxObjective(0).setReaction(val)

    def __getitem__(self, key):
        if key != 0 and key != 1:
            raise ValueError("bad index")

        return [self.reaction, '1' if self.type == "maximize" else '-1'][key]

    def __setitem__(self, key, value):
        if key != 0:
            raise ValueError("bad index")

        self.reaction = value[0]
        self.type = "maximize" if int(value[1]) > 0 else "minimize"
