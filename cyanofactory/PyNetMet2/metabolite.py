from .util import python_2_unicode_compatible, create_sid

import libsbml

@python_2_unicode_compatible
class Metabolite(object):
    '''
    This class defines a chemical reaction object. The input should be an
    string containing a reaction in OptGene format.
    '''
    def __init__(self, model, met=None):
        if met is None:
            self._sbml = libsbml.Species(3, 1)
            model.getListOfSpecies().append(self._sbml)
        else:
            self._sbml = met
            if not self.meta_id:
                self.meta_id = self.id

    def get_name_or_id(self):
        return self.name or self.meta_id

    @property
    def id(self):
        return self._sbml.getId()

    @id.setter
    def id(self, val):
        self._sbml.setId(create_sid(val))

    @property
    def meta_id(self):
        return self._sbml.getId()

    @meta_id.setter
    def meta_id(self, val):
        self._sbml.setMetaId(create_sid(val))

    @property
    def name(self):
        return self._sbml.getName()

    @name.setter
    def name(self, val):
        self._sbml.setName(val)

    @property
    def external(self):
       return self._sbml.getBoundaryCondition()

    @external.setter
    def external(self, val):
        self._sbml.setBoundaryCondition(val)

    @property
    def compartment(self):
        return self._sbml.getCompartment()

    @compartment.setter
    def compartment(self, val):
        self._sbml.setCompartment(val)
