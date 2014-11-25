'''
Whole-cell knowledge base haystack indices

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17

Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
'''
from datetime import datetime

from haystack import indexes

import cyano.models as cmodels

class EntryIndex(indexes.SearchIndex):
    text = indexes.CharField(document=True, use_template=True)
    wid = indexes.CharField(model_attr="wid")
    model_type = indexes.CharField(model_attr='model_type__model_name', faceted=True)
    #author = indexes.CharField(model_attr='detail__user', faceted=True)

    def index_queryset(self, using=None):
        return self.get_model().objects.filter(detail__date__lte=datetime.now())

    class Meta:
        abstract = True
        pass

class SpeciesComponentIndex(EntryIndex, indexes.Indexable):
    species_id = indexes.MultiValueField()
    species_wid = indexes.MultiValueField()

    def prepare_species_id(self, obj):
        return obj.species.pk

    def prepare_species_wid(self, obj):
        return obj.species.wid

    def get_model(self):
        return cmodels.SpeciesComponent

    class Meta:
        pass
