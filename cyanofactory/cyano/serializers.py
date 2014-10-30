from rest_framework import serializers
import cyano.models as cmodels

# via https://gist.github.com/dbrgn/4e6fc1fe5922598592d6
class DynamicFieldsMixin(object):
    """
    A serializer mixin that takes an additional `fields` argument that controls
    which fields should be displayed.

    Usage::

        class MySerializer(DynamicFieldsMixin, serializers.HyperlinkedModelSerializer):
            class Meta:
                model = MyModel

    """
    def __init__(self, *args, **kwargs):
        super(DynamicFieldsMixin, self).__init__(*args, **kwargs)
        fields = self.context['request'].QUERY_PARAMS.get('fields')
        if fields:
            fields = fields.split(',')
            # Drop any fields that are not specified in the `fields` argument.
            allowed = set(fields)
            existing = set(self.fields.keys())
            for field_name in existing - allowed:
                self.fields.pop(field_name)


class CyanoSerializer(DynamicFieldsMixin, serializers.ModelSerializer):
    pass


class Entry(CyanoSerializer):
    class Meta:
        model = cmodels.Entry
        fields = cmodels.Entry._meta.field_list


class SpeciesComponent(CyanoSerializer):
    class Meta:
        model = cmodels.SpeciesComponent
        fields = cmodels.SpeciesComponent._meta.field_list


class Genome(CyanoSerializer):
    gc_content = serializers.Field(source='get_gc_content')

    class Meta:
        model = cmodels.Genome
        fields = cmodels.Genome._meta.field_list + ['gc_content']


class Chromosome(CyanoSerializer):
    class Meta:
        model = cmodels.Chromosome
        fields = cmodels.Chromosome._meta.field_list


class Plasmid(CyanoSerializer):
    class Meta:
        model = cmodels.Plasmid
        fields = cmodels.Plasmid._meta.field_list


class ChromosomeFeature(CyanoSerializer):
    class Meta:
        model = cmodels.ChromosomeFeature
        fields = cmodels.ChromosomeFeature._meta.field_list


class Compartment(CyanoSerializer):
    class Meta:
        model = cmodels.Compartment
        fields = cmodels.Compartment._meta.field_list


class Gene(CyanoSerializer):
    class Meta:
        model = cmodels.Gene
        fields = cmodels.Gene._meta.field_list


class Metabolite(CyanoSerializer):
    class Meta:
        model = cmodels.Metabolite
        fields = cmodels.Metabolite._meta.field_list


class Parameter(CyanoSerializer):
    class Meta:
        model = cmodels.Parameter
        fields = cmodels.Parameter._meta.field_list


class Pathway(CyanoSerializer):
    class Meta:
        model = cmodels.Pathway
        fields = cmodels.Pathway._meta.field_list


class ProteinComplex(CyanoSerializer):
    class Meta:
        model = cmodels.ProteinComplex
        fields = cmodels.ProteinComplex._meta.field_list


class ProteinMonomer(CyanoSerializer):
    class Meta:
        model = cmodels.ProteinMonomer
        fields = cmodels.ProteinMonomer._meta.field_list


class Reaction(CyanoSerializer):
    class Meta:
        model = cmodels.Reaction
        fields = cmodels.Reaction._meta.field_list


class State(CyanoSerializer):
    class Meta:
        model = cmodels.State
        fields = cmodels.State._meta.field_list


class Stimulus(CyanoSerializer):
    class Meta:
        model = cmodels.Stimulus
        fields = cmodels.Stimulus._meta.field_list


class TranscriptionUnit(CyanoSerializer):
    class Meta:
        model = cmodels.TranscriptionUnit
        fields = cmodels.TranscriptionUnit._meta.field_list


class TranscriptionalRegulation(CyanoSerializer):
    class Meta:
        model = cmodels.TranscriptionalRegulation
        fields = cmodels.TranscriptionalRegulation._meta.field_list


class Type(CyanoSerializer):
    class Meta:
        model = cmodels.Type
        fields = cmodels.Type._meta.field_list


class PublicationReference(CyanoSerializer):
    class Meta:
        model = cmodels.PublicationReference
        fields = cmodels.PublicationReference._meta.field_list


class MassSpectrometryJob(CyanoSerializer):
    class Meta:
        model = cmodels.MassSpectrometryJob
        fields = cmodels.MassSpectrometryJob._meta.field_list


class Peptide(CyanoSerializer):
    class Meta:
        model = cmodels.Peptide
        fields = cmodels.Peptide._meta.field_list


class MassSpectrometryProtein(CyanoSerializer):
    class Meta:
        model = cmodels.MassSpectrometryProtein
        fields = cmodels.MassSpectrometryProtein._meta.field_list
