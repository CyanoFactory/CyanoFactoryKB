from django.core.exceptions import ObjectDoesNotExist
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

        if 'request' in self.context:
            fields = self.context['request'].QUERY_PARAMS.get('fields')
            if fields:
                fields = fields.split(',')
                # Drop any fields that are not specified in the `fields` argument.
                allowed = set(fields)
                existing = set(self.fields.keys())
                for field_name in existing - allowed:
                    self.fields.pop(field_name)


class WidField(serializers.SlugRelatedField):
    class WidOnlyObject(object):
        def __init__(self, wid):
            self.pk = 0
            self.wid = wid

    def __init__(self, *args, **kwargs):
        kwargs['slug_field'] = 'wid'
        super(WidField, self).__init__(*args, **kwargs)

    def get_attribute(self, instance):

        obj = instance.__class__.objects.filter\
            (pk=instance.pk).values_list("%s__wid" % self.field_name, flat=True)

        if len(obj) == 0:
            return WidField.WidOnlyObject(None)

        return WidField.WidOnlyObject(obj[0])


class WidManyRelatedField(serializers.ManyRelatedField):
    def to_representation(self, iterable):
        return list(iterable.values_list("wid", flat=True))


class WidFieldMany(serializers.RelatedField):
    def to_internal_value(self, data):
        pass

    def to_representation(self, value):
        pass

    @classmethod
    def many_init(cls, *args, **kwargs):
        list_kwargs = {'child_relation': cls(*args, **kwargs)}
        kwargs.pop("queryset")
        for key in kwargs.keys():
            list_kwargs[key] = kwargs[key]
        return WidManyRelatedField(**list_kwargs)


class CyanoSerializer(DynamicFieldsMixin, serializers.ModelSerializer):
    pass


class Entry(CyanoSerializer):
    class Meta:
        model = cmodels.Entry
        fields = cmodels.Entry._meta.field_list

additional_fields = ['species', 'url']


class BasketComponent(CyanoSerializer):
    component = WidField(read_only=True)
    species = WidField(read_only=True)

    class Meta:
        model = cmodels.BasketComponent
        fields = ['component', 'species']


class Basket(CyanoSerializer):
    components = BasketComponent(many=True, read_only=True)

    class Meta:
        model = cmodels.Basket
        fields = ['id', 'name', 'components']


class SpeciesComponent(CyanoSerializer):
    species = WidField(read_only=True)
    type = WidFieldMany(allow_null=True, many=True, queryset=cmodels.Type.objects.all(), required=False)
    url = serializers.URLField(source='get_absolute_url')

    class Meta:
        model = cmodels.SpeciesComponent
        fields = additional_fields


class Genome(SpeciesComponent):
    gc_content = serializers.ReadOnlyField(source='get_gc_content')

    class Meta:
        model = cmodels.Genome
        fields = cmodels.Genome._meta.field_list + SpeciesComponent.Meta.fields + ['gc_content']


class Chromosome(SpeciesComponent):
    class Meta:
        model = cmodels.Chromosome
        fields = cmodels.Chromosome._meta.field_list + SpeciesComponent.Meta.fields


class Plasmid(SpeciesComponent):
    class Meta:
        model = cmodels.Plasmid
        fields = cmodels.Plasmid._meta.field_list + SpeciesComponent.Meta.fields


class ChromosomeFeature(SpeciesComponent):
    class Meta:
        model = cmodels.ChromosomeFeature
        fields = cmodels.ChromosomeFeature._meta.field_list + SpeciesComponent.Meta.fields


class Compartment(SpeciesComponent):
    class Meta:
        model = cmodels.Compartment
        fields = cmodels.Compartment._meta.field_list + SpeciesComponent.Meta.fields


class Gene(SpeciesComponent):
    chromosome = WidField(read_only=True)

    def get_chromosome(self, gene):
        return cmodels.Genome.objects.filter(pk=gene.chromosome_id).values_list('wid', flat=True)[0]

    class Meta:
        model = cmodels.Gene
        fields = cmodels.Gene._meta.field_list + SpeciesComponent.Meta.fields


class Metabolite(SpeciesComponent):
    class Meta:
        model = cmodels.Metabolite
        fields = cmodels.Metabolite._meta.field_list + SpeciesComponent.Meta.fields


class Parameter(SpeciesComponent):
    class Meta:
        model = cmodels.Parameter
        fields = cmodels.Parameter._meta.field_list + SpeciesComponent.Meta.fields


class Pathway(SpeciesComponent):
    class Meta:
        model = cmodels.Pathway
        fields = cmodels.Pathway._meta.field_list + SpeciesComponent.Meta.fields


class Process(SpeciesComponent):
    class Meta:
        model = cmodels.Process
        fields = cmodels.Process._meta.field_list + SpeciesComponent.Meta.fields


class ProteinComplex(SpeciesComponent):
    class Meta:
        model = cmodels.ProteinComplex
        fields = cmodels.ProteinComplex._meta.field_list + SpeciesComponent.Meta.fields


class ProteinMonomer(SpeciesComponent):
    class Meta:
        model = cmodels.ProteinMonomer
        fields = cmodels.ProteinMonomer._meta.field_list + SpeciesComponent.Meta.fields


class Reaction(SpeciesComponent):
    class Meta:
        model = cmodels.Reaction
        fields = cmodels.Reaction._meta.field_list + SpeciesComponent.Meta.fields


class State(SpeciesComponent):
    class Meta:
        model = cmodels.State
        fields = cmodels.State._meta.field_list + SpeciesComponent.Meta.fields


class Stimulus(SpeciesComponent):
    class Meta:
        model = cmodels.Stimulus
        fields = cmodels.Stimulus._meta.field_list + SpeciesComponent.Meta.fields


class TranscriptionUnit(SpeciesComponent):
    class Meta:
        model = cmodels.TranscriptionUnit
        fields = cmodels.TranscriptionUnit._meta.field_list + SpeciesComponent.Meta.fields


class TranscriptionalRegulation(SpeciesComponent):
    class Meta:
        model = cmodels.TranscriptionalRegulation
        fields = cmodels.TranscriptionalRegulation._meta.field_list + SpeciesComponent.Meta.fields


class Type(SpeciesComponent):
    class Meta:
        model = cmodels.Type
        fields = cmodels.Type._meta.field_list + SpeciesComponent.Meta.fields


class PublicationReference(SpeciesComponent):
    class Meta:
        model = cmodels.PublicationReference
        fields = cmodels.PublicationReference._meta.field_list + SpeciesComponent.Meta.fields


class MassSpectrometryJob(SpeciesComponent):
    class Meta:
        model = cmodels.MassSpectrometryJob
        fields = cmodels.MassSpectrometryJob._meta.field_list + SpeciesComponent.Meta.fields


class Peptide(SpeciesComponent):
    class Meta:
        model = cmodels.Peptide
        fields = cmodels.Peptide._meta.field_list + SpeciesComponent.Meta.fields


class MassSpectrometryProtein(SpeciesComponent):
    class Meta:
        model = cmodels.MassSpectrometryProtein
        fields = cmodels.MassSpectrometryProtein._meta.field_list + SpeciesComponent.Meta.fields
