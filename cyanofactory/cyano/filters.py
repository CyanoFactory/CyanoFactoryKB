import django_filters

import cyano.models as cmodels

class CyanoFilter(django_filters.FilterSet):
    pass


class Entry(CyanoFilter):
    class Meta:
        model = cmodels.Entry
        fields = cmodels.Entry._meta.field_list


class SpeciesComponent(CyanoFilter):
    species = django_filters.CharFilter(name="species__wid")
    parent = django_filters.CharFilter(name="parent__wid")
    type = django_filters.CharFilter(name="type__wid")
    publication_references = django_filters.CharFilter(name="publication_references__wid")
    cross_references = django_filters.CharFilter(name="cross_references__wid")

    class Meta:
        model = cmodels.SpeciesComponent
        fields = ['species', 'parent']


class Genome(SpeciesComponent):
    class Meta:
        model = cmodels.Genome
        fields = cmodels.Genome._meta.field_list + SpeciesComponent.Meta.fields


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
    chromosome = django_filters.CharFilter(name="chromosome__wid")
    amino_acid = django_filters.CharFilter(name="amino_acid__wid")

    class Meta:
        model = cmodels.Gene
        fields = cmodels.Gene._meta.field_list + SpeciesComponent.Meta.fields


class Metabolite(SpeciesComponent):
    class Meta:
        model = cmodels.Metabolite
        fields = cmodels.Metabolite._meta.field_list + SpeciesComponent.Meta.fields


class Parameter(SpeciesComponent):
    state = django_filters.CharFilter(name="state__wid")
    process = django_filters.CharFilter(name="process__wid")
    reaction = django_filters.CharFilter(name="reactions__wid")
    molecules = django_filters.CharFilter(name="molecules__wid")

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


class Protein(SpeciesComponent):
    class Meta:
        model = cmodels.Protein
        fields = cmodels.Protein._meta.field_list + SpeciesComponent.Meta.fields


class ProteinComplex(SpeciesComponent):
    formation_process = django_filters.CharFilter(name="formation_process__wid")

    class Meta:
        model = cmodels.ProteinComplex
        fields = cmodels.ProteinComplex._meta.field_list + SpeciesComponent.Meta.fields


class ProteinMonomer(SpeciesComponent):
    gene = django_filters.CharFilter(name="gene__wid")
    signal_sequence = django_filters.CharFilter(name="signal_sequence__wid")

    class Meta:
        model = cmodels.ProteinMonomer
        fields = cmodels.ProteinMonomer._meta.field_list + SpeciesComponent.Meta.fields


class Reaction(SpeciesComponent):
    processes = django_filters.CharFilter(name="processes__wid")
    states = django_filters.CharFilter(name="states__wid")

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
    genes = django_filters.CharFilter(name="genes__wid")

    class Meta:
        model = cmodels.TranscriptionUnit
        fields = cmodels.TranscriptionUnit._meta.field_list + SpeciesComponent.Meta.fields


class TranscriptionalRegulation(SpeciesComponent):
    transcription_unit = django_filters.CharFilter(name="transcription_unit__wid")
    transcription_factor = django_filters.CharFilter(name="transcription_factor__wid")

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
