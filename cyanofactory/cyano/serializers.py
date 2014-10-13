from rest_framework import serializers
import cyano.models as cmodels


class Entry(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Entry
        fields = cmodels.Entry._meta.field_list


class SpeciesComponent(serializers.ModelSerializer):
    class Meta:
        model = cmodels.SpeciesComponent
        fields = cmodels.SpeciesComponent._meta.field_list


class Genome(serializers.ModelSerializer):
    gc_content = serializers.Field(source='get_gc_content')

    class Meta:
        model = cmodels.Genome
        fields = cmodels.Genome._meta.field_list + ['gc_content']


class Chromosome(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Chromosome
        fields = cmodels.Chromosome._meta.field_list


class Plasmid(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Plasmid
        fields = cmodels.Plasmid._meta.field_list


class ChromosomeFeature(serializers.ModelSerializer):
    class Meta:
        model = cmodels.ChromosomeFeature
        fields = cmodels.ChromosomeFeature._meta.field_list


class Compartment(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Compartment
        fields = cmodels.Compartment._meta.field_list


class Gene(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Gene
        fields = cmodels.Gene._meta.field_list


class Metabolite(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Metabolite
        fields = cmodels.Metabolite._meta.field_list


class Parameter(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Parameter
        fields = cmodels.Parameter._meta.field_list


class Pathway(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Pathway
        fields = cmodels.Pathway._meta.field_list


class ProteinComplex(serializers.ModelSerializer):
    class Meta:
        model = cmodels.ProteinComplex
        fields = cmodels.ProteinComplex._meta.field_list


class ProteinMonomer(serializers.ModelSerializer):
    class Meta:
        model = cmodels.ProteinMonomer
        fields = cmodels.ProteinMonomer._meta.field_list


class Reaction(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Reaction
        fields = cmodels.Reaction._meta.field_list


class State(serializers.ModelSerializer):
    class Meta:
        model = cmodels.State
        fields = cmodels.State._meta.field_list


class Stimulus(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Stimulus
        fields = cmodels.Stimulus._meta.field_list


class TranscriptionUnit(serializers.ModelSerializer):
    class Meta:
        model = cmodels.TranscriptionUnit
        fields = cmodels.TranscriptionUnit._meta.field_list


class TranscriptionalRegulation(serializers.ModelSerializer):
    class Meta:
        model = cmodels.TranscriptionalRegulation
        fields = cmodels.TranscriptionalRegulation._meta.field_list


class Type(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Type
        fields = cmodels.Type._meta.field_list


class PublicationReference(serializers.ModelSerializer):
    class Meta:
        model = cmodels.PublicationReference
        fields = cmodels.PublicationReference._meta.field_list


class MassSpectrometryJob(serializers.ModelSerializer):
    class Meta:
        model = cmodels.MassSpectrometryJob
        fields = cmodels.MassSpectrometryJob._meta.field_list


class Peptide(serializers.ModelSerializer):
    class Meta:
        model = cmodels.Peptide
        fields = cmodels.Peptide._meta.field_list


class MassSpectrometryProtein(serializers.ModelSerializer):
    class Meta:
        model = cmodels.MassSpectrometryProtein
        fields = cmodels.MassSpectrometryProtein._meta.field_list
