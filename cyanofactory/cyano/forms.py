'''
Whole-cell knowledge base forms

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17
'''

from django import forms
import cyano.helpers as chelpers
import cyano.models as cmodels

class SearchForm(forms.Form):
	q = forms.CharField(
		required = True,
		widget = forms.TextInput,
		label = 'query',
		help_text = 'Enter search term(s)',
		initial = ''
		)

class ExportDataForm(forms.Form):
	FORMAT_CHOICES = (
		('bib', 'BibTex'),
		('xlsx', 'Excel'),
		('html', 'HTML'),
		('json', 'JSON'),		
		('pdf', 'PDF'),
		('xml', 'XML'),
	)
	
	species = forms.ChoiceField(
		required = False,
		widget = forms.RadioSelect, 
		label = 'species',
		help_text = 'Select species to export'
		)
	model_type = forms.MultipleChoiceField(
		required = False,
		widget = forms.CheckboxSelectMultiple, 
		label = 'type',
		help_text = 'Select entry types to export'
		)
	all_model_types = forms.ChoiceField(
		choices = (
			('True', 'All'),
		),
		required = False,
		widget = forms.RadioSelect, 
		label = 'All',
		help_text = 'Select all entry types'
	)
	format = forms.ChoiceField(
		choices = FORMAT_CHOICES, 
		initial = FORMAT_CHOICES[0][0],
		required = True,
		widget = forms.RadioSelect, 		
		label = 'format',
		help_text = 'Select an output format'
		)
	
	def __init__(self, *args, **kwargs):
		super(ExportDataForm, self).__init__(*args, **kwargs)
		
		choices = []
		for species in cmodels.Species.objects.values('wid', 'name').all():
			choices.append((species['wid'], species['name'], ))
		self.fields['species'].choices = choices
		
		model_types = chelpers.getObjectTypes()
		models = chelpers.getModels()
		choices = []
		for model_type in model_types:
			choices.append((model_type, models[model_type]._meta.verbose_name_plural, ))
		self.fields['model_type'].choices = choices
		self.fields['model_type'].initial = model_types

class ImportDataForm(forms.Form):
	species = forms.ChoiceField(
		required = False,
		widget = forms.Select, 
		label = "species",
		help_text = "Select species to export"
		)
	data_type = forms.ChoiceField(
		required = True,
		choices = [
				#("fasta", "FASTA"),
				#("fastagene", "FASTA (Gene list)"),
				("genbank", "GenBank"),
				#("optgene", "OptGene / BioOpt"),
				("sbml", "System Biology Markup Language (SBML)"),
				("proopdb", "ProOpDB Operon Prediction")
				],
		widget = forms.RadioSelect,
		label = 'file format'
		)
	file = forms.FileField(
		required = True,
		widget = forms.ClearableFileInput, 
		label = 'file',
		help_text = 'Select a file (Excel or FASTA) to import'
		)
	reason = forms.CharField(
		required = True,
		widget = forms.TextInput,
		label = "reason",
		help_text = "Enter a summary for this import")
		
	def __init__(self, *args, **kwargs):
		super(ImportDataForm, self).__init__(*args, **kwargs)
		
		choices = []
		for species in cmodels.Species.objects.values('wid', 'name').all():
			choices.append((species['wid'], species['name'], ))
		self.fields['species'].choices = choices
		#self.fields['data_type'].choices = ImportDataForm.data_type.choices

class ImportSpeciesForm(forms.Form):
	new_species = forms.CharField(
		required = True,
		widget = forms.TextInput,
		label = "species",
		help_text = "Provide a new name"
		)
	new_wid = forms.SlugField(
		required = True,
		widget = forms.TextInput,
		label = "species",
		help_text = "Enter a new unique identifier"
	)
	reason = forms.CharField(
		required = True,
		widget = forms.TextInput,
		label = "reason",
		help_text = "Enter a summary for this import")
		
	def __init__(self, *args, **kwargs):
		super(ImportSpeciesForm, self).__init__(*args, **kwargs)
