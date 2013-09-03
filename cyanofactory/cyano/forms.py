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
		label = 'a file to import'
		)
	reason = forms.CharField(
		required = True,
		widget = forms.TextInput,
		label = "Enter a short summary")
	
	# GenBank specific
	chromosome = forms.CharField(
		required = False,
		widget = forms.TextInput,
		label = "Enter name of chromosome/plasmid")

	chromosome_wid = forms.SlugField(
		required = False,
		widget = forms.TextInput,
		label = "And a new identifier")

	def __init__(self, *args, **kwargs):
		super(ImportDataForm, self).__init__(*args, **kwargs)
		
		choices = []
		for species in cmodels.Species.objects.values('wid', 'name').all():
			choices.append((species['wid'], species['name'], ))
		self.fields['species'].choices = choices
		#self.fields['data_type'].choices = ImportDataForm.data_type.choices
	
	def clean(self):
		cleaned_data = super(ImportDataForm, self).clean()

		data_type = cleaned_data.get("data_type")
		
		if data_type:
			if data_type == "genbank":
				chromosome = cleaned_data.get("chromosome")
				chromosome_wid = cleaned_data.get("chromosome_wid")
				
				if not chromosome:
					self._errors["chromosome"] = self.error_class(["Chromosome required for GenBank"])
					del cleaned_data["chromosome"]
				
				if not chromosome_wid:
					self._errors["chromosome_wid"] = self.error_class(["Identifier required for GenBank"])
					del cleaned_data["chromosome_wid"]
		
		return cleaned_data

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
		label = "And a new unique identifier"
	)
	reason = forms.CharField(
		required = True,
		widget = forms.TextInput,
		label = "Enter a short summary"
	)
		
	def __init__(self, *args, **kwargs):
		super(ImportSpeciesForm, self).__init__(*args, **kwargs)
	
	def clean(self):
		cleaned_data = super(ImportSpeciesForm, self).clean()
		
		new_wid = cleaned_data.get("new_wid")
		
		if new_wid:
			if cmodels.Species.objects.filter(wid = new_wid).exists():
				self._errors["new_wid"] = self.error_class(['The identifier specified is already in use'])
				del cleaned_data["new_wid"]

		return cleaned_data