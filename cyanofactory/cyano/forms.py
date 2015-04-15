'''
Whole-cell knowledge base forms

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17

Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
'''

from django import forms
import cyano.helpers as chelpers
from cyano.layout import Panel, Modal, ModalBody, ModalFooter, ModalHeader, CancelButton
import cyano.models as cmodels
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Fieldset, ButtonHolder, Submit, Div, Button, HTML


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
    data_type = forms.ChoiceField(
        required=True,
        choices=[
                #("fasta", "FASTA"),
                #("fastagene", "FASTA (Gene list)"),
                ("genbank", "GenBank"),
                #("optgene", "OptGene / BioOpt"),
                ("sbml", "System Biology Markup Language (SBML)"),
                ("proopdb", "ProOpDB Operon Prediction")
                ],
        widget=forms.RadioSelect,
        label='File format',
        help_text='Select the format of your file from the list'
        )
    file = forms.FileField(
        required=True,
        widget=forms.ClearableFileInput,
        label='File',
        help_text='Please provide a file in your selected format'
        )
    reason = forms.CharField(
        required=True,
        widget=forms.TextInput,
        label="Summary",
        help_text="Sum up what you did in a few words"
    )
    
    # GenBank specific
    chromosome = forms.CharField(
        required=False,
        widget=forms.TextInput,
        label="Chromosome/Plasmid name",
        help_text="Provide a new name"
    )
    chromosome_wid = forms.SlugField(
        required=False,
        widget=forms.TextInput,
        label="Unique identifier",
        help_text="Identifiers are visible in the URL and may only contain alphanumeric, - and _"
    )

    def __init__(self, *args, **kwargs):
        super(ImportDataForm, self).__init__(*args, **kwargs)
        
        self.helper = FormHelper()
        self.helper.layout = Layout(
            Panel(
                'Select file format',
                'data_type',
            ),
            Panel(
                'Provide additional information for GenBank',
                'chromosome',
                'chromosome_wid',
                css_id="genbank"
            ),
            Panel(
                'Select a file to import',
                'file'
            ),
            Panel(
                'Summary',
                'reason'
            ),
            ButtonHolder(
                Submit('submit', 'Submit')
            )
        )
    
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
        required=True,
        widget=forms.TextInput,
        label="Species name",
        help_text="Provide a new name"
        )
    new_wid = forms.SlugField(
        required=True,
        widget=forms.TextInput,
        label="Unique identifier",
        help_text="Identifiers are visible in the URL and may only contain alphanumeric, - and _"
    )
    reason = forms.CharField(
        required=True,
        widget=forms.TextInput,
        label="Summary",
        help_text="Sum up what you did in a few words"
    )
        
    def __init__(self, *args, **kwargs):
        super(ImportSpeciesForm, self).__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.layout = Layout(
            Panel(
                'Species',
                'new_species',
                'new_wid'
            ),
            Panel(
                'Summary',
                'reason'
            ),
            ButtonHolder(
                Submit('submit', 'Submit')
            )
        )
    
    def clean(self):
        cleaned_data = super(ImportSpeciesForm, self).clean()
        
        new_wid = cleaned_data.get("new_wid")
        
        if new_wid:
            if cmodels.Species.objects.filter(wid = new_wid).exists():
                self._errors["new_wid"] = self.error_class(['The identifier specified is already in use'])
                del cleaned_data["new_wid"]

        return cleaned_data


class DeleteForm(forms.Form):
    reason = forms.CharField(
        required=False,
        widget=forms.TextInput,
        label="Summary",
        help_text="Explain why you deleted this entry"
    )

    def __init__(self, *args, **kwargs):
        super(DeleteForm, self).__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.layout = Layout(
            ModalHeader(
                'Delete {{ object.wid }}',
            ),
            ModalBody(
                HTML("Are you sure you want to delete {{ model_verbose_name|lower }} {{ object.wid }}?"),
                Fieldset(
                    None,
                    'reason'
                ),
            ),
            ModalFooter(
                ButtonHolder(
                    Submit('delete-submit', 'Delete', css_class='btn btn-danger'),
                    CancelButton('delete-cancel', 'Cancel')
                )
            ),
        )

    def clean(self):
        cleaned_data = super(DeleteForm, self).clean()

        return cleaned_data


class CreateBasketForm(forms.Form):
    name = forms.CharField(
        required=True,
        widget=forms.TextInput,
        label="Name",
        help_text="Provide a name for the new basket"
    )

    def __init__(self, *args, **kwargs):
        super(CreateBasketForm, self).__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.layout = Layout(
            ModalHeader(
                'Create New Basket',
            ),
            ModalBody(
                'name'
            ),
            ModalFooter(
                ButtonHolder(
                    Button('basket-create-submit', 'Create Basket', css_class='btn-primary'),
                    CancelButton('delete-cancel', 'Cancel')
                )
            ),
        )


class RenameBasketForm(forms.Form):
    name = forms.CharField(
        required=True,
        widget=forms.TextInput,
        label="Name",
        help_text="Provide a new name for the basket"
    )

    def __init__(self, *args, **kwargs):
        super(RenameBasketForm, self).__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.layout = Layout(
            ModalHeader(
                'Rename Basket',
            ),
            ModalBody(
                'name'
            ),
            ModalFooter(
                ButtonHolder(
                    Button('basket-create-submit', 'Rename Basket', css_class='btn-primary'),
                    CancelButton('delete-cancel', 'Cancel')
                )
            ),
        )