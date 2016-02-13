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
from django.contrib.auth import authenticate, login
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist

import cyano.helpers as chelpers
from cyano.layout import Panel, Modal, ModalBody, ModalFooter, ModalHeader, CancelButton
import cyano.models as cmodels
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Fieldset, Submit, Button, HTML
from crispy_forms.bootstrap import FormActions
from django.contrib.auth.forms import AuthenticationForm

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
            FormActions(
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
            FormActions(
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
                FormActions(
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
                FormActions(
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
                FormActions(
                    Button('basket-create-submit', 'Rename Basket', css_class='btn-primary'),
                    CancelButton('delete-cancel', 'Cancel')
                )
            ),
        )


class LoginForm(AuthenticationForm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.fields['username'].label = "Username or Mail"

    def clean(self):
        username = self.cleaned_data.get('username')
        password = self.cleaned_data.get('password')

        if username and password:
            if "@" in username:
                try:
                    u = User.objects.get(email=username)
                    username = u.username
                except ObjectDoesNotExist:
                    pass

            self.user_cache = authenticate(username=username,
                                            password=password)

            if self.user_cache is None:
                raise forms.ValidationError(
                    self.error_messages['invalid_login'],
                    code='invalid_login',
                    params={'username': self.username_field.verbose_name},
                )
            else:
                self.confirm_login_allowed(self.user_cache)

        return self.cleaned_data


class ChangeProfileForm(forms.Form):
    username = forms.SlugField(
        required=True,
        widget=forms.TextInput,
        label="Username",
        max_length=20
    )

    first_name = forms.CharField(
        widget=forms.TextInput,
        required=False,
        label="First name",
        max_length=254
    )

    last_name = forms.CharField(
        widget=forms.TextInput,
        required=False,
        label="Last name",
        max_length=254
    )

    email = forms.EmailField(
        required=True,
        widget=forms.EmailInput,
        label="E-Mail address",
        max_length=254
    )

    affiliation = forms.CharField(
        required=False,
        widget=forms.TextInput,
        label="Affiliation",
        max_length=254
    )

    password = forms.CharField(
        label="Your Password",
        widget=forms.PasswordInput,
        help_text="Enter your current password to verify your changes",
        max_length=254
    )

    def __init__(self, request=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.request = request
        user = self.request.user

        self.fields["username"].initial = user.username
        self.fields["first_name"].initial = user.first_name
        self.fields["last_name"].initial = user.last_name
        self.fields["email"].initial = user.email
        self.fields["affiliation"].initial = user.profile.affiliation

    def clean(self):
        password = self.cleaned_data.get('password')

        if password:
            auth = authenticate(username=self.request.user.username,
                                password=password)

            if auth is None:
                raise forms.ValidationError(
                    "Wrong password",
                    code='invalid_login',
                    params={'password': self.fields["password"].label},
                )

        username = self.cleaned_data.get('username')

        if username:
            if username != self.request.user.username:
                if User.objects.filter(username=username).exists():
                    raise forms.ValidationError(
                        "Username already in use",
                        params={'username': "username"})

        email = self.cleaned_data.get("email")
        if email:
            if email != self.request.user.email:
                if User.objects.filter(email=email).exists():
                    raise forms.ValidationError(
                        "E-Mail address already in use",
                        params={'email': "email"})

        return self.cleaned_data

    def save(self):
        u = self.request.user

        u.username = self.cleaned_data.get('username')
        u.email = self.cleaned_data.get('email')
        u.first_name = self.cleaned_data.get('first_name')
        u.last_name = self.cleaned_data.get('last_name')
        u.profile.affiliation = self.cleaned_data.get('affiliation')

        u.profile.save()
        u.save()


class CreateProfileForm(forms.Form):
    username = forms.SlugField(
        required=True,
        widget=forms.TextInput,
        label="Username",
        max_length=20
    )

    first_name = forms.CharField(
        widget=forms.TextInput,
        required=False,
        label="First name",
        max_length=254
    )

    last_name = forms.CharField(
        widget=forms.TextInput,
        required=False,
        label="Last name",
        max_length=254
    )

    email = forms.EmailField(
        required=True,
        widget=forms.EmailInput,
        label="E-Mail address",
        max_length=254
    )

    email2 = forms.EmailField(
        required=True,
        widget=forms.EmailInput,
        label="Confirm E-Mail address",
        max_length=254
    )

    affiliation = forms.CharField(
        required=False,
        widget=forms.TextInput,
        label="Affiliation",
        max_length=254
    )

    password = forms.CharField(
        label="Password",
        widget=forms.PasswordInput,
        max_length=254
    )

    password2 = forms.CharField(
        label="Confirm Password",
        widget=forms.PasswordInput,
        max_length=254
    )

    captcha = forms.CharField(
        label="Challenge",
        widget=forms.TextInput,
        required=True,
        help_text="Type the reverse complement of ATCGAAGG"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.layout = Layout(
            Panel(
                'Account information',
                'username',
                'email',
                'email2',
                'password',
                'password2'
            ),
            Panel(
                'Account details',
                'first_name',
                'last_name',
                'affiliation'
            ),
            Panel(
                'Challenge',
                'captcha'
            ),
            FormActions(
                Submit('submit', 'Register')
            )
        )

    def clean(self):
        username = self.cleaned_data.get("username")
        if username:
            if User.objects.filter(username=username).exists():
                raise forms.ValidationError(
                    "Username already in use",
                    params={'username': "Username"})

        password = self.cleaned_data.get("password")
        password2 = self.cleaned_data.get("password2")

        if password and password2:
            if password != password2:
                raise forms.ValidationError(
                    "The passwords mismatch",
                    params={'password': "Password"})

        email = self.cleaned_data.get("email")
        email2 = self.cleaned_data.get("email2")

        if email and email2:
            if email != email2:
                raise forms.ValidationError(
                    "The mail addresses mismatch",
                    params={'email': "email"})

            if User.objects.filter(email=email).exists():
                raise forms.ValidationError(
                    "E-Mail address already in use",
                    params={'email': "email"})

        captcha = self.cleaned_data.get("captcha")

        if captcha:
            if captcha.strip().lower() != "ccttcgat":
                raise forms.ValidationError(
                    "You failed in biology :(",
                    params={'captcha': "captcha"})

        return self.cleaned_data

    def save(self, request):
        user = User.objects.create_user(
            username=self.cleaned_data["username"],
            email=self.cleaned_data["email"],
            password=self.cleaned_data["password"]
        )

        if self.cleaned_data.get("first_name"):
            user.first_name = self.cleaned_data["first_name"]

        if self.cleaned_data.get("last_name"):
            user.last_name = self.cleaned_data["last_name"]

        if self.cleaned_data.get("affiliation"):
            user.profile.affiliation = self.cleaned_data["affiliation"]
            user.profile.save()

        user.save()

        login_user = authenticate(username=user.username, password=self.cleaned_data["password"])
        login(request, login_user)