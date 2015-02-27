"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, ButtonHolder, Button, Fieldset, Submit

from django import forms
from django.core.exceptions import ValidationError
from bioparser.optgene import OptGeneParser
from cyano.layout import ModalHeader, ModalBody, ModalFooter, CancelButton


class UploadModelForm(forms.Form):
    name = forms.CharField(
        required=True,
        widget=forms.TextInput,
        label="Name"
    )
    file = forms.FileField(
        required=True,
        widget=forms.ClearableFileInput,
        label="Model File",
        help_text="Please provide a model in BioOpt format."
    )
    
    def __init__(self, *args, **kwargs):
        super(UploadModelForm, self).__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.form_action = "cyano-design-save-upload-form"
        self.helper.layout = Layout(
            ModalHeader(
                'Create new model',
            ),
            ModalBody(
                'name',
                'file'
            ),
            ModalFooter(
                ButtonHolder(
                    Button('upload-submit', 'Create Model', css_class='btn-primary'),
                    CancelButton('upload-cancel', 'Cancel')
                )
            ),
        )

    def clean_file(self):
        f = self.cleaned_data['file']

        o = OptGeneParser(f)

        if len(o.errors) > 0:
            raise ValidationError(", ".join(o.errors))

