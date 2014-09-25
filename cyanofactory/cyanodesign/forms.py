"""
Copyright (c) 2014 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from django import forms


class UploadModelForm(forms.Form):
    name = forms.CharField(
        required=True,
        widget=forms.TextInput,
        label="name",
        help_text="Name"
    )
    file = forms.FileField(
        required=True,
        widget=forms.ClearableFileInput,
        label='Model to upload'
    )
