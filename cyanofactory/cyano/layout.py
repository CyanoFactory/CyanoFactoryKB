from crispy_forms.layout import Fieldset, Div, Button


class Panel(Fieldset):
    """
    Layout object. It wraps fields in a <div class="panel panel-default">
    """
    template = "cyano/layout/panel.html"

    def __init__(self, legend, *fields, **kwargs):
        super(Panel, self).__init__(legend, *fields, **kwargs)


class Modal(Fieldset):
    """
    Layout object. It wraps fields in a <div class="panel panel-default">
    """
    template = "cyano/layout/modal.html"

    def __init__(self, legend, *fields, **kwargs):
        super(Modal, self).__init__(legend, *fields, **kwargs)


class ModalHeader(Fieldset):
    """
    Layout object. It wraps fields in a <div class="modal-body">
    """
    template = "cyano/layout/modal_header.html"

    def __init__(self, *fields, **kwargs):
        kwargs["css_class"] = "modal-header"
        super(ModalHeader, self).__init__(*fields, **kwargs)


class ModalBody(Div):
    """
    Layout object. It wraps fields in a <div class="modal-body">
    """
    def __init__(self, *fields, **kwargs):
        kwargs["css_class"] = "modal-body"
        super(ModalBody, self).__init__(*fields, **kwargs)


class ModalFooter(Div):
    """
    Layout object. It wraps fields in a <div class="modal-footer">
    """
    def __init__(self, *fields, **kwargs):
        kwargs["css_class"] = "modal-footer"
        super(ModalFooter, self).__init__(*fields, **kwargs)


class CancelButton(Button):
    """
    Used to create a Bootstrap Cancel input descriptor for the {% crispy %} template tag::

        button = CancelButton('Button 1', 'Press Me!')

    .. note:: The first argument is also slugified and turned into the id for the button.
    """
    template = "cyano/layout/cancel_button.html"