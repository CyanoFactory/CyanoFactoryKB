from django import template

register = template.Library()

@register.filter
def index(value, arg):
    """Access an array (or hash) with the given value"""
    try:
        return value[arg]
    except KeyError:
        from django.conf import settings

        value = settings.TEMPLATE_STRING_IF_INVALID
