from __future__ import absolute_import

# This will make sure the app is always imported when
# Django starts so that shared_task will use this app.
from .celery import app as celery_app

def patch_guardian():
    # Make guardian 1.8 compatible until an update is out
    from django.contrib.contenttypes.models import ContentType
    from django.db.models.base import Model
    import guardian.utils
    import django

    def get_obj_perms_model_patched(obj, base_cls, generic_cls):
        if isinstance(obj, Model):
            obj = obj.__class__
        ctype = ContentType.objects.get_for_model(obj)
        for attr in obj._meta.get_all_related_objects():
            if django.VERSION < (1, 8):
                model = getattr(attr, 'model', None)
            else:
                model = getattr(attr, 'related_model', None)
            if (model and issubclass(model, base_cls) and
                    model is not generic_cls):
                # if model is generic one it would be returned anyway
                if not model.objects.is_generic():
                    # make sure that content_object's content_type is same as
                    # the one of given obj
                    fk = model._meta.get_field_by_name('content_object')[0]
                    if ctype == ContentType.objects.get_for_model(fk.rel.to):
                        return model
        return generic_cls

    guardian.utils.get_obj_perms_model = get_obj_perms_model_patched

patch_guardian()
