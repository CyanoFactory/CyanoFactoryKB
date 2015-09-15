"""
Copyright (c) 2013 Gabriel Kind <gkind@hs-mittweida.de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""
from django.contrib.contenttypes.models import ContentType

from django.db.models.fields.related import ForeignKey, ReverseSingleRelatedObjectDescriptor,\
    ManyToManyField, ReverseManyRelatedObjectsDescriptor, ManyToOneRel
from django.db import router
from django.db.models.manager import Manager

import json

class HistoryReverseSingleRelatedObjectDescriptor(ReverseSingleRelatedObjectDescriptor):
    def __get__(self, instance, instance_type=None):
        from cyano.models import Revision
        if instance is None:
            return self
        try:
            rel_obj = getattr(instance, self.cache_name)
        except AttributeError:
            val = self.field.get_local_related_value(instance)
            if None in val:
                rel_obj = None
            else:
                params = {
                    rh_field.attname: getattr(instance, lh_field.attname)
                    for lh_field, rh_field in self.field.related_fields}
                qs = self.get_queryset(instance=instance)
                extra_filter = self.field.get_extra_descriptor_filter(instance)
                if isinstance(extra_filter, dict):
                    params.update(extra_filter)
                    qs = qs.filter(**params)
                else:
                    qs = qs.filter(extra_filter, **params)
                # Assuming the database enforces foreign keys, this won't fail.
                rel_obj = qs.get()

                if hasattr(instance, "detail_history"):

                    ct_id = ContentType.objects.get_for_model(instance._meta.concrete_model).pk

                    revisions = Revision.objects.filter(
                        object_id=instance.pk,
                        content_type__pk=ct_id,
                        detail_id__lte=instance.detail_history).order_by("-detail")

                    data = {}
                    for revision in revisions:
                        new_data = json.loads(revision.new_data)
                        data.update(new_data)

                    ##print "FK", self.field.name
                    rel_obj.detail_history = instance.detail_history
                    for field in rel_obj._meta.fields:
                        if isinstance(field, ManyToOneRel):
                            pass
                        elif isinstance(field, ManyToManyField):
                            pass
                        elif isinstance(field, ForeignKey):
                            pass
                        else:
                            #print "revisioning: " + str(rel_obj) + " - " + field.name
                            if field.name in data:
                                new_value = field.to_python(data[field.name])
                                setattr(rel_obj, field.name, new_value)
                            else:
                                ##print "Not in", field.name
                                pass

                if not self.field.rel.multiple:
                    setattr(rel_obj, self.field.rel.get_cache_name(), instance)
            setattr(instance, self.cache_name, rel_obj)
        if rel_obj is None and not self.field.null:
            raise self.RelatedObjectDoesNotExist(
                "%s has no %s." % (self.field.model.__name__, self.field.name)
            )
        else:
            return rel_obj

"""Creates a manager that subclasses 'superclass' (which is a Manager)
and adds behavior for many-to-many related objects."""

def create_history_many_related_manager(superclass):
    class HistoryManyRelatedManager(superclass):             
        def get_queryset(self):
            if hasattr(self.instance, "detail_history"):
                from cyano.models import Revision

                #print "m2m history"

                try:
                    return self.instance._prefetched_objects_cache[self.prefetch_cache_name]
                except (AttributeError, KeyError):
                    db = self._db or router.db_for_read(self.instance.__class__, instance=self.instance)

                    ct_id = ContentType.objects.get_for_model(self.model).pk

                    queryset = Manager.get_queryset(self).using(db)._next_is_sticky().all()#.filter(**self.core_filters)

                    revision = Revision.objects.filter(
                        object_id__in=queryset,
                        content_type__pk=ct_id,
                        detail_id__lte=self.instance.detail_history).order_by("-detail").first()

                    data = []

                    if revision is not None:
                        new_data = json.loads(revision.new_data)
                        if self.prefetch_cache_name in new_data:
                            data = new_data[self.prefetch_cache_name]

                    ##print "History:", self.instance.detail_history, "Data:", data

                    queryset = queryset.filter(pk__in=data)

                    for item in queryset:
                        setattr(item, "detail_history", self.instance.detail_history)

                    #queryset = queryset.filter(pk__in = vals)
                    
                    ##print "Queryset:", queryset#, self.core_filters, self.target_field_name, self.source_field_name, self._fk_val
                    
                    return queryset
            else:
                queryset = super(HistoryManyRelatedManager, self).get_queryset()

            #print queryset
            return queryset

    return HistoryManyRelatedManager

class HistoryReverseManyRelatedObjectsDescriptor(ReverseManyRelatedObjectsDescriptor):
    def __get__(self, instance, instance_type=None):        
        manager = super(HistoryReverseManyRelatedObjectsDescriptor, self).__get__(instance, instance_type=None)
        
        manager.__class__ = create_history_many_related_manager(manager.__class__)

        return manager

class HistoryForeignKey(ForeignKey):
    def contribute_to_class(self, cls, name):
        super(HistoryForeignKey, self).contribute_to_class(cls, name)
        setattr(cls, self.name, HistoryReverseSingleRelatedObjectDescriptor(self))

class HistoryManyToManyField(ManyToManyField):
    def contribute_to_class(self, cls, name):
        super(HistoryManyToManyField, self).contribute_to_class(cls, name)
        setattr(cls, self.name, HistoryReverseManyRelatedObjectsDescriptor(self))
