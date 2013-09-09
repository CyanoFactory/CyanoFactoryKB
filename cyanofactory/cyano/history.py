from django.db.models.fields.related import ForeignKey, ReverseSingleRelatedObjectDescriptor,\
    ManyToManyField, ReverseManyRelatedObjectsDescriptor
from django.db.models.related import RelatedObject
from django.db import router
from django.db.models.manager import Manager

class HistoryReverseSingleRelatedObjectDescriptor(ReverseSingleRelatedObjectDescriptor):
    def __get__(self, instance, instance_type=None):
        from cyano.models import Revision, TableMetaColumn
        if instance is None:
            return self
        try:
            rel_obj = getattr(instance, self.cache_name)
        except AttributeError:
            #if hasattr(instance, "detail_history"):
            #    print "redirecting to " + str(instance.detail_history)
            val = getattr(instance, self.field.attname)
            if val is None:
                rel_obj = None
            else:
                other_field = self.field.rel.get_related_field()
                if other_field.rel:
                    params = {'%s__%s' % (self.field.rel.field_name, other_field.rel.field_name): val}
                else:
                    params = {'%s__exact' % self.field.rel.field_name: val}
                qs = self.get_query_set(instance=instance)
                # Assuming the database enforces foreign keys, this won't fail.
                rel_obj = qs.get(**params)
                if hasattr(instance, "detail_history"):
                    rel_obj.detail_history = instance.detail_history
                    for field in rel_obj._meta.fields:
                        if isinstance(field, RelatedObject):
                            pass
                        elif isinstance(field, ManyToManyField):
                            pass
                        elif isinstance(field, ForeignKey):
                            pass
                        else:
                            column = TableMetaColumn.get_by_field(field)
                            history_obj = Revision.objects.filter(current = rel_obj, detail__lte = instance.detail_history, column = column).order_by("-detail")
                            if history_obj.exists():
                                #print "revisioning: " + str(rel_obj) + " - " + field.name
                                new_value = field.to_python(history_obj[0].new_value)        
                                setattr(rel_obj, field.name, new_value)

                if not self.field.rel.multiple:
                    setattr(rel_obj, self.field.related.get_cache_name(), instance)
            setattr(instance, self.cache_name, rel_obj)
        if rel_obj is None and not self.field.null:
            raise self.field.rel.to.DoesNotExist
        else:
            return rel_obj

"""Creates a manager that subclasses 'superclass' (which is a Manager)
and adds behavior for many-to-many related objects."""

def create_history_many_related_manager(superclass):
    class HistoryManyRelatedManager(superclass):             
        def get_query_set(self):
            from cyano.models import RevisionManyToMany, TableMetaManyToMany
            if hasattr(self.instance, "detail_history"):
                #print "m2m history"

                try:
                    return self.instance._prefetched_objects_cache[self.prefetch_cache_name]
                except (AttributeError, KeyError):
                    db = self._db or router.db_for_read(self.instance.__class__, instance=self.instance)
                    
                    table = TableMetaManyToMany.get_by_m2m_model_name(self.through._meta.object_name)
                    
                    vals = RevisionManyToMany.objects.filter(current = self.instance.pk, detail__lte = self.instance.detail_history, table = table).values_list('new_value', flat = True)

                    queryset = Manager.get_query_set(self).using(db)._next_is_sticky().filter(**self.core_filters)
                    queryset = queryset.filter(pk__in = vals)
                    
                    #print queryset, self.core_filters, self.target_field_name, self.source_field_name, self._fk_val
                    
                    return queryset
            else:
                queryset = super(HistoryManyRelatedManager, self).get_query_set()
            
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

