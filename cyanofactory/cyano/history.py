from django.db.models.fields.related import ForeignKey, ReverseSingleRelatedObjectDescriptor,\
    ManyToManyField
from django.utils import six
from django.db.models.related import RelatedObject

class HistoryReverseSingleRelatedObjectDescriptor(ReverseSingleRelatedObjectDescriptor):
    def __get__(self, instance, instance_type=None):
        from cyano.models import Revision, TableMetaColumn
        if instance is None:
            return self
        try:
            rel_obj = getattr(instance, self.cache_name)
        except AttributeError:
            if hasattr(instance, "detail_history"):
                print "redirecting to " + str(instance.detail_history)
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
                                print "revisioning: " + str(rel_obj) + " - " + field.name
                                new_value = field.to_python(history_obj[0].new_value)        
                                setattr(rel_obj, field.name, new_value)

                if not self.field.rel.multiple:
                    setattr(rel_obj, self.field.related.get_cache_name(), instance)
            setattr(instance, self.cache_name, rel_obj)
        if rel_obj is None and not self.field.null:
            raise self.field.rel.to.DoesNotExist
        else:
            return rel_obj

class HistoryForeignKey(ForeignKey):
    def contribute_to_class(self, cls, name):
        super(ForeignKey, self).contribute_to_class(cls, name)
        setattr(cls, self.name, HistoryReverseSingleRelatedObjectDescriptor(self))
        if isinstance(self.rel.to, six.string_types):
            target = self.rel.to
        else:
            target = self.rel.to._meta.db_table
        cls._meta.duplicate_targets[self.column] = (target, "o2m")
