class WarehouseRouter(object):
    def db_for_read(self, model, **hints):
        if model._meta.app_label == 'biowarehouse':
            return 'cyano'
        return 'default'

    def db_for_write(self, model, **hints):
        if model._meta.app_label == 'biowarehouse':
            return 'biowarehouse'
        return 'default'

    def allow_relation(self, obj1, obj2, **hints):
        """
        Allow relations if a model in the auth app is involved.
        """
        if obj1._meta.app_label == 'biowarehouse' or \
           obj2._meta.app_label == 'biowarehouse':
           return True
        return None

    def allow_syncdb(self, db, model):
        """
        Make sure the auth app only appears in the 'auth_db'
        database.
        """
        if db == 'cyano':
            return model._meta.app_label == 'biowarehouse'
        elif model._meta.app_label == 'biowarehouse':
            return False
        return None
