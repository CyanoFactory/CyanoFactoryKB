class WarehouseRouter(object):
    def db_for_read(self, model, **hints):
        #if model._meta.app_label == 'biowarehouse':
        #    return 'cyano'
        
        #if model._meta.app_label == 'djcelery':
        #    return "djcelery"

        if model._meta.app_label == 'cyanointeraction':
        #    return 'stringdb'
            if model._meta.db_interaction_div == 'string' :
                return 'stringdb'
            if model._meta.db_interaction_div == 'stitch':
                return 'stitchdb'
        return None

    def db_for_write(self, model, **hints):
        #if model._meta.app_label == 'biowarehouse':
        #    return 'cyano'
        if model._meta.app_label == 'cyanointeraction':
            return 'stringdb'

        if model._meta.app_label == 'djcelery':
            return "djcelery"
        
        return None

    def allow_relation(self, obj1, obj2, **hints):
        """
        Allow relations if a model in the auth app is involved.
        """
        #if obj1._meta.app_label == 'biowarehouse' or \
        #   obj2._meta.app_label == 'biowarehouse':
        #    return True
        if obj1._meta.app_label == 'cyanointeraction' or \
           obj2._meta.app_label == 'cyanointeraction':
            return True


        return None

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        """
        Make sure the auth app only appears in the 'auth_db'
        database.
        """
        if db == 'cyano':
            return app_label == 'biowarehouse'
        elif app_label == 'biowarehouse':
            return False

        if db == 'stringdb':
            return app_label == 'cyanointeraction'
        elif app_label == 'cyanointeraction':
            return False

        return None
