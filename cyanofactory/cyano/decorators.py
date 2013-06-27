import cyano.models as models
from cyano.helpers import render_queryset_to_response_error
from django.utils.functional import wraps
from django.contrib.auth.models import User
from django.db.models.loading import get_model
from settings import DEBUG
from django.http.response import Http404
from django.core.exceptions import ObjectDoesNotExist, MultipleObjectsReturned

def __get_and_delete(kw, key):
    if not key in kw:
        return False
    
    res = kw[key]
    del kw[key]
    return res
    
def __assign_if_not_false(kw, key, obj):
    if obj != False:
        kw[key] = obj

def resolve_to_objects(function):
    @wraps(function)
    def wrapper(request, *args, **kw):
        # False instead of None, None is default value for passed but not set
        species_wid = __get_and_delete(kw, "species_wid")
        model_type = __get_and_delete(kw, "model_type")
        wid = __get_and_delete(kw, "wid")

        species = False
        model = False
        item = False

        if species_wid:
            species = models.Species.objects.filter(wid = species_wid)
            try:
                species = species.get()
            except (ObjectDoesNotExist, MultipleObjectsReturned):
                return render_queryset_to_response_error(
                    request,
                    error = 404,
                    msg = "The requested species \"{}\" was not found.".format(species_wid))
        if model_type:
            model = get_model("cyano", model_type)
            if model == None:
                return render_queryset_to_response_error(
                    request,
                    error = 404,
                    msg = "The requested species \"{}\" has no model \"{}\".".format(species_wid, model_type))
        if model and wid:
            item = model.objects.filter(wid = wid, species = species)
            try:
                item = item.get()
            except (ObjectDoesNotExist, MultipleObjectsReturned):
                return render_queryset_to_response_error(
                    request,
                    error = 404,
                    msg = "The requested species \"{}\" has no item \"{}\".".format(species_wid, wid))

        # Prepare keyword arguments
        __assign_if_not_false(kw, "species", species)
        __assign_if_not_false(kw, "model", model)
        __assign_if_not_false(kw, "item", item)
        
        return function(request, *args, **kw)

    return wrapper

def permission_required(permission):
    """Checks whether a user has a permission for the specified entry
    """
    def decorator(function):
        @wraps(function)
        def wrapper(request, *args, **kw):            
            species = kw.get("species", False)
            model = kw.get("model", False)
            item = kw.get("item", False)
            
            user = request.user
            
            if not user.is_authenticated():
                # Special handling for guests
                user = User.objects.get(username = "guest")
    
            perm = models.Permission.objects.get(name = permission)
            
            profile = user.profile
            
            # Admins always have full access
            if profile.is_admin():
                return function(request, *args, **kw)

            allow_item = None
            deny_item = None

            if item:
                allow_item, deny_item = profile.has_permission(item, perm)
            
            # allow or deny are not None if an item had permission assigned
            # Mix item permissions with species permissions now
            allow_species, deny_species = profile.has_permission(species, perm)
            
            if (allow_species or allow_item) and not (deny_species or deny_item):
                # Has permission
                return function(request, *args, **kw)
            else:
                msg = "You don't have permission to access this item"
                extra = None
                
                if DEBUG:
                    #allow_mask, deny_mask = profile.get_permission_mask(species)
                    extra = "DEBUG: Permissions are (allow, deny, needed):<br>"
                    
                    allow_perms, deny_perms = profile.get_permissions(species)
                    if not allow_perms or not deny_perms:
                        allow_perms = []
                        deny_perms = []
                    
                    if item:
                        allow_perms_item, deny_perms_item = profile.get_permissions(item)
                        if allow_perms_item and deny_perms_item:
                            allow_perms += allow_perms_item
                            deny_perms += deny_perms_item

                    perm_list = [0 for x in range(8)]
                    perm_allow = list(perm_list)
                    perm_deny = list(perm_list)
                    perm_needed = list(perm_list)
                    
                    for i in range(8):
                        cur_perm = models.Permission.get_by_pk(i + 1)
                        perm_allow[i] = 1 if cur_perm in allow_perms else 0
                        perm_deny[i] = 1 if cur_perm in deny_perms else 0
                        perm_needed[i] = 1 if cur_perm == perm else 0
   
                    extra += "<pre>{0}</pre><pre>{1}</pre><pre>{2}</pre>".format(
                            "".join(str(x) for x in perm_allow),
                            "".join(str(x) for x in perm_deny),
                            "".join(str(x) for x in perm_needed))
                    #print extra
                            
                return render_queryset_to_response_error(
                    request,
                    species = species,
                    error = 403,
                    msg = msg,
                    msg_debug = extra)
        return wrapper
    return decorator
