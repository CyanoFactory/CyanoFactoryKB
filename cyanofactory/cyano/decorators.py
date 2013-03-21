import cyano.models as models
from cyano.helpers import get_queryset_object_or_404, render_queryset_to_response_error
from django.utils.functional import wraps
from django.contrib.auth.models import User
from django.db.models.loading import get_model
import django.http as http
from settings import DEBUG

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
            species = get_queryset_object_or_404(models.Species.objects.filter(wid = species_wid))
        if model_type:
            model = get_model("cyano", model_type)
            if model == None:
                raise http.Http404
        if model and wid:
            item = get_queryset_object_or_404(model.objects.filter(wid = wid))

        # Prepare keyword arguments
        __assign_if_not_false(kw, "species", species)
        __assign_if_not_false(kw, "model", model)
        __assign_if_not_false(kw, "item", item)
        
        return function(request, *args, **kw)

    return wrapper

def permission_required(permission):
    """Checks whether a user has a read_permission for the specified entry
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
    
            profile = user.profile
            allow, deny = profile.has_permission(species, permission)        
            if allow and not deny:
                # Has permission
                return function(request, *args, **kw)
            else:
                msg = "You don't have permission to access this item"
                extra = None
                
                if DEBUG:
                    allow_mask, deny_mask = profile.get_permission_mask(species)
                    extra = "DEBUG: Permissions are (allow, deny, needed):<br>"
                    extra += "<pre>{0}</pre><pre>{1}</pre><pre>{2}</pre>".format(
                            bin(allow_mask or 0)[2:].rjust(8, '0'),
                            bin(deny_mask or 0)[2:].rjust(8, '0'),
                            bin(permission)[2:].rjust(8, '0'))
                    print extra
                            
                return render_queryset_to_response_error(
                    request,
                    species = species,
                    error = 403,
                    msg = msg,
                    msg_debug = extra)
        return wrapper
    return decorator
