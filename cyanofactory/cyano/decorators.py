import cyano.models as model
from cyano.helpers import get_queryset_object_or_404, render_queryset_to_response_error
from django.utils.functional import wraps
from django.contrib.auth.models import User
from settings import DEBUG

def read_permission_required(function):
    """Checks whether a user has a read_permission for the specified entry
    """
    @wraps(function)
    def wrapper(request, species_wid = None, wid = None, *args, **kw):
        if species_wid != None:
            species = get_queryset_object_or_404(model.Species.objects.filter(wid = species_wid))
        
        user = request.user
        
        if not user.is_authenticated():
            # Special handling for guests
            user = User.objects.get(username = "guest")

        profile = user.profile
        allow, deny = profile.has_permission(species, model.PermissionEnum.READ_NORMAL)        
        if allow and not deny:
            # Has permission
            return function(request, species_wid, *args, **kw)
        else:
            msg = "You don't have permission to access this item"
            extra = None
            
            if DEBUG:
                allow_mask, deny_mask = profile.get_permission_mask(species)
                extra = "Permissions are (allow, deny, needed)<br>"
                extra += "<pre>{0}</pre><pre>{1}</pre><pre>{2}</pre>".format(
                        bin(allow_mask or 0)[2:].rjust(8, '0'),
                        bin(deny_mask or 0)[2:].rjust(8, '0'),
                        bin(model.PermissionEnum.READ_NORMAL)[2:].rjust(8, '0'))
                print extra
                        
            return render_queryset_to_response_error(
                request,
                species = species,
                error = 403,
                msg = msg,
                msg_debug = extra)

    return wrapper
